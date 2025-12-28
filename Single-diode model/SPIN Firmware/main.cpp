/*
 * ====================================================================
 * PV EMULATOR APPLICATION
 * * This code implements a real-time PV (Photovoltaic) panel emulator 
 * using a buck converter.
 *
 * It operates in three main modes:
 * 1. IDLE: PWM is disabled.
 * 2. POWER: Acts as a standard buck converter, regulating to a fixed
 * output voltage (voltageReferenceP).
 * 3. EMULATOR: Simulates a PV panel. It finds the intersection of the
 * panel's I-V curve and the attached load's resistance line, then
 * uses a PID controller to regulate the output to that
 * intersection voltage.
 *
 * The core logic involves two main parts:
 * 1. Setup-Time Solver (64-bit double):
 * - Runs once at startup (in setup_routine).
 * - Solves a complex 5-equation system using the Levenberg-Marquardt
 * algorithm to find the 5 parameters (Iph, Is0, A, Rs, Rp) of
 * the single-diode model.
 * - This calculation is done in 64-bit (double) for high precision.
 * * 2. Run-Time Emulator (32-bit float):
 * - Uses the 5 parameters (cast to 32-bit float) to run the emulator.
 * - Pre-calculates a piecewise linear model of the panel's I-V curve.
 * - In real-time, it finds the intersection of the load line and
 * the I-V curve to determine the operating point (V*, I*).
 * - A PID controller in the critical task loop adjusts the PWM
 * duty cycle to force the buck converter's output to match V*.
 * ====================================================================
 */

#include "ShieldAPI.h"
#include "TaskAPI.h"
#include "SpinAPI.h"
#include "pid.h"
#include "zephyr/console/console.h"
#include <math.h> 
#include <stdlib.h> 
#include <time.h> 
#include <float.h> // For DBL_MAX
#include <stdbool.h>
#include <string.h> 

//====================================================================
// Global Declarations
//====================================================================

// --- EMULATOR BEHAVIOR CONSTANTS ---
#define HISTORY_SIZE 2                // Size of history buffer for checking stability
#define N_POINTS 26                   // Number of points to pre-calculate for the I-V curve
#define LOAD_RESISTANCE_THRESHOLD 0.02f // (2%) Relative change in R to trigger recalibration
#define DUTY_CYCLE_CHANGE_THRESHOLD 0.05f // (5%) Change in duty to detect manual override
#define CURRENT_STABILITY_THRESHOLD 0.1f  // (10%) Allowed current variance to be "stable"
#define REQUIRED_PERSISTENCE_COUNT 2  // How many cycles a load change must be detected

// Structure for a single linear segment of the I-V curve (I = a*V + b)
struct Segment {
    float32_t a; // Slope
    float32_t b; // Y-intercept
    float32_t Vmin; // Start voltage of the segment
    float32_t Vmax; // End voltage of the segment
};

// --- Task Prototypes ---
void loop_application_task(); 
void loop_communication_task(); 
void loop_critical_task(); 
void setup_routine(); 

// --- PV MODULE REFERENCE PARAMETERS (32-bit float) ---
// These define the panel to be modeled (e.g., CS6P-250P)
const float32_t ns = 60; // Number of series cells
const float32_t np = 1;  // Number of parallel branches

// Parameters at Standard Test Conditions (STC)
const float32_t Vmp_mod_ref = 30.1f; // Voltage at max power (V)
const float32_t Imp_mod_ref = 8.30f; // Current at max power (A)
const float32_t Voc_mod_ref = 37.2f; // Open-circuit voltage (V)
const float32_t Isc_mod_ref = 8.87f; // Short-circuit current (A)

const float32_t Tref = 25.0f + 273.15f; // Reference Temperature (K)
const float32_t Sref = 1000.0f; // Reference Irradiance (W/m^2)

const float32_t alpha = 0.00065f; // Temperature coefficient of Isc (%/K)
const float32_t beta  = -0.0037f; // Temperature coefficient of Voc (%/K)

// --- CURRENT OPERATING CONDITIONS (32-bit float) ---
static float32_t T = 44.5f + 273.15f; // Current module temperature (K)
static float32_t S = 765.0f;          // Current irradiance (W/m^2)

// Initial estimation of Voc under the new temperature (this value of Voc is simply an approximation)
// Formula: Voc(T) = Voc_ref * (1 + beta * (T - Tref))
static float32_t Voc_estimation = Voc_mod_ref * (1.0f + beta * (T - Tref));

// Duty cycles and Voltage references
float32_t dutyCycle = 0.0f; // Current PWM duty cycle
static float32_t dutyCycleP = 0.8f; // Stored duty cycle for Power mode
static float32_t voltageReferenceP = 0.8f * 50.0f; // Target voltage for Power mode
static float32_t voltageReferenceE_test_point = 1.05*Voc_estimation; // Initial target for Emulator mode (the "Test Point")
static float32_t voltageReferenceE = 0.0f; // Calculated target for Emulator mode (the "Intersection")

// --- 5-PARAMETER MODEL (SINGLE DIODE MODEL) ---
// These globals will be populated by the 64-bit solver at startup
static float32_t Iph_ref = 0.0f; // Photocurrent (A) at STC
static float32_t Is0_ref = 0.0f; // Diode saturation current (A) at STC
static float32_t A       = 0.0f; // Diode ideality factor
static float32_t Rs      = 0.0f; // Series resistance (Ohms)
static float32_t Rp      = 0.0f; // Parallel resistance (Ohms)

// --- PHYSICAL CONSTANTS (32-bit float for runtime) ---
const float32_t q = 1.60217662e-19f; // Elementary charge (C)
const float32_t k = 1.38064852e-23f; // Boltzmann constant (J/K)

// --- MATERIAL CONSTANTS (Silicon) ---
const float32_t E_G0 = 1.166f; // Band gap energy at 0K (eV)
const float32_t k1 = 4.73e-4f;   // Coefficient k1 (eV/K)
const float32_t k2 = 636.0f;   // Coefficient k2 (K)

// --- SOLVER CONSTANTS (64-bit double for high precision setup) ---
const double SOLVER_q_d = 1.60217662e-19; 
const double SOLVER_k_d = 1.38064852e-23; 
// Cell-level parameters converted to double for the solver
const double SOLVER_ISC_MOD_REF_d = (double)Isc_mod_ref;
const double SOLVER_VOC_MOD_REF_d = (double)Voc_mod_ref;
const double SOLVER_VMP_MOD_REF_d = (double)Vmp_mod_ref;
const double SOLVER_IMP_MOD_REF_d = (double)Imp_mod_ref;
const double SOLVER_TREF_K_d = (double)Tref;

// --- SOLVER NUMERICAL TOLERANCES (64-bit double) ---
const double JACOBIAN_EPSILON_d = 1e-8; // Step size for finite difference
const double MIN_ABS_STEP_d = 1e-12; // Minimum step
const double CONVERGENCE_TOL_d = 1e-10; // Target for residual norm
const double STEP_TOL_d = 1e-10; // Target for step size

// --- CONTROL & STATE VARIABLES (32-bit float) ---

// Operating modes
enum modeMenu { MODE_IDLE = 0, MODE_POWER, MODE_EMULATOR };
uint8_t mode = MODE_IDLE; // Current operating mode

// Measured values from sensors
static float32_t lowVoltage1, lowVoltage2; 
static float32_t lowCurrent1, lowCurrent2; 
static float32_t highCurrent, highVoltage; 

static uint32_t controlTaskPeriod = 100; // Period of the critical task (us)
static bool pwmEnable = false; // Flag to track PWM state
uint8_t receivedChar; // For serial communication

// --- PID Controller ---
static float32_t kp = 0.000215f; 
static float32_t Ti = 7.5175e-5f; 
static float32_t Td = 0.0f; 
static float32_t Np = 0.0f; 
static float32_t upperBound = 1.0f; 
static float32_t lowerBound = 0.0f; 
static float32_t Ts = controlTaskPeriod * 1e-6f; // PID sample time (s)
static PidParams pidParams(Ts, kp, Ti, Td, Np, lowerBound, upperBound);
static Pid pid;

// --- Measurement Accumulators ---
// Used to average measurements in the application task
static volatile float32_t sumLowCurrent1 = 0.0f, sumLowVoltage1 = 0.0f;
static volatile float32_t sumLowCurrent2 = 0.0f, sumLowVoltage2 = 0.0f;
static volatile float32_t sumHighCurrent = 0.0f, sumHighVoltage = 0.0f;
static volatile uint32_t measurementCount = 0;

// --- Emulator State ---
static float32_t currentHistory[HISTORY_SIZE] = {0.0f, 0.0f};
static bool emulatorSteadyState = false; // True after first intersection is found
static bool waitingSteadyState = false; // True if waiting for system to stabilize

Segment segments[N_POINTS - 1]; // Array to hold the linear I-V curve segments

static float32_t lastSteadyVoltage = 0.0f;
static float32_t lastSteadyCurrent = 0.0f;
static float32_t lastDutyCycle = 0.0f;

static uint32_t loadChangeCounter = 0; // Counter to detect persistent load changes


//====================================================================
// --- START: (5-Parameter Solver) Implementations (64-bit Double) ---
// These functions run ONCE at startup to find the 5 model parameters.
//====================================================================

// --- Solver-local Typedefs (double) ---
typedef double Vector5[5];
typedef double Matrix5x5[5][5];

// --- Utility Functions (double) ---

bool is_finite_double(double val) {
    return isfinite(val); 
}

// Caps the argument of exp() to prevent overflow
double cap_exp_arg_double(double arg) {
    const double cap = 700.0; 
    if (arg > cap) return cap;
    if (arg < -cap) return -cap;
    return arg;
}

// Matrix transpose (5x5)
void transposeMatrix_d(const Matrix5x5 A, Matrix5x5 B, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            B[j][i] = A[i][j];
        }
    }
}

// Matrix multiplication (5x5)
void matrixMultiply_d(const Matrix5x5 A, const Matrix5x5 B, Matrix5x5 C, int r1, int c1, int c2) {
    for (int i = 0; i < r1; ++i) {
        for (int j = 0; j < c2; ++j) {
            double sum = 0.0;
            for (int k = 0; k < c1; ++k) {
                sum += A[i][k] * B[k][j];
            }
            C[i][j] = sum;
        }
    }
}

// Matrix-vector multiplication (5x1)
void matrixVectorMultiply_d(const Matrix5x5 A, const Vector5 v, Vector5 result, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        double sum = 0.0;
        for (int j = 0; j < cols; ++j) {
            sum += A[i][j] * v[j];
        }
        result[i] = sum;
    }
}

// Matrix inversion using LUP Decomposition (double)
int invertMatrix_d(const Matrix5x5 A, Matrix5x5 A_inv, int n) {
    if (n != 5) return 1;
    Matrix5x5 L, U;
    int P[5];
    memcpy(U, A, sizeof(Matrix5x5));
    memset(L, 0, sizeof(Matrix5x5));
    for (int i = 0; i < n; ++i) P[i] = i;

    for (int k = 0; k < n; ++k) {
        int pivot_row = k;
        double max_val = 0.0;
        for (int i = k; i < n; ++i) {
            if (fabs(U[i][k]) > max_val) {
                max_val = fabs(U[i][k]);
                pivot_row = i;
            }
        }
        if (max_val < 1e-15) { 
            return 1; // Singular matrix
        }

        if (pivot_row != k) {
            int temp_p = P[k]; P[k] = P[pivot_row]; P[pivot_row] = temp_p;
            double temp_row[5];
            memcpy(temp_row, U[k], n * sizeof(double));
            memcpy(U[k], U[pivot_row], n * sizeof(double));
            memcpy(U[pivot_row], temp_row, n * sizeof(double));
            for (int j = 0; j < k; ++j) {
                double temp_l = L[k][j]; L[k][j] = L[pivot_row][j]; L[pivot_row][j] = temp_l;
            }
        }

        L[k][k] = 1.0;
        if (fabs(U[k][k]) < 1e-15) return 1;
        double pivot_inv = 1.0 / U[k][k];
        if (!is_finite_double(pivot_inv)) return 1;

        for (int i = k + 1; i < n; ++i) {
            L[i][k] = U[i][k] * pivot_inv;
            U[i][k] = 0.0;
            for (int j = k + 1; j < n; ++j) {
                U[i][j] -= L[i][k] * U[k][j];
            }
        }
    }

    for (int j_col = 0; j_col < n; ++j_col) {
        Vector5 b, y, x;
        for (int i = 0; i < n; ++i) b[i] = (P[i] == j_col) ? 1.0 : 0.0;

        for (int i = 0; i < n; ++i) {
            y[i] = b[i];
            for (int k = 0; k < i; ++k) {
                y[i] -= L[i][k] * y[k];
            }
            if (!is_finite_double(y[i])) return 1;
        }

        for (int i = n - 1; i >= 0; --i) {
            x[i] = y[i];
            for (int k = i + 1; k < n; ++k) {
                x[i] -= U[i][k] * x[k];
            }
            if (fabs(U[i][i]) < 1e-15) return 1;
            x[i] /= U[i][i];
            if (!is_finite_double(x[i])) return 1;
        }
        for (int i = 0; i < n; ++i) A_inv[i][j_col] = x[i];
    }
    return 0; // Success
}


// Un-scales the 'x' vector from the solver into physical parameters
// x_scaled is [Iph_ref, log10(Is0_ref), A, log10(Rs), Rp]
void unscale_params_d(const Vector5 x_scaled, double *Iph_param, double *Is0_param, double *A_param, double *Rs_param, double *Rp_param) {
    *Iph_param = x_scaled[0];
    *A_param   = x_scaled[2];
    *Rp_param  = x_scaled[4];

    // Convert from log10 scale back to linear scale
    *Rs_param = pow(10.0, x_scaled[3]);
    
    // Clamp log10(Is0) to a reasonable range and convert back
    double clamped_log10_Is0 = fmax(-20.0, fmin(x_scaled[1], -6.0)); 
    *Is0_param = pow(10.0, clamped_log10_Is0);
    if (*Is0_param < 1e-25) *Is0_param = 1e-25; // Prevent underflow
}


/**
 * @brief Calculates the residual vector F(x) for the 5-parameter model.
 * The goal of the solver is to find an 'x' that makes F(x) = [0,0,0,0,0].
 * * @param x_scaled The scaled parameter vector: [Iph, log10(Is0), A, log10(Rs), Rp]
 * @param F        The output residual vector (5x1).
 */
void residuals_2_20_double(const Vector5 x_scaled, Vector5 F) {
    double Iph_param, Is0_param, A_param, Rs_param, Rp_param;
    unscale_params_d(x_scaled, &Iph_param, &Is0_param, &A_param, &Rs_param, &Rp_param);

    double A_safe = fmax(0.1, A_param);
    double Rp_safe = fmax(0.01, Rp_param);

    // --- Intermediate Terms ---
    double C_denominator = A_safe * SOLVER_k_d * SOLVER_TREF_K_d;
    if (fabs(C_denominator) < 1e-35) C_denominator = 1e-35;
    double C = SOLVER_q_d / C_denominator; 
    double cap = 700.0; 

    // Arguments for the exponentials
    double arg_sc = cap_exp_arg_double(C * Rs_param * SOLVER_ISC_MOD_REF_d); 
    double arg_oc = cap_exp_arg_double(C * SOLVER_VOC_MOD_REF_d); 
    double arg_mp = cap_exp_arg_double(C * (Rs_param * SOLVER_IMP_MOD_REF_d + SOLVER_VMP_MOD_REF_d));
    double arg_eq5 = cap_exp_arg_double(C * SOLVER_ISC_MOD_REF_d); 

    double exp_sc = exp(arg_sc);
    double exp_oc = exp(arg_oc);
    double exp_mp = exp(arg_mp);
    double exp_eq5 = exp(arg_eq5); // Exponential for F[3]

    if (!is_finite_double(exp_sc) || !is_finite_double(exp_oc) || !is_finite_double(exp_mp) || !is_finite_double(exp_eq5)) {
        for (int i = 0; i < 5; ++i) F[i] = DBL_MAX / 10.0; 
        return;
    }

    // --- Residual Equations ---

    double term0_Rp = (Rs_param * SOLVER_ISC_MOD_REF_d) / Rp_safe;
    F[0] = SOLVER_ISC_MOD_REF_d - (Iph_param - Is0_param * (exp_sc - 1.0) - term0_Rp);

    double term1_Rp = (SOLVER_VOC_MOD_REF_d) / Rp_safe;
    F[1] = 0 - (Iph_param - Is0_param * (exp_oc - 1.0) - term1_Rp);

    double term2_Rp = (Rs_param * SOLVER_IMP_MOD_REF_d + SOLVER_VMP_MOD_REF_d) / Rp_safe;
    F[2] = SOLVER_IMP_MOD_REF_d - (Iph_param - Is0_param * (exp_mp - 1.0) - term2_Rp);

    double term4_coeff = (SOLVER_q_d * Is0_param * Rp_safe * (Rs_param - Rp_safe)) / C_denominator; 
    F[3] = Rs_param + term4_coeff * exp_eq5; 

    double term5_bracket_coeff = 1.0 + C * (SOLVER_VMP_MOD_REF_d - Rs_param * SOLVER_IMP_MOD_REF_d); 
    double term5_inside_bracket = term5_bracket_coeff * exp_mp; 
    double term5_Rp = (2.0 * SOLVER_VMP_MOD_REF_d) / Rp_safe;
    
    F[4] = Iph_param - term5_Rp + Is0_param - Is0_param * term5_inside_bracket;

    for (int i = 0; i < 5; ++i) {
        if (!is_finite_double(F[i])) {
            F[i] = DBL_MAX / 10.0; 
        }
    }
}

/**
 * @brief Computes the Jacobian matrix J (partial derivatives) of F(x).
 * * Uses the finite difference method: J[i][j] = (F_i(x + h) - F_i(x)) / h
 * This is numerically sensitive and benefits greatly from 64-bit 'double'.
 */
void compute_jacobian_double(const Vector5 x_scaled, Matrix5x5 J) {
    static Vector5 f_x, f_x_h; // F(x) and F(x+h)
    Vector5 x_temp; 
    const int n_params = 5;
    const double min_safe_val = 1e-35; 

    residuals_2_20_double(x_scaled, f_x); // Calculate F(x) once
    for (int i = 0; i < n_params; ++i) {
        if (!is_finite_double(f_x[i]) || fabs(f_x[i]) > DBL_MAX / 10.0) {
            // F(x) is invalid, return identity matrix as fallback
            for (int r = 0; r < n_params; ++r) for (int c = 0; c < n_params; ++c) J[r][c] = (r == c) ? 1.0 : 0.0;
            return;
        }
    }

    // Iterate over each parameter j to compute the j-th column of J
    for (int j = 0; j < n_params; ++j) {
        memcpy(x_temp, x_scaled, sizeof(Vector5));

        // Calculate step size 'h'
        double h = JACOBIAN_EPSILON_d * (fabs(x_temp[j]) > 1.0 ? fabs(x_temp[j]) : 1.0);
        if (j == 1 || j == 3) { // Special handling for log-scaled parameters
             h = JACOBIAN_EPSILON_d; 
        }
        
        if (fabs(h) < MIN_ABS_STEP_d) {
            h = (x_temp[j] >= 0.0) ? MIN_ABS_STEP_d : -MIN_ABS_STEP_d;
            if (fabs(x_temp[j]) < MIN_ABS_STEP_d) h = MIN_ABS_STEP_d;
        }
        if (fabs(h) < min_safe_val) {
            for (int i = 0; i < n_params; ++i) J[i][j] = 0.0;
            continue;
        }

        x_temp[j] += h; // Take the step
        residuals_2_20_double(x_temp, f_x_h); // Calculate F(x+h)
        double h_inv = 1.0 / h;

        // Compute the derivative column
        for (int i = 0; i < n_params; ++i) {
            if (is_finite_double(f_x_h[i]) && is_finite_double(f_x[i]) && fabs(f_x_h[i]) < DBL_MAX/10.0 && fabs(f_x[i]) < DBL_MAX/10.0) {
                J[i][j] = (f_x_h[i] - f_x[i]) * h_inv;
                if (!is_finite_double(J[i][j])) {
                    J[i][j] = 0.0;
                }
            } else {
                J[i][j] = 0.0;
            }
        }
    }
}

/**
 * @brief Core Levenberg-Marquardt (LM) solver.
 * * This function iteratively tries to find a solution 'x' that
 * minimizes the norm of the residual vector F(x).
 */
double solve_parameters_levenberg_marquardt(const Vector5 initial_x_scaled, Vector5 final_x_scaled) {
    Vector5 x_scaled;
    memcpy(x_scaled, initial_x_scaled, sizeof(Vector5));

    static Vector5 F, g, delta_x_scaled; 
    static Matrix5x5 J, J_T, H_approx, H_inv;
    const int max_iterations = 250; 
    double lambda = 0.01; // Initial damping factor
    double nu = 2.0;
    double current_norm_F = DBL_MAX;
    const int n_params = 5;
    const double max_lambda = 1e10; 

    // Bounds for scaled parameters: [Iph_ref, log10(Is0_ref), A, log10(Rs), Rp]
    const double lower_bounds_scaled[5] = {0.8 * SOLVER_ISC_MOD_REF_d, -25.0, 0.8*ns, log10(1e-4), 10.0};
    const double upper_bounds_scaled[5] = {1.2 * SOLVER_ISC_MOD_REF_d, -5.0, 2.5*ns, log10(10), 100000.0};

    for (int iter = 0; iter < max_iterations; ++iter) {
        // 1. Calculate Residuals F(x)
        residuals_2_20_double(x_scaled, F);
        double current_norm_F_sq = 0.0;
        bool residuals_ok = true;
        for (int i = 0; i < n_params; ++i) {
            if (!is_finite_double(F[i]) || fabs(F[i]) > DBL_MAX / 10.0) {
                residuals_ok = false;
                break;
            }
            current_norm_F_sq += F[i] * F[i];
        }

        if (residuals_ok && current_norm_F_sq < DBL_MAX / 10.0) {
            current_norm_F = sqrt(current_norm_F_sq);
        } else {
            residuals_ok = false;
            current_norm_F = DBL_MAX;
        }

        // 2. Check for convergence
        if (current_norm_F < CONVERGENCE_TOL_d) break;

        // 3. Calculate Jacobian J(x) and Gradient g = J^T * F
        compute_jacobian_double(x_scaled, J);
        transposeMatrix_d(J, J_T, n_params, n_params);
        matrixVectorMultiply_d(J_T, F, g, n_params, n_params);

        bool hessian_solved = false;
        int attempts = 0;
        const int max_attempts = 5; 

        // 4. Try to solve the LM step: (J^T*J + lambda*I) * delta_x = -g
        while (!hessian_solved && attempts < max_attempts && lambda < max_lambda) {
            attempts++;
            // 5. Approximate Hessian H = J^T*J and apply damping H_LM = H + lambda*I
            matrixMultiply_d(J_T, J, H_approx, n_params, n_params, n_params);
            for (int i = 0; i < n_params; ++i) {
                H_approx[i][i] += lambda;
                for (int j = 0; j < n_params; ++j) if (!is_finite_double(H_approx[i][j])) goto increase_lambda_inner;
            }

            // 6. Invert H_LM
            if (invertMatrix_d(H_approx, H_inv, n_params) == 0) {
                bool inv_ok = true;
                for (int i = 0; i < n_params; ++i) for (int j = 0; j < n_params; ++j) if (!is_finite_double(H_inv[i][j])) inv_ok = false;

                if (inv_ok) {
                    Vector5 neg_gradient;
                    bool grad_ok = true;
                    for (int i = 0; i < n_params; ++i) {
                        if (!is_finite_double(g[i])) grad_ok = false;
                        neg_gradient[i] = -g[i];
                    }

                    // 7. Calculate step delta_x = H_LM_inv * (-g)
                    if (grad_ok) {
                        matrixVectorMultiply_d(H_inv, neg_gradient, delta_x_scaled, n_params, n_params);
                        hessian_solved = true;
                    } else {
                        goto increase_lambda_inner;
                    }
                } else {
                    goto increase_lambda_inner;
                }
            } else {
                goto increase_lambda_inner;
            }
            continue; 

        // If inversion fails, increase damping and retry
        increase_lambda_inner:
            lambda *= nu;
            nu = fmin(nu * 2.0, 10.0);
        }

        if (!hessian_solved) break; // Failed to solve step

        double step_norm_sq = 0.0; double x_norm_sq = 0.0; bool step_ok = true;
        for (int i = 0; i < n_params; ++i) {
            if (!is_finite_double(delta_x_scaled[i])) { step_ok = false; break; }
            step_norm_sq += delta_x_scaled[i] * delta_x_scaled[i];
            x_norm_sq += x_scaled[i] * x_scaled[i];
        }
        double step_norm = sqrt(step_norm_sq); double x_norm = sqrt(x_norm_sq);

        // Check for step size convergence
        if (iter > 0 && step_norm < STEP_TOL_d * (x_norm + STEP_TOL_d)) break;

        // 8. Calculate candidate x_new = x + delta_x (with bounds)
        Vector5 x_candidate_scaled;
        for (int i = 0; i < n_params; ++i) {
            x_candidate_scaled[i] = x_scaled[i] + delta_x_scaled[i];
            x_candidate_scaled[i] = fmax(lower_bounds_scaled[i], fmin(x_candidate_scaled[i], upper_bounds_scaled[i]));
        }

        // 9. Evaluate gain ratio (rho) to see if the step was good
        Vector5 F_new;
        residuals_2_20_double(x_candidate_scaled, F_new);
        double new_norm_F_sq = 0.0; bool new_residuals_ok = true;
        for (int i = 0; i < n_params; ++i) {
            if (!is_finite_double(F_new[i]) || fabs(F_new[i]) > DBL_MAX / 10.0) { new_residuals_ok = false; break; }
            new_norm_F_sq += F_new[i] * F_new[i];
        }

        double rho = -1.0;
        if (new_residuals_ok) {
            double actual_reduction = current_norm_F_sq - new_norm_F_sq;
            double predicted_denominator = 0.0;
            for (int i = 0; i < n_params; ++i) {
                predicted_denominator += delta_x_scaled[i] * (lambda * delta_x_scaled[i] - g[i]);
            }

            if (predicted_denominator > 1e-10) {
                rho = actual_reduction / predicted_denominator;
            } else if (actual_reduction > 1e-10) {
                rho = 10.0;
            }
        }

        // 10. Adjust damping factor lambda based on rho
        if (rho > 0.0) {
            // Good step: accept x_new, decrease damping
            memcpy(x_scaled, x_candidate_scaled, sizeof(Vector5));
            lambda = lambda * fmax(1.0 / 3.0, 1.0 - pow(2.0 * rho - 1.0, 3));
            nu = 2.0;
        } else {
            // Bad step: reject x_new, increase damping
            lambda *= nu;
            nu = fmin(nu * 2.0, 10.0);
        }
        lambda = fmax(1e-10, fmin(lambda, max_lambda));
    }

    // Return the final solution and its residual norm
    memcpy(final_x_scaled, x_scaled, sizeof(Vector5));
    residuals_2_20_double(final_x_scaled, F);
    double final_norm_sq = 0.0; bool final_residuals_ok = true;
    for (int i = 0; i < 5; ++i) {
        if (!is_finite_double(F[i])) final_residuals_ok = false;
        final_norm_sq += F[i] * F[i];
    }

    if (!final_residuals_ok || final_norm_sq > DBL_MAX / 10.0) return DBL_MAX;
    return sqrt(final_norm_sq);
}


/**
 * @brief Main function to run the 5-parameter solver.
 *
 * This sets up the initial guess (x0) and calls the
 * Levenberg-Marquardt solver. On success, it populates the
 * global 32-bit float parameters (Iph_ref, Is0_ref, etc.).
 */
void solve_parameters() {
    // Initial guess (x0), in double
    Vector5 initial_x_scaled_d;
    initial_x_scaled_d[0] = (double)Isc_mod_ref; // Iph_ref ~ Isc_cell_ref
    initial_x_scaled_d[1] = log10(1e-9);          // log10(Is0_ref)
    initial_x_scaled_d[2] = ns;                  // A (a common guess)
    initial_x_scaled_d[3] = log10(0.05);         // log10(Rs)
    initial_x_scaled_d[4] = 500.0;                // Rp

    printk("Starting 5-parameter model calculation (double-precision).\n");

    Vector5 best_params_scaled;
    double min_residual_norm = solve_parameters_levenberg_marquardt(initial_x_scaled_d, best_params_scaled);

    // Check if solver failed
    if (!is_finite_double(min_residual_norm) || min_residual_norm > 1e-5) { // 1e-5 is a reasonable success tolerance
        printk("\nERROR: Double-precision solver failed to converge (Norm = %.3e).\n", min_residual_norm);
        return;
    }

    // On success, unscale the parameters
    double final_Iph_ref, final_Is0_ref, final_A, final_Rs, final_Rp;
    unscale_params_d(best_params_scaled, &final_Iph_ref, &final_Is0_ref, &final_A, &final_Rs, &final_Rp);

    printk(">>> Updating global parameters from double-precision solver results.\n");
    // Cast the 64-bit (double) results to 32-bit (float) for runtime use
    Iph_ref = (float32_t)final_Iph_ref;
    Is0_ref = (float32_t)final_Is0_ref;
    A = (float32_t)final_A;
    Rs = (float32_t)final_Rs;
    Rp = (float32_t)final_Rp;
}


//====================================================================
// --- START: (Emulator Logic) Implementations (32-bit Float) ---
// These functions run continuously to operate the emulator.
//====================================================================

// Helper to update history buffer (FIFO)
void updateHistory(float32_t *hist, float32_t val) {
    hist[0] = hist[1];
    hist[1] = val;
}

// Pre-calculates the Voltage points 
// Concentrates points around the knee (Vmp) for high precision.
void computeVoltagePoints(float32_t V[N_POINTS]) {
    
    // Adjust Vmp proportionally to the drop in Voc
    float32_t Vmp_estimation = Vmp_mod_ref * (Voc_estimation / Voc_mod_ref);

    // 1. Define the Factors Vector (25 points)
    // Dense distribution around 1.0 (Vmp) to capture the knee curve
    const float32_t factors[25] = {
        0.00f, 0.20f, 0.40f, 0.60f, 0.70f, 0.80f, 0.85f, 0.88f, 
        0.90f, 0.92f, 0.93f, 0.94f, 0.95f, 0.96f, 0.97f, 0.98f, 
        0.99f, 1.00f, 1.01f, 1.02f, 1.03f, 1.04f, 1.05f, 1.08f, 
        1.12f
    };

    // 2. Generate points and apply filtering
    int index = 0;
    
    // Loop through factors to calculate V_mod candidates
    for (int i = 0; i < 25; i++) {
        float32_t candidate_V = factors[i] * Vmp_estimation;

        // Ensure strictly increasing and below Voc_est
        if (candidate_V >= 0.0f && candidate_V < Voc_estimation) {
            // Safety check to ensure monotonicity (avoid numeric noise)
            if (index == 0 || candidate_V > V[index - 1]) {
                V[index++] = candidate_V;
            }
        }
    }

    // 3. Fill the remaining slots with Voc_estimation
    // Duplicates at the end (Voc, Voc) produce segments with slope 0, which is handled safely by computeSegments.
    while (index < N_POINTS) {
        V[index++] = Voc_estimation;
    }
}

/**
 * @brief Solves the single-diode model equation for I, given V.
 * * This function calculates the I-V curve at the *current*
 * operating conditions (S, T) using the 5 parameters found
 * by the solver. It uses a Newton-Raphson iterative solver.
 * * @param V The cell voltage (V_module).
 * @return float32_t The cell current (I_module).
 */
float32_t solve_I_V_final_model(float32_t V) {
    // 1. Calculate Photons Current (Iph) at current S, T
    float32_t Iph = Iph_ref * (S / Sref) * (1.0f + alpha * (T - Tref));

    // 2. Calculate Saturation Current (Is) at current T
    float32_t Eg_T_eV = E_G0 - (k1 * T * T) / (T + k2);
    float32_t temp_diff = (1.0f / Tref) - (1.0f / T);
    float32_t exponent_term = (ns*q / (A * k)) * Eg_T_eV * temp_diff;

    if (!isfinite(exponent_term)) {
        exponent_term = 0.0f;
    }

    float32_t Is = Is0_ref * powf(T / Tref, 3.0f) * expf(exponent_term);

    // 3. Iterative Newton-Raphson solver for I
    // Equation: f(I) = Iph - Is*(exp(q*(V+I*Rs)/(A*k*T)) - 1) - (V+I*Rs)/Rp - I = 0
    float32_t I = Iph; // Initial guess
    const float32_t cap = 80.0f; 

    for (int it = 0; it < 30; it++) {
        float32_t V_diode = V + I * Rs;
        
        float32_t A_safe = fmaxf(A, 0.1f); // Avoid division by zero
        float32_t V_T = A_safe * k * T / q; // Thermal voltage
        if (V_T < 1e-9f) V_T = 1e-9f; 
        
        float32_t arg_exp = fminf(V_diode / V_T, cap);
        float32_t expo = expf(arg_exp);

        float32_t Rp_safe = fmaxf(Rp, 0.01f); // Avoid division by zero
        
        // Calculate f(I) and df/dI
        float32_t f = Iph - Is * (expo - 1.0f) - V_diode / Rp_safe - I;
        float32_t df = -Is * expo * (Rs / V_T) - Rs / Rp_safe - 1.0f;

        if (fabsf(df) < 1e-9f) df = -1.0f; // Avoid division by zero
        
        // Newton-Raphson step
        float32_t dI = -f / df;
        I = I + dI;

        if (fabsf(dI) < 1e-4f) { // Check for convergence
            break;
        }
    }
    return I;
}

// Calculates the 11 current points corresponding to the 11 voltage points
void computeCurrentPoints(float32_t I[N_POINTS], float32_t V[N_POINTS]) {
    for (int i = 0; i < N_POINTS; i++) {
        float32_t V_cell = V[i]; 
        float32_t I_cell = solve_I_V_final_model(V_cell); 
        I[i] = I_cell * np; 
    }
}

// Calculates the (a, b) coefficients for the 10 linear segments
void computeSegments(struct Segment segments[N_POINTS - 1], float32_t V[N_POINTS], float32_t I[N_POINTS]) {
    for (int i = 0; i < N_POINTS - 1; i++) {
        float32_t deltaV = V[i + 1] - V[i];
        if (fabsf(deltaV) < 1e-9f) deltaV = 1e-9f; // Avoid division by zero
        
        segments[i].a = (I[i + 1] - I[i]) / deltaV; // Slope
        segments[i].b = I[i] - segments[i].a * V[i]; // Intercept
        segments[i].Vmin = V[i]; 
        segments[i].Vmax = V[i + 1]; 
    }
}

/**
 * @brief Finds the intersection of the load line and the I-V curve.
 * @param segments The piecewise linear I-V curve.
 * @param n        Number of segments (N_POINTS - 1).
 * @param Vmes     The measured "Test Point" voltage.
 * @param Imes     The measured "Test Point" current.
 * @param Vint     Output: The intersection voltage (V*).
 * @param Iint     Output: The intersection current (I*).
 */
void findIntersection(struct Segment segments[N_POINTS - 1], int n,
                      float32_t Vmes, float32_t Imes,
                      float32_t &Vint, float32_t &Iint) {
    
    float32_t a_mes = 0.0f; // Load line slope (1/R_load)
    if (Vmes > 1e-3f) { // Avoid division by zero
        a_mes = Imes / Vmes; 
    } else {
        // Handle short-circuit case
        Vint = 0.0f; 
        Iint = segments[0].b; // Intersects at V=0 (Isc)
        return;
    }

    // Find which segment the load line intersects
    for (int i = 0; i < n; i++) {
        float32_t delta_a = segments[i].a - a_mes;
        if (fabsf(delta_a) < 1e-9f) continue; // Parallel lines

        // Vint = (b_mes - b_pv) / (a_pv - a_mes)
        // Since b_mes = 0 (load line passes through origin), Vint = -b_pv / delta_a
        Vint = -segments[i].b / delta_a;
        Iint = a_mes * Vint;

        // Check if the intersection is within this segment's voltage range
        if (Vint >= segments[i].Vmin && Vint <= segments[i].Vmax) {
            return; // Found it
        }
    }
    // If no intersection found (e.g., open circuit), default to Voc
    Vint = segments[n - 1].Vmax; 
    Iint = segments[n - 1].a * Vint + segments[n - 1].b;
}

//====================================================================
// --- START: (Main Program Logic) ---
//====================================================================

/**
 * @brief Runs once at startup.
 */
void setup_routine() {
    // 1. Initialize hardware
    shield.power.initBuck(ALL);

    // 2. Configure sensors
    shield.sensors.enableDefaultTwistSensors();
    shield.sensors.setConversionParametersLinear(I1_LOW, 0.0056f, -11.537f);
    shield.sensors.setConversionParametersLinear(V1_LOW, 0.0456f, -92.69f);
    shield.sensors.setConversionParametersLinear(I2_LOW, 0.0046f, -9.3977f);
    shield.sensors.setConversionParametersLinear(V2_LOW, 0.0453f, -92.061f);
    shield.sensors.setConversionParametersLinear(I_HIGH, 0.0046f, -9.1739f);
    shield.sensors.setConversionParametersLinear(V_HIGH, 0.0299f, 0.2921f);

    // 3. Initialize PID controller
    pid.init(pidParams);

    // 4. *** RUN THE 64-BIT SOLVER ***
    // This is the most critical part of the setup.
    // It populates the global 32-bit parameters (Iph_ref, Is0_ref, etc.)
    printk("Starting 5-parameter model calculation (double-precision)...\n");
    solve_parameters(); 
    printk("Parameter calculation finished.\n");

    // 5. Pre-calculate the I-V curve segments for the emulator
    float32_t Vp[N_POINTS], Ip[N_POINTS];
    computeVoltagePoints(Vp);
    computeCurrentPoints(Ip, Vp); 
    computeSegments(segments, Vp, Ip);

    // 6. Print the final solved parameters
    printk("\n");
    printk("| --- Calculated PV Model Parameters --- |\n");
    printk("| Iph_ref = %.6f A |\n", Iph_ref);
    printk("| Is0_ref = %.2e A |\n", Is0_ref);
    printk("| A       = %.6f |\n", A);
    printk("| Rs      = %.6f Ohms |\n", Rs);
    printk("| Rp      = %.6f Ohms |\n", Rp);
    printk("| Conditions: S = %.1f W/m^2, T = %.1f K |\n", S, T - 273.15f); // Display T in Celsius
    printk("\n");

    // 7. Create and start all tasks
    uint8_t appTaskId = task.createBackground(loop_application_task);
    uint8_t commTaskId = task.createBackground(loop_communication_task);
    task.createCritical(loop_critical_task, controlTaskPeriod);
    task.startBackground(appTaskId);
    task.startBackground(commTaskId);
    task.startCritical();
}

/**
 * @brief Background task for serial communication (Menu Handler).
 */
void loop_communication_task() {
    while (1) {
        receivedChar = console_getchar(); // Wait for user input
        switch (receivedChar) {
            case 'h': // Print help menu
                printk(" ________________________________________\n");
                printk("| ---- Buck Voltage Mode MENU ---- |\n");
                printk("| Press 'i' for idle mode |\n");
                printk("| Press 'p' for power mode |\n");
                printk("| Press 'e' for emulator mode |\n");
                printk("| Press 'u' to increase voltage ref|\n");
                printk("| Press 'd' to decrease voltage ref|\n");
                printk("| Press 's' to show solver solution|\n"); // Added 's'
                printk(" ________________________________________\n");
                break;
            case 'i': // IDLE mode
                printk("Idle mode\n");
                mode = MODE_IDLE;
                break;
            case 'p': // POWER mode
                printk("Power mode\n");
                mode = MODE_POWER;
                dutyCycle = dutyCycleP;
                break;
            case 'e': // EMULATOR mode
                printk("Emulator mode\n");
                dutyCycleP = dutyCycle; // Save current duty
                mode = MODE_EMULATOR;
                // Reset the emulator state machine
                dutyCycle = pid.calculateWithReturn(voltageReferenceE_test_point, lowVoltage1);
                shield.power.setDutyCycle(ALL, dutyCycle);
                currentHistory[0] = currentHistory[1] = 0.0f; 
                emulatorSteadyState = false; 
                waitingSteadyState = false; 
                loadChangeCounter = 0; 
                break;
            case 'u': // Manually increase duty cycle (for testing)
                dutyCycle = fminf(dutyCycle + 0.1f, 1.0f);
                break;
            case 'd': // Manually decrease duty cycle (for testing)
                dutyCycle = fmaxf(dutyCycle - 0.1f, 0.0f);
                break;
            
            // --- NEW CASE ADDED ---
            case 's': // Show solver solution
                printk("\n| --- Current Solver Solution --- |\n");
                printk("| Iph_ref = %.6f A |\n", Iph_ref);
                printk("| Is0_ref = %.2e A |\n", Is0_ref);
                printk("| A       = %.6f |\n", A);
                printk("| Rs      = %.6f Ohms |\n", Rs);
                printk("| Rp      = %.6f Ohms |\n", Rp);
                printk(" --------------------------------- \n");
                break;
        }
    }
}

/**
 * @brief Background task for high-level logic and printing.
 *
 * This task runs at a lower priority. It averages measurements
 * and manages the emulator's state machine.
 */
void loop_application_task() {
    static uint32_t elapsed = 0;
    elapsed += 100; // Timekeeping

    // Turn LED on if not in IDLE mode
    if (mode == MODE_IDLE) spin.led.turnOff();
    else spin.led.turnOn();

    uint32_t period = 500; // Print period (us)
    
    // This block runs every 'period' (500us)
    if (elapsed >= period && measurementCount > 0) {
        
        // Calculate average values from the accumulators
        float32_t avgLowI1 = sumLowCurrent1 / measurementCount; 
        float32_t avgLowV1 = sumLowVoltage1 / measurementCount; 
        float32_t avgLowI2 = sumLowCurrent2 / measurementCount; 
        float32_t avgLowV2 = sumLowVoltage2 / measurementCount; 
        float32_t avgHighV = sumHighVoltage / measurementCount; 
        float32_t avgLowI = avgLowI1 + avgLowI2; // Total low-side current
        float32_t avgLowV = (avgLowV1 + avgLowV2) * 0.5f; // Avg low-side voltage

        // --- Print debug info for POWER mode ---
        if (mode == MODE_POWER) {
            printk("lowCurrent [Ilow1 + Ilow2]: %f A\n", avgLowI);
            printk("lowVoltage [(Vlow1 + Vlow2)/2]: %f V\n", avgLowV);
            printk("highVoltage: %f V\n", avgHighV);
        }

        // --- EMULATOR STATE MACHINE ---

        // STATE 1: Initial calibration phase
        if (mode == MODE_EMULATOR && !emulatorSteadyState) {
            updateHistory(currentHistory, avgLowI); // Update current history
            
            // Wait for current to be stable
            if (fabsf(currentHistory[0]) > 0.01f && (fabsf(currentHistory[1] - currentHistory[0]) / currentHistory[0]) < CURRENT_STABILITY_THRESHOLD) {
                
                // --- STABLE: Find the intersection ---
                float32_t Vint, Iint;
                // Vmes = avgLowV, Imes = currentHistory[1] (stable current)
                findIntersection(segments, N_POINTS - 1, avgLowV, currentHistory[1], Vint, Iint);
                
                // Set the *new* PID target to the intersection voltage
                voltageReferenceE = Vint; 
                
                float32_t loadLineSlope = (avgLowV != 0.0f) ? (avgLowI / avgLowV) : 0.0f; 

                printk("Test Point: V_test = %f V, I_test = %f A\n", avgLowV, avgLowI);
                printk("Load Line Equation: I = %f * V\n", loadLineSlope);
                printk("Intersection Found: Vo* = %f V, Io* = %f A\n", Vint, Iint);
                
                // Store this as the last known good state
                lastSteadyVoltage = Vint; 
                lastSteadyCurrent = Iint; 
                
                // --- Transition to STATE 2 ---
                emulatorSteadyState = true; 
                loadChangeCounter = 0; 
            } else {
                printk("Measurements unstable for initial measurement phase.\n");
            }
        }

        // STATE 2: Monitoring phase (steady state)
        if (mode == MODE_EMULATOR && emulatorSteadyState) {
            float32_t candidateI = avgLowI; 
            float32_t candidateV = avgLowV; 
            float32_t currentDuty = dutyCycle; 

            // Check if user is manually changing the duty cycle
            if (fabsf(currentDuty - lastDutyCycle) > DUTY_CYCLE_CHANGE_THRESHOLD) {
                waitingSteadyState = true;
                loadChangeCounter = 0;
                lastSteadyVoltage = candidateV;
                lastSteadyCurrent = candidateI;
            }
            // If we were waiting for stability after a manual change
            else if (waitingSteadyState) {
                if (fabsf(currentHistory[0]) > 0.01f && (fabsf(currentHistory[1] - currentHistory[0]) / currentHistory[0]) < CURRENT_STABILITY_THRESHOLD) {
                    waitingSteadyState = false; // Now stable
                    lastSteadyVoltage = candidateV;
                    lastSteadyCurrent = candidateI;
                }
            }
            // Normal monitoring
            else {
                // Check if the load resistance has changed
                float32_t Rlast = (lastSteadyCurrent > 0.0f) ? (lastSteadyVoltage / lastSteadyCurrent) : 0.0f; 
                float32_t Rcand = (candidateI > 0.0f) ? (candidateV / candidateI) : 0.0f; 
                float32_t relErr = (Rlast > 0.0f) ? fabsf(Rcand - Rlast) / Rlast : 0.0f; 
                
                if (relErr > LOAD_RESISTANCE_THRESHOLD) {
                    loadChangeCounter++; // Increment persistence counter
                    
                    // If change is persistent, trigger recalibration
                    if (loadChangeCounter >= REQUIRED_PERSISTENCE_COUNT) {
                        loadChangeCounter = 0;
                        pid.init(pidParams); // Reset PID
                        // Set target back to the "Test Point"
                        dutyCycle = pid.calculateWithReturn(voltageReferenceE_test_point, lowVoltage1);
                        shield.power.setDutyCycle(ALL, dutyCycle);
                        
                        // --- Transition back to STATE 1 ---
                        emulatorSteadyState = false; 
                        printk("Load change detected. Recalibrating...\n");
                    }
                } else {
                    // No change, reset counter and update last known state
                    loadChangeCounter = 0; 
                    lastSteadyVoltage = candidateV;
                    lastSteadyCurrent = candidateI;
                }
            }
            lastDutyCycle = currentDuty; 
        }

        // Reset accumulators for the next 500us cycle
        sumLowCurrent1 = sumLowVoltage1 = sumLowCurrent2 = sumLowVoltage2 =
            sumHighCurrent = sumHighVoltage = 0.0f;
        measurementCount = 0;
        elapsed = 0;
    }
    task.suspendBackgroundMs(100); // Suspend task for 100us
}

/**
 * @brief Critical (real-time) task.
 *
 * This task runs every 'controlTaskPeriod' (100us).
 * It reads sensors, runs the PID controller, and sets the PWM.
 * It must be very fast.
 */
void loop_critical_task() {
    // 1. Read all sensors
    lowCurrent1 = shield.sensors.getLatestValue(I1_LOW);
    lowVoltage1 = shield.sensors.getLatestValue(V1_LOW);
    lowVoltage2 = shield.sensors.getLatestValue(V2_LOW);
    lowCurrent2 = shield.sensors.getLatestValue(I2_LOW);
    highCurrent = shield.sensors.getLatestValue(I_HIGH);
    highVoltage = shield.sensors.getLatestValue(V_HIGH);

    float32_t currentVoltage;
    
    // --- Mode-based Control ---
    
    if (mode == MODE_IDLE) {
        if (pwmEnable) {
            shield.power.stop(ALL);
        }
        pwmEnable = false;
    
    } else if (mode == MODE_POWER) {
        // --- Standard Buck Converter Control ---
        currentVoltage = (lowVoltage1 + lowVoltage2) * 0.5f;
        
        // Run PID: Calculate new duty cycle to match voltageReferenceP
        dutyCycle = pid.calculateWithReturn(voltageReferenceP, currentVoltage);
        shield.power.setDutyCycle(ALL, dutyCycle);

        // Accumulate measurements for the background task
        sumLowCurrent1 += lowCurrent1;
        sumLowVoltage1 += lowVoltage1;
        sumLowCurrent2 += lowCurrent2;
        sumLowVoltage2 += lowVoltage2;
        sumHighCurrent += highCurrent;
        sumHighVoltage += highVoltage;
        measurementCount++;

        if (!pwmEnable) {
            pwmEnable = true;
            shield.power.start(ALL);
        }
    
    } else { // MODE_EMULATOR
        // --- Emulator Control ---
        currentVoltage = (lowVoltage1 + lowVoltage2) * 0.5f;
        
        // Select the correct PID target based on the state
        float32_t targetVoltage = emulatorSteadyState ? voltageReferenceE : voltageReferenceE_test_point;

        // Run PID: Calculate new duty cycle to match the target voltage
        dutyCycle = pid.calculateWithReturn(targetVoltage, currentVoltage);
        shield.power.setDutyCycle(ALL, dutyCycle);

        // Accumulate measurements for the background task
        sumLowCurrent1 += lowCurrent1;
        sumLowVoltage1 += lowVoltage1;
        sumLowCurrent2 += lowCurrent2;
        sumLowVoltage2 += lowVoltage2;
        sumHighCurrent += highCurrent;
        sumHighVoltage += highVoltage;
        measurementCount++;

        if (!pwmEnable) {
            pwmEnable = true;
            shield.power.start(ALL);
        }
    }
}
//====================================================================
// Entry Point
//====================================================================

int main(void) {
    setup_routine(); // Run setup
    return 0; // setup_routine() starts tasks, main() exits
}
