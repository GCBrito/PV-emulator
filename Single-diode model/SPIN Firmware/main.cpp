#include "ShieldAPI.h"
#include "TaskAPI.h"
#include "SpinAPI.h"
#include "pid.h"
#include "zephyr/console/console.h"
#include <math.h>

//====================================================================
// Global Declarations
//====================================================================

// Loop task declarations
void loop_application_task(); // Main task: averages calculation and application logic
void loop_communication_task(); // Serial communication task
void loop_critical_task(); // Sensor reading and PWM control task

// Parameters for I-V interpolation
#define HISTORY_SIZE 2 // History size to ensure steady-state current
#define N_POINTS 11 // Number of points for I-V curve approximation

// PV module parameters (datasheet)
const float32_t ns = 36;     // Number of series cells
const float32_t np = 1;      // Number of parallel cells

const float32_t Vmp_mod_ref = 18.2f; // Module voltage at maximum power point [V]
const float32_t Imp_mod_ref =  2.2f; // Module current at maximum power point [A]
const float32_t Voc_mod_ref = 22.0f; // Module open-circuit voltage [V]
const float32_t Isc_mod_ref =  2.37f; // Module short-circuit current [A]

const float32_t Tref = 25.0f + 273.15f; // Reference temperature [K]
const float32_t Sref = 1000.0f; // Reference irradiance [W/m²]

const float32_t alpha = 0.0024885f; // Temperature coefficient of Isc (%/K)

// Operating conditions
static float32_t T = 25.0f + 273.15f; // Current temperature [K]
static float32_t S = 1000.0f;       // Current irradiance [W/m²]

// PV cell parameters 
double Vmp_cell_ref = Vmp_mod_ref / ns; // Cell voltage at maximum power point [V]
double Imp_cell_ref = Imp_mod_ref / np; // Cell currrent at maximum power point [A]
double Voc_cell_ref = Voc_mod_ref / ns; // Cell open-circuit voltage [V]
double Isc_cell_ref = Isc_mod_ref / np; // Cell short-circuit current [A]

// Complete model parameters (determined by the MATLAB code)
static float32_t Iph_ref = 2.370265f; // Reference photogenerated current [A]
static float32_t Is0_ref = 7.99e-08f;  // Reference saturation current [A]
static float32_t A = 1.383022f;      // Ideality factor [-]
static float32_t Rs = 0.004203f;     // Series resistance [Ohms]
static float32_t Rp = 37.572742f;     // Parallel resistance [Ohms]

// Physical constants
const float32_t q = 1.60217662e-19f; // Elementary charge [C]
const float32_t k = 1.38064852e-23f; // Boltzmann constant [J/K]

// --- MATERIAL CONSTANTS ---
const float32_t E_G0 = 1.166f; // Band gap energy at 0K (eV)
const float32_t k1 = 4.73e-4f; // Coefficient k1 (eV/K)
const float32_t k2 = 636.0f;   // Coefficient k2 (K)

// Operating modes
enum modeMenu { MODE_IDLE = 0, MODE_POWER, MODE_EMULATOR };
uint8_t mode = MODE_IDLE;

// Duty cycles and Voltage refereces
float32_t dutyCycle = 0.0; // Control variable
static float32_t dutyCycleP = 0.8f; // Duty cycle in POWER mode
static float32_t voltageReferenceP = 0.8f*50; // Voltage reference in POWER mode [V]
static float32_t voltageReferenceE_test_point = 1.1f*Voc_mod_ref; // Test point voltage reference in EMULATOR mode [V]
static float32_t voltageReferenceE = 0.0; // Voltage reference in EMULATOR mode [V]

// Thresholds and persistence
#define LOAD_RESISTANCE_THRESHOLD 0.02f // Relative variation of R > 2%
#define DUTY_CYCLE_CHANGE_THRESHOLD 0.05f // Duty cycle variation > 0.05
#define CURRENT_STABILITY_THRESHOLD 0.1f // Relative variation of Is of 10% for steady-state
#define REQUIRED_PERSISTENCE_COUNT 2 // Number of cycles to validate an event

// Measured values
static float32_t lowVoltage1, lowVoltage2; // Measured low-side voltages [V]
static float32_t lowCurrent1, lowCurrent2; // Measured low-side currents [A]
static float32_t highCurrent, highVoltage; // Measured high-side voltage [V] and current [A]
static float measData; // Raw read value

static uint32_t controlTaskPeriod = 100; // Critical task period [µs]
static bool pwmEnable = false;
uint8_t receivedChar; // Character received via serial

// PID parameters
static float32_t kp = 0.000215f; // Proportional gain
static float32_t Ti = 7.5175e-5f; // Integral time
static float32_t Td = 0.0f; // Derivative time
static float32_t Np = 0.0f; // Derivative filter coefficient
static float32_t upperBound = 1.0f; // Upper bound for PID output
static float32_t lowerBound = 0.0f; // Lower bound for PID output
static float32_t Ts = controlTaskPeriod * 1e-6f; // PID sampling period [s]
static PidParams pidParams(Ts, kp, Ti, Td, Np, lowerBound, upperBound);
static Pid pid;

// Accumulators for average calculation
static volatile float32_t sumLowCurrent1 = 0.0f, sumLowVoltage1 = 0.0f;
static volatile float32_t sumLowCurrent2 = 0.0f, sumLowVoltage2 = 0.0f;
static volatile float32_t sumHighCurrent = 0.0f, sumHighVoltage = 0.0f;
static volatile uint32_t measurementCount = 0;

// History and flags for EMULATOR MODE
static float32_t currentHistory[HISTORY_SIZE] = {0.0f, 0.0f};
static bool emulatorSteadyState = false;
static bool waitingSteadyState = false;

// Structure representing a linear segment of the I-V curve
struct Segment {
    float32_t a; // Slope of the segment (dI/dV)
    float32_t b; // Y-intercept
    float32_t Vmin; // Minimum valid voltage for the segment [V]
    float32_t Vmax; // Maximum valid voltage for the segment [V]
};
Segment segments[N_POINTS - 1]; // Array of segments

// Last stable point (intersection)
static float32_t lastSteadyVoltage = 0.0f;
static float32_t lastSteadyCurrent = 0.0f;
static float32_t lastDutyCycle = 0.0f;

// Load change detection counter
static uint32_t loadChangeCounter = 0;

//====================================================================
// Utility Functions
//====================================================================

// Updates a circular history of the last two measurements.
//
// Inputs:
// - hist: Pointer to the history array.
// - val: The new value to add to the history.
//
// Returns:
// - void
void updateHistory(float32_t *hist, float32_t val) {
    hist[0] = hist[1];
    hist[1] = val;
}

// Generates an array of voltage points for I-V curve interpolation.
// The points are distributed from 0 to Voc_mod_ref based on fractions of Vmp_mod_ref.
//
// Inputs:
// - V: Array to store the calculated voltage points [V].
//
// Returns:
// - void
void computeVoltagePoints(float32_t V[N_POINTS]) {
    V[0] = 0.0f;
    V[1] = 0.2f * Vmp_cell_ref * ns;
    V[2] = 0.4f * Vmp_cell_ref * ns;
    V[3] = 0.6f * Vmp_cell_ref * ns;
    V[4] = 0.8f * Vmp_cell_ref * ns;
    V[5] = 0.9f * Vmp_cell_ref * ns;
    V[6] = Vmp_cell_ref * ns;
    V[7] = (Vmp_cell_ref * ns + Voc_cell_ref*ns) / 2.0f;
    V[8] = 0.9f * Voc_cell_ref*ns;
    V[9] = 0.95f * Voc_cell_ref*ns;
    V[10] = Voc_cell_ref*ns;
}

// Solves the I-V equation using the final model with temperature-dependent Eg and alpha.
//
// Inputs:
// - V: The voltage for which to solve the current [V].
//
// Returns:
// - The calculated current [A].
float32_t solve_I_V_final_model(float32_t V) {
    // 1. Photocurrent with T and S dependency, using alpha
    float32_t Iph = Iph_ref * (S / Sref) * (1.0f + alpha * (T - Tref));

    // 2. Temperature-dependent Band Gap Energy (Eg) (Varshni's equation)
    // Eg is calculated in eV and then converted to Joules
    float32_t Eg_T_eV = E_G0 - (k1 * T * T) / (T + k2);
    float32_t temp_diff = (1.0f / Tref) - (1.0f / T);
    float32_t exponent_term = (q / (A * k)) * Eg_T_eV * temp_diff;

    // 3. Saturation Current (Is) with T and Eg(T) dependency
    // Using a robust form of the equation
    float32_t Is = Is0_ref * powf(T / Tref, 3.0f) * expf(exponent_term);

    // 4. Iterative solver (Newton-Raphson) for the diode equation
    float32_t I = Iph; // Initial guess

    const float32_t cap = 700.0f;

    for (int it = 0; it < 30; it++) {
        float32_t V_diode = V + I * Rs;
        float32_t arg_exp = fminf(q * V_diode / (A * k * T), cap);
        float32_t expo = expf(arg_exp);

        // Diode equation function f(I) = 0
        float32_t f = Iph - Is * (expo - 1.0f) - V_diode / Rp - I;

        // Derivative df/dI
        float32_t df = -Is * expo * (q * Rs / (A * k * T)) - Rs / Rp - 1.0f;

        // Newton-Raphson step
        float32_t dI = -f / df;
        I = I + dI;

        // Convergence criterion
        if (fabsf(dI) < 1e-6f) {
            break;
        }
    }
    return I;
}

// Calculates the currents corresponding to the voltage points for the I-V curve
// using the final 5-parameter model.
//
// Inputs:
// - I: Array to store the calculated current points [A].
// - V: Array of voltage points [V].
//
// Returns:
// - void
void computeCurrentPoints(float32_t I[N_POINTS], float32_t V[N_POINTS]) {
    for (int i = 0; i < N_POINTS; i++) {
        float32_t V_cell = V[i] / ns; // Voltage per cell [V]
        float32_t I_cell = solve_I_V_final_model(V_cell); // Current per cell [A]
        I[i] = I_cell * np; // Total module current [A]
    }
}

// Generates linear segments (lines) to approximate the I-V curve.
// Each segment covers a range [Vmin, Vmax] with a slope 'a' and intercept 'b'.
//
// Inputs:
// - segments: Array of Segment structures to store the calculated segments.
// - V: Array of voltage points [V].
// - I: Array of current points [A].
//
// Returns:
// - void
void computeSegments(Segment segments[N_POINTS - 1], float32_t V[N_POINTS], float32_t I[N_POINTS]) {
    for (int i = 0; i < N_POINTS - 1; i++) {
        segments[i].a = (I[i + 1] - I[i]) / (V[i + 1] - V[i]); // Slope (dI/dV)
        segments[i].b = I[i] - segments[i].a * V[i]; // Y-intercept
        segments[i].Vmin = V[i]; // Minimum voltage for this segment [V]
        segments[i].Vmax = V[i + 1]; // Maximum voltage for this segment [V]
    }
}

// Finds the intersection between the load line and the I-V segments.
// This function iterates through the pre-calculated I-V curve segments to find
// where the load line (represented by Vmes and Imes) intersects one of the segments.
//
// Inputs:
// - segments: Array of Segment structures representing the I-V curve.
// - n: The number of segments.
// - Vmes: Measured voltage for the load line [V].
// - Imes: Measured current for the load line [A].
// - Vint: Reference to store the intersection voltage [V].
// - Iint: Reference to store the intersection current [A].
//
// Returns:
// - void (updates Vint and Iint by reference)
void findIntersection(Segment segments[N_POINTS - 1], int n,
                      float32_t Vmes, float32_t Imes,
                      float32_t &Vint, float32_t &Iint) {
    float32_t a_mes = (Vmes != 0.0f) ? (Imes / Vmes) : 0.0f; // Slope of the load line (I/V)
    for (int i = 0; i < n; i++) {
        Vint = (0.0f - segments[i].b) / (segments[i].a - a_mes); // Calculate intersection voltage [V]
        Iint = a_mes * Vint; // Calculate intersection current [A]
        if (Vint >= segments[i].Vmin && Vint <= segments[i].Vmax) {
            return; // Valid intersection found
        }
    }
    Vint = Iint = -1.0f; // No intersection found
}

//====================================================================
// Configuration and Initialization
//====================================================================

// Initializes hardware, sensors, PID, and creates tasks.
//
// Inputs:
// - void
//
// Returns:
// - void
void setup_routine() {
    // Hardware initialization
    shield.power.initBuck(ALL);

    // Sensor configuration
    shield.sensors.enableDefaultTwistSensors();
    shield.sensors.setConversionParametersLinear(I1_LOW, 0.0056f, -11.537f);
    shield.sensors.setConversionParametersLinear(V1_LOW, 0.0456f, -92.69f);
    shield.sensors.setConversionParametersLinear(I2_LOW, 0.0046f, -9.3977f);
    shield.sensors.setConversionParametersLinear(V2_LOW, 0.0453f, -92.061f);
    shield.sensors.setConversionParametersLinear(I_HIGH, 0.0046f, -9.1739f);
    shield.sensors.setConversionParametersLinear(V_HIGH, 0.0299f, 0.2921f);

    // PID initialization
    pid.init(pidParams);

    // Calculate I-V segments using the complete model
    float32_t Vp[N_POINTS], Ip[N_POINTS];
    computeVoltagePoints(Vp);
    computeCurrentPoints(Ip, Vp); // Uses the final model
    computeSegments(segments, Vp, Ip);

    // Display model parameters
    printk("\n");
    printk("|  --- Final PV Model Parameters ---  |\n");
    printk("| Iph_ref = %.6f A                    |\n", Iph_ref);
    printk("| Is0_ref = %.2e A                    |\n", Is0_ref);
    printk("| A       = %.6f                      |\n", A);
    printk("| Rs      = %.6f Ohms                 |\n", Rs);
    printk("| Rp      = %.6f Ohms                 |\n", Rp);
    printk("| Conditions: S = %.1f W/m², T = %.1f K |\n", S, T);
    printk("\n");

    // Task creation and startup
    uint8_t appTaskId = task.createBackground(loop_application_task);
    uint8_t commTaskId = task.createBackground(loop_communication_task);
    task.createCritical(loop_critical_task, controlTaskPeriod);
    task.startBackground(appTaskId);
    task.startBackground(commTaskId);
    task.startCritical();
}

//====================================================================
// Serial Communication Loop
//====================================================================

// Background task managing the serial interface.
// - Waits for a character input.
// - Updates the operating mode.
// - Allows modification of PV model parameters.
//
// Inputs:
// - void
//
// Returns:
// - void
void loop_communication_task() {
    while (1) {
        receivedChar = console_getchar(); // Read a character from the console
        switch (receivedChar) {
            case 'h':
                // Display the menu
                printk(" ________________________________________\n");
                printk("| ---- Buck Voltage Mode MENU ---- |\n");
                printk("| Press 'i' for idle mode          |\n");
                printk("| Press 'p' for power mode         |\n");
                printk("| Press 'e' for emulator mode      |\n");
                printk("| Press 'u' to increase voltage ref|\n");
                printk("| Press 'd' to decrease voltage ref|\n");
                printk(" ________________________________________\n");
                break;
            case 'i':
                // Idle mode
                printk("Idle mode\n");
                mode = MODE_IDLE;
                break;
            case 'p':
                // Buck (POWER) mode
                printk("Power mode\n");
                mode = MODE_POWER;
                dutyCycle = dutyCycleP;
                break;
            case 'e':
                // Emulator mode
                printk("Emulator mode\n");
                dutyCycleP = dutyCycle;
                mode = MODE_EMULATOR;
                // PID calculates initial duty cycle based on test point voltage
                dutyCycle = pid.calculateWithReturn(voltageReferenceE_test_point, lowVoltage1);
                shield.power.setDutyCycle(ALL, dutyCycle);
                currentHistory[0] = currentHistory[1] = 0.0f; // Reset current history
                emulatorSteadyState = false; // Reset steady-state flag
                waitingSteadyState = false; // Reset waiting for steady-state flag
                loadChangeCounter = 0; // Reset load change counter
                break;
            case 'u':
                // Increase duty cycle
                dutyCycle += 0.1f;
                break;
            case 'd':
                // Decrease duty cycle
                dutyCycle -= 0.1f;
                break;
        }
    }
}

//====================================================================
// Main Application Task
//====================================================================

// In emulator mode:
// It checks current stability and calculates the intersection between the load line
// and the I-V characteristic to determine the voltageReferenceE value. It displays
// the test point, the load line equation, and the intersection point. It also
// handles dynamic load changes by recalibrating the system and recalculating the
// intersection point.
//
// Inputs:
// - void
//
// Returns:
// - void
void loop_application_task() {
    static uint32_t elapsed = 0;
    elapsed += 100; // Increment elapsed time by 100 µs (critical task period)

    // LED indicator based on mode
    if (mode == MODE_IDLE) spin.led.turnOff();
    else spin.led.turnOn();

    // Averaging calculation period (in us)
    uint32_t period = 500; // Average over 500 us
    if (elapsed >= period && measurementCount > 0) {
        // Calculate averages
        float32_t avgLowI1 = sumLowCurrent1 / measurementCount; // Average of low-side current 1 [A]
        float32_t avgLowV1 = sumLowVoltage1 / measurementCount; // Average of low-side voltage 1 [V]
        float32_t avgLowI2 = sumLowCurrent2 / measurementCount; // Average of low-side current 2 [A]
        float32_t avgLowV2 = sumLowVoltage2 / measurementCount; // Average of low-side voltage 2 [V]
        float32_t avgHighI = sumHighCurrent / measurementCount; // Average of high-side current [A]
        float32_t avgHighV = sumHighVoltage / measurementCount; // Average of high-side voltage [V]
        float32_t avgLowI = avgLowI1 + avgLowI2; // Total average low-side current [A]
        float32_t avgLowV = (avgLowV1 + avgLowV2) * 0.5f; // Average low-side voltage [V]

        if (mode == MODE_POWER) {
            printk("Averages (last %uus):\n", period);
            printk("lowCurrent1: %f A\n", avgLowI1);
            printk("lowVoltage1: %f V\n", avgLowV1);
            printk("lowCurrent2: %f A\n", avgLowI2);
            printk("lowVoltage2: %f V\n", avgLowV2);
            printk("lowCurrent [Ilow1 + Ilow2]: %f A\n", avgLowI);
            printk("lowVoltage [(Vlow1 + Vlow2)/2]: %f V\n", avgLowV);
            printk("highCurrent: %f A\n", avgHighI);
            printk("highVoltage: %f V\n", avgHighV);
        }

        // Emulator: static load
        if (mode == MODE_EMULATOR && !emulatorSteadyState) {
            updateHistory(currentHistory, avgLowI); // Update current history
            // Check for current stability (relative change less than threshold)
            if (fabsf(currentHistory[0]) > 0.01f && (fabsf(currentHistory[1] - currentHistory[0]) / currentHistory[0]) < CURRENT_STABILITY_THRESHOLD) {
                float32_t Vint, Iint;
                // Find the intersection between the load line and the I-V curve
                findIntersection(segments, N_POINTS - 1, avgLowV, currentHistory[1], Vint, Iint);
                voltageReferenceE = Vint; // Set the emulator voltage reference to the intersection voltage
                // Display: test point, load line equation, intersection
                printk("Test Point: V_test = %f V, I_test = %f A\n", avgLowV, avgLowI);
                float32_t loadLineSlope = (avgLowV != 0.0F) ? (avgLowI / avgLowV) : 0.0F; // Slope of the load line (I/V)
                printk("Load Line Equation: I = %f * V\n", loadLineSlope);
                printk("Intersection Found: Vo* = %f V, Io* = %f A\n", Vint, Iint);
                lastSteadyVoltage = Vint; // Store the last steady voltage [V]
                lastSteadyCurrent = Iint; // Store the last steady current [A]
                emulatorSteadyState = true; // Set steady-state flag
                loadChangeCounter = 0; // Reset load change counter
            } else {
                printk("Measurements unstable for initial measurement phase.\n");
            }
        }

        // Emulator: dynamic load adaptation
        if (mode == MODE_EMULATOR && emulatorSteadyState) {
            float32_t candidateI = avgLowI; // Current candidate [A]
            float32_t candidateV = avgLowV; // Voltage candidate [V]
            float32_t currentDuty = dutyCycle; // Current duty cycle

            // If duty cycle changes significantly, re-enter waiting state
            if (fabsf(currentDuty - lastDutyCycle) > DUTY_CYCLE_CHANGE_THRESHOLD) {
                waitingSteadyState = true;
                loadChangeCounter = 0;
                lastSteadyVoltage = candidateV;
                lastSteadyCurrent = candidateI;
            }
            // If waiting for steady state, check current stability
            else if (waitingSteadyState) {
                if (fabsf(currentHistory[0]) > 0.01f && (fabsf(currentHistory[1] - currentHistory[0]) / currentHistory[0]) < CURRENT_STABILITY_THRESHOLD) {
                    waitingSteadyState = false; // Steady state reached
                    lastSteadyVoltage = candidateV;
                    lastSteadyCurrent = candidateI;
                }
            }
            // If already in steady state, check for load resistance change
            else {
                float32_t Rlast = (lastSteadyCurrent > 0.0f) ? (lastSteadyVoltage / lastSteadyCurrent) : 0.0f; // Last stable resistance [Ohms]
                float32_t Rcand = (candidateI > 0.0f) ? (candidateV / candidateI) : 0.0f; // Candidate resistance [Ohms]
                float32_t relErr = (Rlast > 0.0f) ? fabsf(Rcand - Rlast) / Rlast : 0.0f; // Relative error in resistance
                if (relErr > LOAD_RESISTANCE_THRESHOLD) {
                    loadChangeCounter++; // Increment load change counter
                    if (loadChangeCounter >= REQUIRED_PERSISTENCE_COUNT) {
                        loadChangeCounter = 0;
                        dutyCycle = pid.calculateWithReturn(voltageReferenceE_test_point, lowVoltage1);
                        shield.power.setDutyCycle(ALL, dutyCycle);
                        emulatorSteadyState = false; // Exit steady state
                        printk("Recalibration\n");
                    }
                } else {
                    loadChangeCounter = 0; // Reset load change counter
                    lastSteadyVoltage = candidateV;
                    lastSteadyCurrent = candidateI;
                }
            }
            lastDutyCycle = currentDuty; // Store the current duty cycle
        }

        // Reset accumulators
        sumLowCurrent1 = sumLowVoltage1 = sumLowCurrent2 = sumLowVoltage2 =
            sumHighCurrent = sumHighVoltage = 0.0f;
        measurementCount = 0;
        elapsed = 0;
    }
    task.suspendBackgroundMs(100); // Suspend background task for 100 ms
}

//====================================================================
// Critical Task (Reading + PWM)
//====================================================================

// Real-time task for reading sensors and adjusting PWM.
// - Reads low/high voltages and currents via DataAPI.
// - Selects behavior based on the current mode:
// - MODE_IDLE: PWM is stopped.
// - MODE_POWER: PID regulation based on voltageReferenceP.
// - MODE_EMULATOR: PID regulation initially on voltageReferenceE_test_point,
// then on voltageReferenceE after initial setup.
// - Accumulates measurements for the application task.
//
// Inputs:
// - void
//
// Returns:
// - void
void loop_critical_task() {
    // Read sensors
    measData = shield.sensors.getLatestValue(I1_LOW);
    if (measData != NO_VALUE) lowCurrent1 = measData; // Low-side current 1 [A]

    measData = shield.sensors.getLatestValue(V1_LOW);
    if (measData != NO_VALUE) lowVoltage1 = measData; // Low-side voltage 1 [V]

    measData = shield.sensors.getLatestValue(V2_LOW);
    if (measData != NO_VALUE) lowVoltage2 = measData; // Low-side voltage 2 [V]

    measData = shield.sensors.getLatestValue(I2_LOW);
    if (measData != NO_VALUE) lowCurrent2 = measData; // Low-side current 2 [A]

    measData = shield.sensors.getLatestValue(I_HIGH);
    if (measData != NO_VALUE) highCurrent = measData; // High-side current [A]

    measData = shield.sensors.getLatestValue(V_HIGH);
    if (measData != NO_VALUE) highVoltage = measData; // High-side voltage [V]

    if (mode == MODE_IDLE) {
        // MODE_IDLE: Simply stop PWM, no accumulation
        if (pwmEnable) {
            shield.power.stop(ALL);
        }
        pwmEnable = false;
    } else if (mode == MODE_POWER) {
        // MODE_POWER: PID calculation, application, and accumulation
        dutyCycle = pid.calculateWithReturn(voltageReferenceP, lowVoltage1); // Calculate duty cycle using PID [dimensionless, 0-1]
        shield.power.setDutyCycle(ALL, dutyCycle); // Set PWM duty cycle

        // Accumulate measurements
        sumLowCurrent1 += lowCurrent1;
        sumLowVoltage1 += lowVoltage1;
        sumLowCurrent2 += lowCurrent2;
        sumLowVoltage2 += lowVoltage2;
        sumHighCurrent += highCurrent;
        sumHighVoltage += highVoltage;
        measurementCount++; // Increment measurement count

        if (!pwmEnable) {
            pwmEnable = true;
            shield.power.start(ALL); // Start PWM
        }
    } else { // MODE_EMULATOR
        // EMULATOR: PID regulation on voltageReferenceE_test_point until initial measurement, then on voltageReferenceE
        if (!emulatorSteadyState) {
            dutyCycle = pid.calculateWithReturn(voltageReferenceE_test_point, lowVoltage1);
            shield.power.setDutyCycle(ALL, dutyCycle);
        } else {
            dutyCycle = pid.calculateWithReturn(voltageReferenceE, lowVoltage1);
            shield.power.setDutyCycle(ALL, dutyCycle);
        }

        // Accumulate measurements
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
    setup_routine();
    return 0;
}
