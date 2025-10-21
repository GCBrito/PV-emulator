import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
import warnings

# --- Auxiliary Functions ---

def residuals_2_20(x, Voc, Isc, Vmp, Imp, q, k, Tref):
    """
    Calculates the residuals for the five-point model equations.
    This function is used by least_squares to find the model parameters.
    """
    Iph, Is0, A, Rs, Rp = x

    C = q / (A * k * Tref)
    cap = 700  # Safety cap to prevent overflow in exp()

    # Arguments for the exponential function, capped for stability
    a_sc = min(C * Rs * Isc, cap)
    a_oc = min(C * Voc, cap)
    a_mp = min(C * (Rs * Imp + Vmp), cap)

    F = np.zeros(5)
    F[0] = Isc - (Iph - Is0 * (np.exp(a_sc) - 1) - (Rs * Isc) / Rp)
    F[1] = 0 - (Iph - Is0 * (np.exp(a_oc) - 1) - Voc / Rp)
    F[2] = Imp - (Iph - Is0 * (np.exp(a_mp) - 1) - (Rs * Imp + Vmp) / Rp)
    F[3] = Rs + (q * Is0 * Rp * (Rs - Rp) / (A * k * Tref)) * np.exp(a_sc)
    term = 1 + q * (Vmp - Rs * Imp) / (A * k * Tref)
    F[4] = Iph - 2 * Vmp / Rp - (Is0 * term * (np.exp(a_mp) - 1))
    
    return F

def solve_I_V_2_11(V, Iph_ref, Is0_ref, A, Rs, Rp,
                   q, k, S, Sref, alpha, T, Tref, E_G0, k1, k2):
    """
    Solves for the cell current (I) at a given voltage (V) using Newton-Raphson.
    """
    # 1. Photocurrent with T and S dependency
    Iph = (S / Sref) * (Iph_ref + alpha * (T - Tref))
    
    # 2. Temperature-dependent Band Gap Energy (Varshni's equation)
    Eg_ref_J = (E_G0 - (k1 * Tref**2) / (Tref + k2)) * q
    Eg_J = (E_G0 - (k1 * T**2) / (T + k2)) * q
    
    # 3. Saturation Current with T and Eg(T) dependency
    exponent_term = (Eg_ref_J / (A * k * Tref)) - (Eg_J / (A * k * T))
    Is = Is0_ref * (T / Tref)**3 * np.exp(exponent_term)
    
    # 4. Iterative solver (Newton-Raphson)
    I = Iph  # Initial guess
    for _ in range(30): # Iterate up to 30 times
        V_diode = V + I * Rs
        arg_exp = min(q * V_diode / (A * k * T), 700) # Capped for stability
        expo = np.exp(arg_exp)
        
        f = Iph - Is * (expo - 1) - V_diode / Rp - I
        df = -Is * expo * (q * Rs / (A * k * T)) - Rs / Rp - 1
        
        # Avoid division by zero
        if abs(df) < 1e-9:
            break
            
        dI = -f / df
        I = I + dI
        
        if abs(dI) < 1e-6:
            break
            
    return I

def piecewise_pv_model(V_in, points_V, points_I):
    """
    Calculates the current from a piecewise linear model using interpolation.
    """
    # np.interp performs piecewise linear interpolation.
    # It requires the x-coordinates (points_V) and y-coordinates (points_I).
    # `left` and `right` define values for inputs outside the interpolation range.
    I_out = np.interp(V_in, points_V, points_I, left=0.0, right=0.0)
    return I_out

# --- Main Script ---

## Module Parameters
ns = 60       # Number of series cells
n_p = 1       # Number of parallel branches 

Vmp_mod_ref = 30.1 # Voltage at maximum power point of the module (V)
Imp_mod_ref = 8.30 # Current at maximum power point of the module (A)
Voc_mod_ref = 37.2 # Open-circuit voltage of the module (V)
Isc_mod_ref = 8.87 # Short-circuit current of the module (A)

Tref = 25 + 273.15  # Reference temperature (K)
Sref = 1000         # Reference irradiance (W/m²)
alpha = 0.00065     # temperature coefficient of Isc (%/K)

## Operating Conditions
T = 44.5 + 273.15   # Current temperature (K)
S = 765           # Current irradiance (W/m²)

## PV cell parameters 
Vmp_cell_ref = Vmp_mod_ref / ns
Imp_cell_ref = Imp_mod_ref / n_p 
Voc_cell_ref = Voc_mod_ref / ns
Isc_cell_ref = Isc_mod_ref / n_p 

## Physical Constants
q = 1.60217662e-19  # Elementary charge (C)
k = 1.38064852e-23  # Boltzmann constant (J/K)
E_G0 = 1.166          # Band gap energy at 0K (eV)
k1 = 4.73e-4          # Coefficient k1 (eV/K)
k2 = 636              # Coefficient k2 (K)

## 4) Parameter Estimation via least_squares (more robust than fsolve)
# Initial guess based on MATLAB's x0 for better convergence
x0 = np.array([Isc_cell_ref, 1e-9, 3.0, 0.01, 2.0])

# Define physical bounds for the parameters ([Iph, Is0, A, Rs, Rp])
lb = [0, 0, 0.5, 0, 1e-3]  # Lower bounds
ub = [Isc_cell_ref * 1.2, 1e-6, 5.0, 1.0, 10000] # Upper bounds

result = least_squares(
    residuals_2_20, 
    x0, 
    args=(Voc_cell_ref, Isc_cell_ref, Vmp_cell_ref, Imp_cell_ref, q, k, Tref), 
    bounds=(lb, ub),
    method='trf', # Trust Region Reflective algorithm, good for bounds
    ftol=1e-10,
    xtol=1e-10,
    gtol=1e-10
)

if not result.success:
    warnings.warn(f'least_squares did not converge. Status: {result.status}, Message: {result.message}')

Iph_ref, Is0_ref, A, Rs, Rp = result.x
print(f'Iph_ref = {Iph_ref:.6f} A')
print(f'Is0_ref = {Is0_ref:.2e} A')
print(f'A       = {A:.6f}')
print(f'Rs      = {Rs:.6f} Ohms')
print(f'Rp      = {Rp:.6f} Ohms')

## 5) Characteristic Points (used to build the I-V model)
V_cell = np.array([
    0.0,
    0.2 * Vmp_mod_ref / ns,
    0.4 * Vmp_mod_ref / ns,
    0.6 * Vmp_mod_ref / ns,
    0.8 * Vmp_mod_ref / ns,
    0.9 * Vmp_mod_ref / ns,
    Vmp_mod_ref / ns,
    (Vmp_mod_ref / ns + Voc_mod_ref / ns) / 2.0,
    0.9 * Voc_mod_ref / ns,
    0.95 * Voc_mod_ref / ns,
    Voc_mod_ref / ns
])

# Use a list comprehension, which is a Pythonic way to apply a function to each element
I_cell = np.array([solve_I_V_2_11(v, Iph_ref, Is0_ref, A, Rs, Rp,
                                q, k, S, Sref, alpha, T, Tref, E_G0, k1, k2) for v in V_cell])

points_V = V_cell * ns
points_I = I_cell * n_p # Using the corrected variable n_p

## 6) Piecewise PV Model (Coefficient calculation is no longer needed)
# The model function now uses points_V and points_I directly
modele_pv = lambda V: piecewise_pv_model(V, points_V, points_I)


## 7) Provided Data: (R, V_test, I_test, V*, I*)
mesures = np.array([
    
    # Data from CS6P-250P under S = 556, T = 33
    # [3.312101911, 35.098381, 10.830381, 15.96633, 4.926764],
    # [4.470212766, 35.093315, 8.052373, 21.384974, 4.906911],
    # [5.072649573, 34.928596, 7.081912, 24.119793, 4.890384],
    # [5.428571429, 34.964813, 6.637497, 25.479982, 4.836957],
    # [5.846827133, 34.940865, 6.137167, 27.134899, 4.766093],
    # [6.560465116, 34.992397, 5.522508, 28.611193, 4.515424],
    # [7.335802469, 35.004902, 4.953645, 30.103748, 4.260068],
    # [11.3297491, 35.064438, 3.307633, 31.929382, 3.011903],
    # [20.85, 35.042049, 1.906419, 33.620697, 1.829092],
    # [62.91512915, 35.042969, 0.812353, 34.34491, 0.796171],
    
    
    # Data from CS6P-250P under S = 765, T = 44.5 
     [3.424196018,34.855976,10.283596,22.807545,6.728934],
     [3.684782609,34.869011,9.658541,24.17536,6.696453],
     [4.031847134,34.885811,8.816074,25.775291,6.513734],
     [4.44278607,34.931274,8.078631,27.240286,6.29992],
     [5.227356747,34.942135,6.905494,28.715485,5.674942],
     [5.995125914,34.961349,6.021228,29.939571,5.156351],
     [7.701643489,34.964321,4.749564,30.825687,4.187371],
     [13.45162653,34.977104,2.806052,32.122658,2.577053],
     [68.56557377,35.020119,0.753771,33.615498,0.723538],
])

# Use Python's 0-based indexing
R_data = mesures[:, 0]
V_violet = mesures[:, 1]  # Purple Points (Test Point)
I_violet = mesures[:, 2]  # Purple Points (Test Point)
V_noir = mesures[:, 3]    # Black Points (Intersection)
I_noir = mesures[:, 4]    # Black Points (Intersection)

## 8) Plotting (Figure Generation)
plt.style.use('default') # Use a standard style
fig, ax = plt.subplots(figsize=(12, 8))

# 1. Plot "Load Lines"
Rs_to_plot = np.sort(np.unique(R_data))
V_line_max = np.max(V_violet) * 1.25 if len(V_violet) > 0 else Voc_mod_ref * 1.1
V_line = np.linspace(0, V_line_max, 100)

for r_val in Rs_to_plot:
    if r_val > 0 and np.isfinite(r_val):
        I_line = V_line / r_val
        # Corrected plot call to avoid warning
        ax.plot(V_line, I_line, linestyle='--', linewidth=0.8, color='0.5')

# 2. Plot the I-V Model
V_plot = np.linspace(0, Voc_mod_ref, 500)
I_plot = modele_pv(V_plot)
h_model, = ax.plot(V_plot, I_plot, 'r-', linewidth=2, label='I-V Curve - Single diode model')

# 3. Plot "Test Points" (Purple Points)
h_test_point = None
if V_violet.size > 0:
    h_test_point = ax.scatter(V_violet, I_violet, s=40, c='m', marker='o', label='Test Point', zorder=5)

# 4. Plot "Intersections" (Black Points)
h_intersection = None
if V_noir.size > 0 and I_noir.size > 0:
    h_intersection = ax.scatter(V_noir, I_noir, s=70, c='k', marker='s', label='Intersection', zorder=5)

## 9) Plot Configurations
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Times New Roman'
plt.rcParams['axes.labelweight'] = 'normal'

ax.set_xlabel('Voltage (V)', fontsize=22)
ax.set_ylabel('Current (A)', fontsize=22)
ax.set_xlim(0, 1.01 * Voc_mod_ref)
# Use max of measured current or calculated short-circuit current for robust y-limit
ylim_max = max(np.max(I_violet) if I_violet.size > 0 else 0, Isc_mod_ref)
ax.set_ylim(0, 1.1 * ylim_max)
ax.tick_params(axis='both', which='major', labelsize=20)
ax.grid(False)

# Create the legend
handles = [h_model]
if h_test_point:
    handles.append(h_test_point)
if h_intersection:
    handles.append(h_intersection)

leg = ax.legend(handles=handles, loc='upper right', fontsize=16) # 'NorthWest' in MATLAB is often 'upper left' or 'upper right' in Matplotlib
leg.get_frame().set_edgecolor('k') # Add box around legend
leg.get_frame().set_linewidth(1.0)

plt.tight_layout()
plt.show()

