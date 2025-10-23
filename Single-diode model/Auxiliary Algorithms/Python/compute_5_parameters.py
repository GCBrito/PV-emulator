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
    a_sc = min(C*Rs*Isc, cap)
    a_oc = min(C*Voc, cap)
    a_mp = min(C*(Rs*Imp+Vmp), cap)
    a_eq5 = min(C*Is0, cap)

    F = np.zeros(5)
    F[0] = Isc - ( Iph - Is0*(np.exp(a_sc)-1) - (Rs*Isc)/Rp )
    F[1] = 0   - ( Iph - Is0*(np.exp(a_oc)-1) - (Voc)/Rp )
    F[2] = Imp - ( Iph - Is0*(np.exp(a_mp)-1) - (Rs*Imp+Vmp)/Rp )
    F[3] = Rs + (q*Is0*Rp*(Rs-Rp)/(A*k*Tref))*np.exp(a_eq5)
    term_coeff = 1 + q*(Vmp - Rs*Imp)/(A*k*Tref)
    F[4] = Iph - 2*Vmp/Rp + Is0 - Is0 * term_coeff * np.exp(a_mp)
      
    return F

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

