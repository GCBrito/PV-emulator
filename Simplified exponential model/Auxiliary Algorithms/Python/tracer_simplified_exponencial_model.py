import numpy as np
import matplotlib.pyplot as plt

# Equivalente a 'clc; clear; close all;' do MATLAB
# (limpar console/variáveis não é padrão em scripts Python)
plt.close('all')

# Configura a fonte padrão para 'Times New Roman' para todo o gráfico
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

# %% 1) Panel Parameters
Vmp = 15.0  # Maximum Power Voltage (V)
Imp = 4.1   # Maximum Power Current (A)
Voc = 21.0  # Open-Circuit Voltage (V)
Isc = 5.1   # Short-Circuit Current (A)

# %% 2) Simplified exponential equation constant
# This constant 'c' is derived from a simplified exponential model,
c = - (Voc - Vmp) / np.log(1 - Imp / Isc)

# %% 3) Characteristic Points for the linear interpolation
points_V = np.array([
    0,
    0.2 * Vmp,
    0.4 * Vmp,
    0.6 * Vmp,
    0.8 * Vmp,
    0.9 * Vmp,
    0.95 * Vmp,
    Vmp,
    (Vmp + Voc) / 2,
    0.9 * Voc,
    0.95 * Voc,
    Voc
])
points_I = Isc * (1 - np.exp((points_V - Voc) / c))  # Current derived from the simplified exp model

# %% 4) Linear Interpolation

# Initialization of coefficients for the linear interpolation
coeffs = np.zeros((len(points_V) - 1, 2))
for i in range(len(points_V) - 1):
    a = (points_I[i + 1] - points_I[i]) / (points_V[i + 1] - points_V[i])
    b = points_I[i] - a * points_V[i]
    coeffs[i, :] = [a, b]  # Store [a, b] in each row

def piecewise_pv_model(V_in, points_V, coeffs):
    # Garante que a entrada seja um array numpy para iteração
    V_in = np.atleast_1d(V_in)
    I_out = np.zeros_like(V_in, dtype=float)  # Initialize the output current array
    
    for k in range(V_in.size):  # Iterate over each voltage input
        V = V_in[k]
        found = False
        for i in range(len(points_V) - 1):
            # Use a small tolerance for floating-point comparisons
            tol = 1e-9
            if (V >= (points_V[i] - tol) and V <= (points_V[i + 1] + tol)):
                a = coeffs[i, 0]
                b = coeffs[i, 1]
                I_out[k] = a * V + b
                found = True
                break
        
        if not found:
            # If V is outside all segments, typically 0 for V > Voc and V < 0
            if V < points_V[0] or V > points_V[-1]:
                I_out[k] = 0  # Or handle as needed
            else:
                # This case should ideally not be reached if segments cover the range
                I_out[k] = np.nan
    
    # Retorna um escalar se a entrada era um escalar
    return I_out.item() if V_in.ndim == 0 and I_out.size == 1 else I_out

pv_model = lambda V_in: piecewise_pv_model(V_in, points_V, coeffs)

# %% 5) Provided Data: (R, V_test, I_test, V*, I*)
measurements = np.array([
    # ---------------------KC200GT---------------------
    # [3.699029126,30.056282,7.876191,15.578393,4.082288],
    # [4.289340102,30.050781,7.080119,17.230631,4.059626],
    # [4.693298969,29.872391,6.489758,18.552963,4.030620],
    # [4.913265306,29.879427,6.065750,19.599977,3.978944],
    # [5.397849462,29.889881,5.765186,20.419359,3.938504],
    # [5.911111111,29.917889,5.258528,21.612139,3.798665],
    # [6.507204611,29.941908,4.745583,22.879311,3.626210],
    # [7.993311037,29.952457,3.937661,24.182100,3.179068],
    # [9.636015326,29.963566,3.255278,25.403482,2.759865],
    # [12.20283019,29.978542,2.608978,26.099630,2.271403],
    # [15.66272189,29.987871,2.071892,26.694136,1.844324],
    # [40.19117647,29.994781,0.951167,27.493385,0.871845],
    # [410.2941176,30.012428,0.284112,27.971466, 0.264791],
    # ---------------------KC85TS---------------------
    # [2.528599606,23.775753,9.540161,13.240793,5.312946],
    # [2.694552529,23.927464,8.889205,14.250062,5.293989],
    # [2.948207171,23.923685,8.276764,15.192984,5.256245],
    # [3.116,23.798820,7.758160,15.953911,5.200803],
    # [3.481404959,23.811413,6.968660,17.226923,5.041640],
    # [4.040540541,23.839123,6.009305,18.291187,4.610795],
    # [4.905370844,23.851524,4.968151,19.491007,4.059877],
    # [6.158385093,23.908186,3.995351,20.092628,3.357725],
    # [8.023622047,23.868225,3.134047,20.618151,2.707292],
    # [15.85606061,23.903660,1.681039,21.106844,1.484351],
    # [405.6603774,23.868862,0.325685,21.582323,0.294486],
    # ---------------Uni-Solar ES-62T-----------------
    [2.159645233,22.119232,10.492447,10.140304,4.810140],
    [2.591224018,21.972612,8.881744,11.602366,4.689895],
    [2.891203704,22.107285,7.773019,12.875836,4.527201],
    [3.356968215,21.984802,6.669855,14.143713,4.290987],
    [3.854497354,22.010578,6.045714,14.959622,4.109006],
    [5.061919505,22.015034,4.473690,16.689295,3.391443],
    [6.518518519,22.034122,3.533989,17.921898,2.874441],
    [9.083373964,22.046028,2.588672,18.896595,2.218861],
    [18.08256881,22.052059,1.446311,19.905468,1.305524],
    [58.02857143,22.056940,0.659405,20.491583,0.612608],
])

# Extract columns for easier use
if measurements.size > 0:
    R_data = measurements[:, 0]
    V_magenta = measurements[:, 1]  # V_test -> Magenta Points (Test Point)
    I_magenta = measurements[:, 2]  # I_test -> Magenta Points (Test Point)
    V_black = measurements[:, 3]    # V* -> Black Points (Intersection)
    I_black = measurements[:, 4]    # I* -> Black Points (Intersection)
else:
    R_data, V_magenta, I_magenta, V_black, I_black = [], [], [], [], []

# %% 6) Plotting (Figure Generation)
fig, ax = plt.subplots(figsize=(12, 8))

max_V_plot = np.max([Voc, np.max(V_magenta) if len(V_magenta)>0 else 0, np.max(V_black) if len(V_black)>0 else 0])
max_I_plot = np.max([Isc, np.max(I_magenta) if len(I_magenta)>0 else 0, np.max(I_black) if len(I_black)>0 else 0]) * 1.1

# 1. Plot "Load Lines"
Rs_to_plot = np.unique(R_data)
if Rs_to_plot.size == 0:
    Rs_to_plot = np.sort(np.array([0.1, 0.5] + list(range(1, 16)) + list(range(20, 51, 5)) + \
                                  list(range(60, 201, 10)) + [300, 400, 500, 1000]))

V_line_range = np.linspace(0, max_V_plot * 1.1, 100)
for r_val in Rs_to_plot:
    if np.isinf(r_val): # Open Circuit (I = 0)
        ax.plot(V_line_range, np.zeros_like(V_line_range), '--', linewidth=0.8, color=[0.5, 0.5, 0.5])
    elif r_val == 0: # Short Circuit (V = 0)
        ax.plot(np.zeros_like(V_line_range), np.linspace(0, max_I_plot, 100), '--', linewidth=0.8, color=[0.5, 0.5, 0.5])
    else:
        I_line = V_line_range / r_val
        ax.plot(V_line_range, I_line, '--', linewidth=0.8, color=[0.5, 0.5, 0.5])

# 2. Plot the I-V Model (PV Curve)
V_plot = np.linspace(0, Voc, 500)
I_plot = pv_model(V_plot)
h_model, = ax.plot(V_plot, I_plot, 'r-', linewidth=2, label='I-V Curve - Simplified exponencial model')

# 3. Plot "Test Points" (Magenta Points)
h_test_point = None
if len(V_magenta) > 0:
    h_test_point = ax.scatter(V_magenta, I_magenta, s=40, c='m', marker='o', label='Test Point')

# 4. Plot "Intersections" (Black Points)
h_intersection = None
if len(V_black) > 0 and len(I_black) > 0:
    h_intersection = ax.scatter(V_black, I_black, s=70, c='k', marker='s', label='Intersection')

# Item para a legenda da "Load Line", plotado fora da área visível
h_load_line_legend, = ax.plot([], [], '--', color=[0.5, 0.5, 0.5], label='Load Line')

# %% Plot Configurations
ax.set_xlabel('Voltage (V)', fontsize=30)
ax.set_ylabel('Current (A)', fontsize=30)

# Axis limits
ax.set_xlim(0, 1.2 * Voc)
ax.set_ylim(0, 1.2 * np.max(I_magenta) if len(I_magenta) > 0 else 1.2 * Isc)
ax.tick_params(axis='both', which='major', labelsize=30)
ax.grid(False) # Equivalente a grid off

# Coleta handles e labels para a legenda dinamicamente
handles, labels = ax.get_legend_handles_labels()
leg = ax.legend(handles=handles, labels=labels, loc='best', fontsize=20, frameon=True)

plt.tight_layout()
plt.show()