"""Plot a piecewise approximation of the simplified exponential PV model.

Standalone equivalent of tracer_simplified_exponencial_model.m.
Figures are displayed only; this script does not create files or folders.
"""
import argparse

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager


def configure_plots(font_size=30):
    """Use Times New Roman where installed; otherwise use a serif fallback."""
    fonts = {font.name for font in font_manager.fontManager.ttflist}
    family = "Times New Roman" if "Times New Roman" in fonts else "DejaVu Serif"
    plt.rcParams.update({
        "font.family": family, "font.size": font_size,
        "mathtext.fontset": "stix", "figure.facecolor": "white",
        "axes.labelsize": font_size, "xtick.labelsize": font_size,
        "ytick.labelsize": font_size, "legend.fontsize": 20,
    })


def finish(figures, args):
    """Display interactive plot windows without writing image files."""
    if not args.no_show:
        plt.show()
    else:
        for _, figure in figures:
            plt.close(figure)


def plot_arguments(description):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--no-show", action="store_true", help="Do not open plot windows.")
    return parser


# Editable panel parameters, copied from MATLAB.
VMP = 34.0  # V
IMP = 4.4   # A
VOC = 43.4  # V
ISC = 4.8   # A

# Columns: R_load [ohm], V_test [V], I_test [A], V_emulator [V], I_emulator [A].
# The active rows are kept exactly as provided (including their column order).
MEASUREMENTS = np.array([
    [5.109090909, 45.275391, 8.406238, 25.567003, 4.747],
    [5.365269461, 45.267052, 7.877808, 27.200571, 4.733705],
    [5.730382294, 45.318977, 7.373886, 28.811605, 4.687959],
    [6.126530612, 45.318756, 6.928646, 30.372927, 4.643624],
    [6.752642706, 45.338726, 6.330338, 32.341831, 4.51567],
    [7.41978022, 45.339798, 5.804495, 34.1399, 4.370661],
    [8.985037406, 45.374962, 4.861022, 36.391037, 3.898574],
    [11.06051873, 45.385956, 3.990404, 38.724979, 3.404761],
    [14.33212996, 45.385956, 3.132179, 40.008648, 2.759869],
    [20.28712871, 45.413887, 2.293783, 41.242683, 2.083102],
    [54.08974359, 45.4189, 1.05701, 42.378609, 0.986255],
    [428.4, 45.421684, 0.414125, 42.994045, 0.391991],
], dtype=float)

# Inactive datasets preserved from the original MATLAB comments.
# Select a dataset only with matching panel parameters.
ALTERNATIVE_MEASUREMENTS = {
    'KC200GT': np.array([
        [3.699029126, 30.056282, 7.876191, 15.578393, 4.082288],
        [4.289340102, 30.050781, 7.080119, 17.230631, 4.059626],
        [4.693298969, 29.872391, 6.489758, 18.552963, 4.030620],
        [4.913265306, 29.879427, 6.065750, 19.599977, 3.978944],
        [5.397849462, 29.889881, 5.765186, 20.419359, 3.938504],
        [5.911111111, 29.917889, 5.258528, 21.612139, 3.798665],
        [6.507204611, 29.941908, 4.745583, 22.879311, 3.626210],
        [7.993311037, 29.952457, 3.937661, 24.182100, 3.179068],
        [9.636015326, 29.963566, 3.255278, 25.403482, 2.759865],
        [12.20283019, 29.978542, 2.608978, 26.099630, 2.271403],
        [15.66272189, 29.987871, 2.071892, 26.694136, 1.844324],
        [40.19117647, 29.994781, 0.951167, 27.493385, 0.871845],
        [410.2941176, 30.012428, 0.284112, 27.971466, 0.264791],
    ], dtype=float),
    'KC85TS': np.array([
        [2.528599606, 23.775753, 9.540161, 13.240793, 5.312946],
        [2.694552529, 23.927464, 8.889205, 14.250062, 5.293989],
        [2.948207171, 23.923685, 8.276764, 15.192984, 5.256245],
        [3.116, 23.798820, 7.758160, 15.953911, 5.200803],
        [3.481404959, 23.811413, 6.968660, 17.226923, 5.041640],
        [4.040540541, 23.839123, 6.009305, 18.291187, 4.610795],
        [4.905370844, 23.851524, 4.968151, 19.491007, 4.059877],
        [6.158385093, 23.908186, 3.995351, 20.092628, 3.357725],
        [8.023622047, 23.868225, 3.134047, 20.618151, 2.707292],
        [15.85606061, 23.903660, 1.681039, 21.106844, 1.484351],
        [405.6603774, 23.868862, 0.325685, 21.582323, 0.294486],
    ], dtype=float),
    'Uni-Solar ES-62T': np.array([
        [2.159645233, 22.119232, 10.492447, 10.140304, 4.810140],
        [2.591224018, 21.972612, 8.881744, 11.602366, 4.689895],
        [2.891203704, 22.107285, 7.773019, 12.875836, 4.527201],
        [3.356968215, 21.984802, 6.669855, 14.143713, 4.290987],
        [3.854497354, 22.010578, 6.045714, 14.959622, 4.109006],
        [5.061919505, 22.015034, 4.473690, 16.689295, 3.391443],
        [6.518518519, 22.034122, 3.533989, 17.921898, 2.874441],
        [9.083373964, 22.046028, 2.588672, 18.896595, 2.218861],
        [18.08256881, 22.052059, 1.446311, 19.905468, 1.305524],
        [58.02857143, 22.056940, 0.659405, 20.491583, 0.612608],
    ], dtype=float),
}


def characteristic_points(Vmp=VMP, Imp=IMP, Voc=VOC, Isc=ISC):
    """Original exponential constant, voltage samples and sampled currents."""
    if not (Voc > Vmp > 0 and Isc > Imp > 0):
        raise ValueError("Require Voc > Vmp > 0 and Isc > Imp > 0.")
    c = -(Voc - Vmp) / np.log(1 - Imp / Isc)
    points_V = np.array([0, .2 * Vmp, .4 * Vmp, .6 * Vmp, .8 * Vmp,
                         .9 * Vmp, .95 * Vmp, Vmp, (Vmp + Voc) / 2,
                         .9 * Voc, .95 * Voc, Voc])
    points_I = Isc * (1 - np.exp((points_V - Voc) / c))
    return c, points_V, points_I


def piecewise_coefficients(points_V, points_I):
    """Store [a, b] for I = a*V + b, preserving the original segment order."""
    points_V, points_I = np.asarray(points_V), np.asarray(points_I)
    if len(points_V) != len(points_I) or len(points_V) < 2 or np.any(np.diff(points_V) == 0):
        raise ValueError("Adjacent voltage samples must be distinct and arrays equally sized.")
    a = np.diff(points_I) / np.diff(points_V)
    b = points_I[:-1] - a * points_V[:-1]
    return np.column_stack((a, b))


def piecewise_pv_model(V_in, points_V, coeffs):
    """Original first-matching-segment rule, including 1e-9 tolerance.

    Inputs may be scalars or arrays. Current is zero outside the endpoints,
    and NaN if an interior voltage cannot be assigned to any segment.
    """
    values = np.asarray(V_in, dtype=float)
    result = np.zeros_like(values)
    found = np.zeros(values.shape, dtype=bool)
    for index, (a, b) in enumerate(coeffs):
        mask = ((~found) & (values >= points_V[index] - 1e-9)
                & (values <= points_V[index + 1] + 1e-9))
        result = np.where(mask, a * values + b, result)
        found |= mask
    interior = (values >= points_V[0]) & (values <= points_V[-1])
    result = np.where((~found) & interior, np.nan, result)
    return float(result) if result.ndim == 0 else result


def main(argv=None):
    args = plot_arguments(__doc__).parse_args(argv)
    configure_plots()
    c, points_V, points_I = characteristic_points(VMP, IMP, VOC, ISC)
    coeffs = piecewise_coefficients(points_V, points_I)
    measurements = np.asarray(MEASUREMENTS, dtype=float)
    if measurements.size:
        if measurements.ndim != 2 or measurements.shape[1] != 5:
            raise ValueError("MEASUREMENTS must have five columns.")
        R_data, V_test, I_test, V_emu, I_emu = measurements.T
        max_voltage = max(VOC, V_test.max(), V_emu.max())
        max_current = max(ISC, I_test.max(), I_emu.max()) * 1.1
        upper_current_limit = 1.2 * I_test.max()
        resistance_values = np.unique(R_data)
    else:
        # The original final ylim uses max([]) for empty measurements.
        # Use the model Isc as a fallback to keep that branch usable.
        R_data = V_test = I_test = V_emu = I_emu = np.array([])
        max_voltage, max_current, upper_current_limit = VOC, 1.1 * ISC, 1.2 * ISC
        resistance_values = np.sort(np.r_[.1, .5, np.arange(1, 16),
                                         np.arange(20, 51, 5), np.arange(60, 201, 10),
                                         300, 400, 500, 1000])
    fig, ax = plt.subplots(figsize=(12, 7), layout="constrained")
    V_line = np.linspace(0, max_voltage * 1.1, 100)
    for resistance in resistance_values:
        if np.isinf(resistance):
            x, y = V_line, np.zeros_like(V_line)
        elif resistance == 0:
            x, y = np.zeros_like(V_line), np.linspace(0, max_current, 100)
        else:
            x, y = V_line, V_line / resistance
        ax.plot(x, y, "--", color=(.5, .5, .5), linewidth=.8)
    load_handle, = ax.plot(np.nan, np.nan, "--", color=(.5, .5, .5))
    V_plot = np.linspace(0, VOC, 500)
    model_handle, = ax.plot(V_plot, piecewise_pv_model(V_plot, points_V, coeffs),
                            "r-", linewidth=2)
    handles, labels = [model_handle, load_handle], ["Simplified exp. model I-V curve", "Load line"]
    if measurements.size:
        handles.append(ax.scatter(V_test, I_test, s=40, color="m"))
        labels.append("Test point")
        handles.append(ax.scatter(V_emu, I_emu, s=70, color="k", marker="s"))
        labels.append("Intersection")
    ax.set(xlabel="Voltage (V)", ylabel="Current (A)",
           xlim=(0, 1.2 * VOC), ylim=(0, upper_current_limit))
    ax.grid(False)
    ax.legend(handles, labels, loc="upper left", fontsize=20, frameon=True)
    finish([("simplified_exponential_model", fig)], args)
    return c, points_V, points_I, coeffs


if __name__ == "__main__":
    main()
