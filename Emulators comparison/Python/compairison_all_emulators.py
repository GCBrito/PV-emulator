"""Compare three PV emulators against three measured PV-panel I-V curves.

Standalone equivalent of compairison_all_emulators.m; original filename retained.
Compute RMSE, nRMSE, filtered MAPE and characteristic-point errors; display
three I-V plots, two summary bar charts and a terminal summary. No files saved.
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

from scipy.interpolate import PchipInterpolator, interp1d
from matplotlib.patches import Polygon

# Editable plot and interpolation settings, copied from MATLAB.
FONT_AXIS_LABEL = 30
FONT_TICKS = 26
FONT_LEGEND = 15
FONT_CALLOUT = 15
LINE_WIDTH_EMU = 2.2
MARKER_SIZE_REF = 6
N_COMMON_GRID = 1000

PANELS = [
    {
        'name': 'ME Solar MESM-50W',
        'condition': 'G = 1000 W/m^2, T = 25 degC',
        'V_panel': np.array([0.1345090408816909, 3.9418622020801637, 7.1102974135229875, 10.316096247185497, 13.219249737750683, 14.493349274149194, 15.386339859898216, 16.114930504317343, 16.794948437178412, 17.542220892707384, 18.26333881312798, 18.943356748216104, 19.39172022153349, 19.65700194241112, 19.90733821529171, 20.228665369498877, 20.51636526682619, 20.882528769478625, 21.1702286645789, 21.446719472567864, 21.71200119344551], dtype=float),
        'I_panel': np.array([3.0107092341738317, 2.9969690844678114, 2.9832289352435426, 2.970297029807906, 2.9614063447606873, 2.9533239039838524, 2.9533239035021017, 2.9525156597134687, 2.932309557771381, 2.880581935547086, 2.7868256213795988, 2.643766417317217, 2.503940189854619, 2.3099616083200765, 2.122448979503352, 1.8508789653549917, 1.5663770457709956, 1.1840775911493435, 0.8163265305040968, 0.4598908867537168, 0.030713275530074746], dtype=float),
        "emu": [
            {
                "name": 'Open-source PV Emulator - Single Diode Model', "shortName": 'Single diode',
                "color": (0, 0, 1),
                "lineStyle": '-',
                "V": np.array([0, 5.922083, 12.551384, 13.841137, 15.125275, 15.540461, 15.992638, 16.356771, 16.715565, 17.007088, 17.418226, 17.789234, 18.138689, 18.464882, 18.834526, 19.264559, 19.65144, 19.971132, 20.552086, 21.464041, 21.969181, 22.3], dtype=float),
                "I": np.array([3.03, 3.026813, 3.019652, 3.0118, 2.995558, 2.984576, 2.96906, 2.951598, 2.929032, 2.905429, 2.862687, 2.811599, 2.749611, 2.676852, 2.569205, 2.419445, 2.236233, 2.070245, 1.553807, 0.743124, 0.29408, 0], dtype=float),
            },
            {
                "name": 'Open-source PV Emulator - Simplified Exp. Model', "shortName": 'Simplified exp.',
                "color": (1, 0, 0),
                "lineStyle": '-',
                "V": np.array([0, 8.08962, 11.45164, 14.526568, 16.046844, 16.94203, 17.753181, 18.2904, 19.387537, 20.072838, 20.594084, 21.106194, 21.724693, 22.3], dtype=float),
                "I": np.array([3.03, 3.028738, 3.021295, 2.99427, 2.949897, 2.878464, 2.813736, 2.679999, 2.389157, 2.202046, 1.848619, 1.501389, 0.7471, 0], dtype=float),
            },
            {
                "name": 'Commercial PV emulator', "shortName": 'Commercial',
                "color": (1, 0, 1),
                "lineStyle": '-',
                "V": np.array([1.1, 2.7, 4.6, 6.5, 8.5702, 9.4, 13.06, 13.8, 16.2, 17.6774, 18.83, 19.4, 20.1742, 20.939, 21.12, 21.485, 21.7584, 21.93, 22.1151], dtype=float),
                "I": np.array([3.0261, 3.0203, 3.0129, 3.0024, 2.991, 2.985, 2.95, 2.9382, 2.88, 2.818, 2.59, 2.205, 1.61, 1.033, 0.89, 0.6146, 0.4062, 0.2774, 0.1343], dtype=float),
            },
        ],
    },
    {
        'name': 'Renogy RNG-50DB-H',
        'condition': 'G = 1000 W/m^2, T = 25 degC',
        'V_panel': np.array([0.061238999661659355, 3.0378744826113886, 6.06593628879359, 9.019207316748718, 11.551928051656516, 12.757528558290698, 13.855704699997231, 14.771648303222042, 15.865224726830606, 16.729902361354295, 17.664819810269254, 18.450378783493765, 19.100584026063657, 19.605999152099663, 20.04147065780234, 20.416281613412224, 20.80544927623029, 21.20426937792104, 21.462738787922408, 21.777746419471736, 22.02257696007481, 22.272296443566514, 22.475133634590215, 22.579352135807298], dtype=float),
        'I_panel': np.array([2.9194630872527965, 2.9194630872527965, 2.915548098563706, 2.915548098563706, 2.913870246173136, 2.9149888143223963, 2.9077181208521656, 2.8987695749913653, 2.879753914787145, 2.848993288182294, 2.794742729109542, 2.7125279639717688, 2.6045861294008548, 2.4854586130039786, 2.342281879231173, 2.1851230434259357, 1.9737136465477885, 1.713087248101959, 1.4787472036639384, 1.170022371466326, 0.8747203580599141, 0.5447427290262028, 0.23937360194308832, 0.012304250975298636], dtype=float),
        "emu": [
            {
                "name": 'Open-source PV Emulator - Single Diode Model', "shortName": 'Single diode',
                "color": (0, 0, 1),
                "lineStyle": '-',
                "V": np.array([0, 5.894531, 12.125871, 13.816704, 16.145973, 16.715784, 17.231913, 17.664825, 18.180428, 18.727219, 19.188766, 19.602093, 19.970377, 20.420042, 20.759869, 21.242069, 21.980999, 22.318256, 22.6], dtype=float),
                "I": np.array([2.92, 2.915734, 2.909352, 2.902917, 2.874516, 2.856222, 2.831469, 2.801653, 2.751547, 2.673501, 2.579515, 2.459445, 2.335244, 2.107949, 1.91376, 1.412265, 0.643768, 0.293014, 0], dtype=float),
            },
            {
                "name": 'Open-source PV Emulator - Simplified Exp. Model', "shortName": 'Simplified exp.',
                "color": (1, 0, 0),
                "lineStyle": '-',
                "V": np.array([0, 6.817072, 11.080575, 15.057276, 16.803661, 17.610468, 18.552885, 19.459518, 20.528112, 21.02478, 21.487143, 22.002245, 22.6], dtype=float),
                "I": np.array([2.92, 2.919855, 2.918193, 2.894286, 2.843843, 2.780185, 2.695215, 2.44177, 2.114295, 1.793764, 1.483582, 0.796883, 0], dtype=float),
            },
            {
                "name": 'Commercial PV emulator', "shortName": 'Commercial',
                "color": (1, 0, 1),
                "lineStyle": '-',
                "V": np.array([0, 1, 1.2329, 3.09, 5.12, 8.5471, 10.24, 12.12, 14.44, 15.2639, 16.61, 17.83, 18.35, 19.41, 19.69, 20.02, 20.6976, 21.1701, 21.5176, 21.8835, 22.0688, 22.24, 22.5], dtype=float),
                "I": np.array([2.92, 2.917, 2.9162, 2.9106, 2.9031, 2.8875, 2.8768, 2.8623, 2.8358, 2.8236, 2.795, 2.756, 2.71, 2.51, 2.289, 2.0328, 1.5, 1.127, 0.852, 0.5613, 0.4157, 0.2815, 0.0755], dtype=float),
            },
        ],
    },
    {
        'name': 'Shell Solar SQ150-PC',
        'condition': 'G = 1000 W/m^2, T = 25 degC',
        'V_panel': np.array([0.07780604683747905, 4.2189436898142905, 7.660223341663279, 10.863995324097095, 14.176868703729975, 17.155893827390162, 20.141363211697026, 23.036948208895787, 24.995167557029617, 26.42690527335035, 27.81375668458678, 29.001602344647385, 30.36284434236471, 31.512299834342564, 32.59122368036738, 33.76659988606946, 34.743098058886865, 35.65550884640268, 36.446045441351224, 37.16601846173125, 37.82826853208371, 38.47780170128618, 39.10175339833374, 39.68073097512157, 40.29205756197019, 40.91626468796069, 41.54063145321372, 42.171478402636254, 42.8474591077566, 43.38200934569101], dtype=float),
        'I_panel': np.array([4.804197165221235, 4.795540075268057, 4.798166796403276, 4.789358636164255, 4.791964661586494, 4.79244477341078, 4.786709596059194, 4.780959931956091, 4.775059198186533, 4.773217832998629, 4.756864468979553, 4.734262703388705, 4.699256226304568, 4.642458468761035, 4.56182008609613, 4.442863222890885, 4.284504230031027, 4.096089321638853, 3.878645240109937, 3.6456489653261293, 3.397102567352327, 3.1164364242871767, 2.809864790583181, 2.5115743462965194, 2.1490537237230516, 1.7761746301469223, 1.3618533741500927, 0.9319923422019132, 0.4524079577975444, 0.02045929702601157], dtype=float),
        "emu": [
            {
                "name": 'Open-source PV Emulator - Single Diode Model', "shortName": 'Single diode',
                "color": (0, 0, 1),
                "lineStyle": '-',
                "V": np.array([0, 23.428795, 25.957676, 27.723593, 28.600252, 29.371901, 30.186039, 31.42337, 32.285755, 33.184868, 33.870052, 34.74757, 35.65852, 36.230217, 36.874702, 37.490471, 38.099937, 39.200233, 40.635674, 42.10854, 43.159355, 43.4], dtype=float),
                "I": np.array([4.8, 4.771824, 4.754732, 4.736563, 4.721648, 4.70373, 4.679128, 4.624589, 4.569901, 4.492322, 4.415938, 4.290498, 4.119548, 3.97847, 3.803591, 3.594777, 3.382138, 2.680002, 1.764004, 0.824125, 0.153565, 0], dtype=float),
            },
            {
                "name": 'Open-source PV Emulator - Simplified Exp. Model', "shortName": 'Simplified exp.',
                "color": (1, 0, 0),
                "lineStyle": '-',
                "V": np.array([0, 25.567003, 27.200571, 28.811605, 30.372927, 32.341831, 34.1399, 36.391037, 38.724979, 40.008648, 41.242683, 42.378609, 42.994045, 43.4], dtype=float),
                "I": np.array([4.8, 4.747, 4.733705, 4.687959, 4.643624, 4.51567, 4.370661, 3.898574, 3.404761, 2.759869, 2.083102, 0.986255, 0.391991, 0], dtype=float),
            },
            {
                "name": 'Commercial PV emulator', "shortName": 'Commercial',
                "color": (1, 0, 1),
                "lineStyle": '-',
                "V": np.array([1.55, 3.4, 8.22, 11.75, 14.85, 17.9, 19.52, 21.8, 24.55, 27.47, 30.13, 33.09, 35.9, 36.8505, 37.8, 38.5, 38.8934, 39.3103, 39.5975, 39.99, 40.2553, 40.7371, 41.228, 41.7376, 42.0248, 42.4047, 42.6456], dtype=float),
                "I": np.array([4.7949, 4.7884, 4.7695, 4.7529, 4.735, 4.715, 4.7027, 4.68, 4.65, 4.6095, 4.555, 4.4518, 4.1523, 3.7759, 3.23, 2.8289, 2.6045, 2.3639, 2.1974, 1.9679, 1.8173, 1.5363, 1.252, 0.9905, 0.7912, 0.571, 0.4338], dtype=float),
            },
        ],
    },
]

def prepare_curve(V, I):
    V, I = np.asarray(V, dtype=float).ravel(), np.asarray(I, dtype=float).ravel()
    if len(V) != len(I) or len(V) < 2 or not np.all(np.isfinite(V)) or not np.all(np.isfinite(I)):
        raise ValueError("A curve requires equally sized finite arrays and at least two points.")
    order = np.argsort(V, kind="stable")
    V, I = V[order], I[order]
    V, indices = np.unique(V, return_index=True)
    if len(V) < 2:
        raise ValueError("A curve requires at least two distinct voltages.")
    return V, I[indices]


def estimate_isc_from_curve(V, I):
    V, I = prepare_curve(V, I)
    zeros = np.flatnonzero(np.abs(V) < 1e-9)
    if len(zeros):
        return float(I[zeros[0]])
    return float(interp1d(V, I, kind="linear", fill_value="extrapolate")(0.0))


def estimate_voc_from_curve(V, I):
    V, I = prepare_curve(V, I)
    zeros = np.flatnonzero(np.abs(I) < 1e-9)
    if len(zeros):
        return float(V[zeros[-1]])
    low = I <= 0.20 * np.max(I)
    if np.count_nonzero(low) < 2:
        low = np.zeros(len(I), dtype=bool)
        low[-4:] = True
    Ilow, Vlow = I[low], V[low]
    order = np.argsort(Ilow, kind="stable")
    Ilow, Vlow = Ilow[order], Vlow[order]
    Ilow, indices = np.unique(Ilow, return_index=True)
    Vlow = Vlow[indices]
    if len(Ilow) < 2:
        raise ValueError("Voc extrapolation requires two distinct low-current samples.")
    return float(interp1d(Ilow, Vlow, kind="linear", fill_value="extrapolate")(0.0))


def find_mpp_from_curve(V, I, n_points=3000):
    V, I = prepare_curve(V, I)
    Vgrid = np.linspace(V.min(), V.max(), n_points)
    Igrid = np.maximum(PchipInterpolator(V, I)(Vgrid), 0)
    power = Vgrid * Igrid
    index = np.argmax(power)
    return float(Vgrid[index]), float(Igrid[index]), float(power[index])


def characteristic_points(V, I):
    Vmpp, Impp, Pmpp = find_mpp_from_curve(V, I)
    return np.array([estimate_isc_from_curve(V, I), estimate_voc_from_curve(V, I),
                     Vmpp, Impp, Pmpp])




def augment_curve_with_isc_voc(V, I):
    """Add estimated (0, Isc) and (Voc, 0); keep first samples at duplicates."""
    V, I = prepare_curve(V, I)
    Isc, Voc = estimate_isc_from_curve(V, I), estimate_voc_from_curve(V, I)
    return prepare_curve(np.r_[0, V, Voc], np.r_[Isc, I, 0])


def relative_error_percent(value, reference):
    return float(100 * np.abs(value - reference) / np.abs(reference))


def calculate_panel_metrics(panel, n_common_grid=N_COMMON_GRID):
    """Return all original MATLAB result columns for this panel's emulators.

    Characteristic points use the original curves. Common-grid metrics use
    endpoint-augmented curves, one shared overlap interval, PCHIP and clipped
    nonnegative current. The 10% MAPE masks use strict '>' as in MATLAB.
    """
    if n_common_grid < 2:
        raise ValueError("The common voltage grid requires at least two points.")
    Vref, Iref = prepare_curve(panel["V_panel"], panel["I_panel"])
    reference = characteristic_points(Vref, Iref)
    Isc_ref, Voc_ref, Vmpp_ref, Impp_ref, Pmpp_ref = reference
    reference_aug = augment_curve_with_isc_voc(Vref, Iref)
    emulators = [prepare_curve(emu["V"], emu["I"]) for emu in panel["emu"]]
    augmented = [augment_curve_with_isc_voc(V, I) for V, I in emulators]
    all_augmented = [reference_aug] + augmented
    grid_min = max(V.min() for V, _ in all_augmented)
    grid_max = min(V.max() for V, _ in all_augmented)
    if grid_max <= grid_min:
        raise ValueError(f"Invalid common voltage grid for panel {panel['name']}.")
    grid = np.linspace(grid_min, grid_max, n_common_grid)
    Iref_grid = np.maximum(PchipInterpolator(*reference_aug)(grid), 0)
    Pref_grid = grid * Iref_grid
    valid_I = Iref_grid > .10 * abs(Isc_ref)
    valid_P = Pref_grid > .10 * abs(Pmpp_ref)
    rows = []
    for emu, (Vemu, Iemu), augmented_curve in zip(panel["emu"], emulators, augmented):
        Iemu_grid = np.maximum(PchipInterpolator(*augmented_curve)(grid), 0)
        Pemu_grid = grid * Iemu_grid
        RMSE_I = np.sqrt(np.mean((Iemu_grid - Iref_grid)**2))
        RMSE_P = np.sqrt(np.mean((Pemu_grid - Pref_grid)**2))
        MAPE_I = (100 * np.mean(np.abs((Iemu_grid[valid_I] - Iref_grid[valid_I])
                                       / Iref_grid[valid_I])) if np.any(valid_I) else np.nan)
        MAPE_P = (100 * np.mean(np.abs((Pemu_grid[valid_P] - Pref_grid[valid_P])
                                       / Pref_grid[valid_P])) if np.any(valid_P) else np.nan)
        values = characteristic_points(Vemu, Iemu)
        row = {
            "PanelName": panel["name"], "Condition": panel["condition"],
            "Emulator": emu["shortName"], "Npoints": len(grid),
            "CommonGrid_Vmin_V": float(grid_min), "CommonGrid_Vmax_V": float(grid_max),
        }
        for label, unit, ref, value in zip(("Isc", "Voc", "Vmpp", "Impp", "Pmpp"),
                                          ("A", "V", "V", "A", "W"), reference, values):
            row[f"{label}_realPV_{unit}"] = float(ref)
            row[f"{label}_emu_{unit}"] = float(value)
            row[f"{label}_error_percent"] = relative_error_percent(value, ref)
        row.update({
            "RMSE_I_A": float(RMSE_I), "nRMSE_I_percent_Isc": float(100 * RMSE_I / abs(Isc_ref)),
            "MAPE_I_percent": float(MAPE_I), "RMSE_P_W": float(RMSE_P),
            "nRMSE_P_percent_Pmpp": float(100 * RMSE_P / abs(Pmpp_ref)),
            "MAPE_P_percent": float(MAPE_P), "Npoints_MAPE_I": int(valid_I.sum()),
            "Npoints_MAPE_P": int(valid_P.sum()),
        })
        rows.append(row)
    return rows


def draw_fixed_speech_callout(ax, x_tip, y_tip, x_center, y_top, lines,
                              tail_corner, font_size=FONT_CALLOUT):
    """Original fixed box placement, text measurement and triangular tail."""
    xlims, ylims = ax.get_xlim(), ax.get_ylim()
    xrange, yrange = np.diff(xlims)[0], np.diff(ylims)[0]
    text = "\n".join(lines)
    temp = ax.text(x_center, y_top, text, fontsize=font_size,
                   ha="center", va="center", alpha=0)
    ax.figure.canvas.draw()
    extent = temp.get_window_extent(ax.figure.canvas.get_renderer()).transformed(ax.transData.inverted())
    temp.remove()
    box_w, box_h = 1.18 * extent.width, 1.30 * extent.height
    x_center = np.clip(x_center, xlims[0] + .025*xrange + box_w/2,
                       xlims[1] - .025*xrange - box_w/2)
    y_top = np.clip(y_top, ylims[0] + .06*yrange + box_h, ylims[1] - .035*yrange)
    left, right, top, bottom = x_center-box_w/2, x_center+box_w/2, y_top, y_top-box_h
    tri_w, tri_h = .10 * box_w, .18 * box_h
    corners = {
        "nw": [(left, top-tri_h), (left+tri_w, top)],
        "ne": [(right-tri_w, top), (right, top-tri_h)],
        "sw": [(left, bottom+tri_h), (left+tri_w, bottom)],
        "se": [(right-tri_w, bottom), (right, bottom+tri_h)],
    }
    if tail_corner.lower() not in corners:
        raise ValueError("tail_corner must be nw, ne, sw or se.")
    style = dict(facecolor=(.94, .94, .94), edgecolor=(.45, .45, .45), linewidth=.9)
    ax.add_patch(Polygon(corners[tail_corner.lower()] + [(x_tip, y_tip)], **style, zorder=4))
    ax.add_patch(Polygon([(left,bottom), (right,bottom), (right,top), (left,top)], **style, zorder=5))
    ax.text(x_center, .5*(bottom+top), text, fontsize=font_size,
            ha="center", va="center", color=(.1,.1,.1), zorder=6)


def plot_panel_iv(panel, metrics):
    Vref, Iref = prepare_curve(panel["V_panel"], panel["I_panel"])
    ref = characteristic_points(Vref, Iref)
    Isc_ref, Voc_ref, Vmpp_ref, Impp_ref, _ = ref
    fig, ax = plt.subplots(figsize=(11, 6.5))
    # Fix axes geometry before measuring callout text.
    fig.subplots_adjust(left=.115, right=.98, bottom=.16, top=.97)
    ax.tick_params(labelsize=FONT_TICKS)
    handles, labels = [], []
    for emu in panel["emu"]:
        handle, = ax.plot(emu["V"], emu["I"], emu["lineStyle"],
                          color=emu["color"], linewidth=LINE_WIDTH_EMU)
        handles.append(handle)
        labels.append(emu["name"].replace(" - ", " -\n"))
    handle, = ax.plot(Vref, Iref, "ko", markersize=MARKER_SIZE_REF, markerfacecolor="k")
    handles.append(handle)
    labels.append("Real PV panel data")
    for (x, y), marker, fill in zip(((0,Isc_ref), (Vmpp_ref,Impp_ref), (Voc_ref,0)),
                                    ("s","d","o"), ("y","g","c")):
        ax.plot(x, y, "k"+marker, markersize=9, markerfacecolor=fill)
    all_V = [*Vref, Voc_ref, Vmpp_ref]
    all_I = [*Iref, Isc_ref, Impp_ref]
    for emu, row in zip(panel["emu"], metrics):
        for x, y, marker in ((0, row["Isc_emu_A"], "s"),
                              (row["Vmpp_emu_V"], row["Impp_emu_A"], "d"),
                              (row["Voc_emu_V"], 0, "o")):
            ax.plot(x, y, marker, markersize=8, markeredgecolor=emu["color"],
                    markerfacecolor=emu["color"])
        all_V.extend([*emu["V"], row["Voc_emu_V"], row["Vmpp_emu_V"]])
        all_I.extend([*emu["I"], row["Isc_emu_A"], row["Impp_emu_A"]])
    xmax, ymax = max(all_V), max(all_I)
    ax.set(xlim=(0,1.05*xmax), ylim=(0,1.18*ymax))
    ax.set_xlabel(r"$V_{out}\,[\mathrm{V}]$", fontsize=FONT_AXIS_LABEL)
    ax.set_ylabel(r"$I_{out}\,[\mathrm{A}]$", fontsize=FONT_AXIS_LABEL)
    ax.grid(True)
    ax.legend(handles, labels, loc="lower left", fontsize=FONT_LEGEND)
    txt_isc, txt_voc = [r"$E_{I_{sc}}$"], [r"$E_{V_{oc}}$"]
    txt_mpp = [r"$E_{V_{mpp}}/E_{I_{mpp}}/E_{P_{mpp}}$ [%]"]
    for row in metrics:
        short = row["Emulator"]
        txt_isc.append(f"{short}: {row['Isc_error_percent']:.2f}%")
        txt_voc.append(f"{short}: {row['Voc_error_percent']:.2f}%")
        txt_mpp.append(f"{short}: {row['Vmpp_error_percent']:.2f} / "
                       f"{row['Impp_error_percent']:.2f} / {row['Pmpp_error_percent']:.2f}%")
    draw_fixed_speech_callout(ax,0,Isc_ref,.17*xmax,.78*ymax,txt_isc,"nw")
    draw_fixed_speech_callout(ax,Vmpp_ref,Impp_ref,.58*xmax,.67*ymax,txt_mpp,"ne")
    draw_fixed_speech_callout(ax,Voc_ref,0,.80*xmax,.31*ymax,txt_voc,"se")
    return fig


def print_summary(results):
    print("\n================ VALIDATION SUMMARY ================")
    for panel_name in dict.fromkeys(row["PanelName"] for row in results):
        rows = [row for row in results if row["PanelName"] == panel_name]
        print(f"\n{panel_name}")
        print("  Curve errors [%]")
        print(f"  {'Emulator':<16} {'nRMSE_I':>9} {'MAPE_I':>9} {'nRMSE_P':>9} {'MAPE_P':>9}")
        for row in rows:
            print(f"  {row['Emulator']:<16} {row['nRMSE_I_percent_Isc']:9.2f} "
                  f"{row['MAPE_I_percent']:9.2f} {row['nRMSE_P_percent_Pmpp']:9.2f} {row['MAPE_P_percent']:9.2f}")
        print("\n  Characteristic-point errors [%]")
        print(f"  {'Emulator':<16} {'Isc':>7} {'Voc':>7} {'Vmpp':>7} {'Impp':>7} {'Pmpp':>7}")
        for row in rows:
            values = ' '.join(f"{row[f'{label}_error_percent']:7.2f}" for label in ("Isc","Voc","Vmpp","Impp","Pmpp"))
            print(f"  {row['Emulator']:<16} {values}")
    print("\n====================================================")


def plot_summary(results, quantity):
    suffix, normalizer, ylabel = (("I", "Isc", "Current error metrics [%]") if quantity == "current"
                                  else ("P", "Pmpp", "Power error metrics [%]"))
    positions, width = np.arange(len(results)), .36
    # Taller than the MATLAB window so the original 30 pt y-label and long
    # rotated panel/emulator names remain fully visible.
    fig, ax = plt.subplots(figsize=(16, 9), layout="constrained")
    ax.bar(positions-width/2, [row[f"nRMSE_{suffix}_percent_{normalizer}"] for row in results],
           width, label=rf"$nRMSE_{suffix}$")
    ax.bar(positions+width/2, [row[f"MAPE_{suffix}_percent"] for row in results],
           width, label=rf"$MAPE_{suffix}$")
    ax.set_xticks(positions, [f"{row['PanelName']} - {row['Emulator']}" for row in results],
                  rotation=35, ha="right", fontsize=16)
    ax.tick_params(axis="y", labelsize=16)
    ax.set_ylabel(ylabel, fontsize=FONT_AXIS_LABEL)
    ax.legend(loc="best", fontsize=FONT_LEGEND)
    ax.grid(False)
    return fig


def main(argv=None):
    args = plot_arguments(__doc__).parse_args(argv)
    configure_plots(FONT_TICKS)
    figures, results = [], []
    for index, panel in enumerate(PANELS, 1):
        metrics = calculate_panel_metrics(panel)
        results.extend(metrics)
        figures.append((f"comparison_panel_{index}", plot_panel_iv(panel, metrics)))
    print_summary(results)
    for quantity in ("current", "power"):
        figures.append((f"comparison_{quantity}_metrics", plot_summary(results, quantity)))
    print(f"\nValidation completed for {len(PANELS)} panels and "
          f"{len(PANELS[0]['emu'])} emulators using common voltage grids.")
    finish(figures, args)
    return results


if __name__ == "__main__":
    main()
