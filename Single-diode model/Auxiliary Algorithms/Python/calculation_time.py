"""Plot the embedded measured calculation time versus output current/voltage."""
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

I_raw = np.array([3.02, 2.99, 2.96, 2.88, 2.86, 2.68, 2.58, 2.21, 2.02, 1.76, 1.50, 1.29, 1.11, 0.98, 0.73, 0.52, 0.38, 0.27], dtype=float)
V_raw = np.array([14.78, 16.27, 16.65, 17.20, 17.55, 18.38, 18.75, 19.62, 20.03, 20.33, 20.50, 20.66, 20.86, 21.03, 21.31, 21.49, 21.74, 21.69], dtype=float)
I_capt = np.array([3.000, 2.955, 2.933, 2.888, 2.840, 2.715, 2.613, 2.290, 2.130, 1.910, 1.727, 1.620, 1.456, 1.340, 1.140, 0.995, 0.830, 0.830], dtype=float)
V_capt = np.array([14.8100, 16.2700, 16.6700, 17.1860, 17.5700, 18.3000, 18.7016, 19.5300, 19.8750, 20.1500, 20.3500, 20.4800, 20.6600, 20.7850, 21.0130, 21.1800, 21.3600, 21.3600], dtype=float)
t_ms = np.array([3.11418, 3.74697, 4.15935, 4.79214, 5.21163, 6.05061, 6.67629, 6.88248, 6.88248, 7.08867, 7.09578, 7.08156, 7.09578, 7.09578, 7.09578, 7.08156, 7.08867, 7.08867], dtype=float)

# Same default selection as MATLAB. --raw selects the alternative raw data.
USE_CAPTURED_DATA = True


def main(argv=None):
    parser = plot_arguments(__doc__)
    parser.add_argument("--raw", action="store_true", help="Use raw instead of captured operating points.")
    args = parser.parse_args(argv)
    configure_plots()
    use_captured = USE_CAPTURED_DATA and not args.raw
    I_data, V_data = (I_capt, V_capt) if use_captured else (I_raw, V_raw)
    figures = []
    for name, values, xlabel in (
        ("calculation_time_current", I_data, r"$I_{out}^{OT}\,[\mathrm{A}]$"),
        ("calculation_time_voltage", V_data, r"$V_{out}^{OT}\,[\mathrm{V}]$"),
    ):
        fig, ax = plt.subplots(figsize=(9, 5), layout="constrained")
        ax.plot(values, t_ms, "-ko", linewidth=2, markersize=6, markerfacecolor="k")
        ax.grid(True)
        ax.set(xlabel=xlabel, ylabel=r"$t_{calc}\,[\mathrm{ms}]$")
        figures.append((name, fig))
    finish(figures, args)


if __name__ == "__main__":
    main()
