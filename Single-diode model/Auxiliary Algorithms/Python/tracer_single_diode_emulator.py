"""Plot the single-diode model, resistive load lines and emulator measurements."""
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

from scipy.optimize import least_squares

Q = 1.60217662e-19
K = 1.38064852e-23
E_G0 = 1.166
K1 = 4.73e-4
K2 = 636.0


def residuals_2_20(x, Voc, Isc, Vmp, Imp, Tref):
    """The five original MATLAB equations, in physical parameter coordinates.

    x = [Iph (A), Is0 (A), A, Rs (ohm), Rp (ohm)].
    A is the module effective ideality factor, usually n * ns.
    """
    Iph, Is0, A, Rs, Rp = x
    C = Q / (A * K * Tref)
    exp_sc = np.exp(np.clip(C * Isc * Rs, -700, 700))
    exp_oc = np.exp(np.clip(C * Voc, -700, 700))
    exp_mp = np.exp(np.clip(C * (Vmp + Imp * Rs), -700, 700))
    return np.array([
        Iph - Is0 * (exp_sc - 1) - Isc * Rs / Rp - Isc,
        Iph - Is0 * (exp_oc - 1) - Voc / Rp,
        Iph - Is0 * (exp_mp - 1) - (Vmp + Imp * Rs) / Rp - Imp,
        Iph - 2 * Vmp / Rp - Is0 * ((1 + C * (Vmp - Imp * Rs)) * exp_mp - 1),
        Rs + Is0 * C * Rp * (Rs - Rp) * exp_sc,
    ])


def estimate_parameters(Voc, Isc, Vmp, Imp, ns, Tref=298.15):
    """Solve the original five equations with consistent positive coordinates.

    The MATLAB initial vector mixes log10(Is0), log10(Rs) with physical
    values, although residuals_2_20 uses only physical values. Here ALL
    parameters are optimized in log coordinates and decoded before evaluation.
    Scaling changes conditioning only, not the roots of the equations.
    A solution is accepted only after checking the original residuals.
    """
    if not (Voc > Vmp > 0 and Isc > Imp > 0 and ns > 0 and Tref > 0):
        raise ValueError("Require Voc > Vmp > 0, Isc > Imp > 0, ns > 0 and Tref > 0.")
    Rs0 = (Voc - Vmp) / Imp
    Rp0 = Vmp / (Isc - Imp)
    scale = np.array([Isc, Isc, Isc, Isc, max(Rs0, 0.1)])
    lower = np.log([Isc * 0.5, 1e-30, ns * 0.05, 1e-12, 1e-3])
    upper = np.log([Isc * 2, 1.0, ns * 5, Voc / Isc, 1e9])

    def objective(log_x):
        return residuals_2_20(np.exp(log_x), Voc, Isc, Vmp, Imp, Tref) / scale

    best = None
    for A0 in (ns, 0.5 * ns, 1.5 * ns):
        x0 = np.array([Isc, 1e-9, A0, min(Rs0 * 0.2, Voc / Isc * 0.9), Rp0])
        result = least_squares(objective, np.log(x0), bounds=(lower, upper),
                               ftol=1e-13, xtol=1e-13, gtol=1e-13, max_nfev=2000)
        error = np.max(np.abs(objective(result.x)))
        if best is None or error < best[0]:
            best = (error, np.exp(result.x))
        if result.success and error < 1e-9:
            break
    if best[0] > 1e-8:
        raise RuntimeError(f"Parameter estimation did not converge: scaled residual = {best[0]:.3g}")
    return best[1]


def print_parameters(parameters):
    Iph_ref, Is0_ref, A, Rs, Rp = parameters
    print(f"Iph_ref = {Iph_ref:.6f} A")
    print(f"Is0_ref = {Is0_ref:.2e} A")
    print(f"A       = {A:.6f}")
    print(f"Rs      = {Rs:.6f} Ohms")
    print(f"Rp      = {Rp:.6f} Ohms")


def solve_I_V_2_11(V, Iph_ref, Is0_ref, A, Rs, Rp, G=1000.0,
                   Gref=1000.0, alpha=0.0003, T=298.15, Tref=298.15, ns=72):
    """Original temperature/irradiance law and Newton-Raphson current solve.

    V can be a scalar or array. Temperature arguments are in kelvin;
    alpha is a relative coefficient in 1/K. Negative current is retained,
    as in MATLAB, if the applied voltage exceeds the model's actual Voc.
    """
    if G < 0 or Gref <= 0 or T <= 0 or Tref <= 0:
        raise ValueError("Invalid irradiance or absolute temperature.")
    V = np.asarray(V, dtype=float)
    Iph = Iph_ref * (G / Gref) * (1 + alpha * (T - Tref))
    Eg_T_eV = E_G0 - K1 * T**2 / (T + K2)
    exponent = ns * Q / (A * K) * Eg_T_eV * (1 / Tref - 1 / T)
    Is = Is0_ref * (T / Tref)**3 * np.exp(exponent)
    current = np.full_like(V, Iph)
    active = np.ones(V.shape, dtype=bool)
    for _ in range(30):
        diode_voltage = V + current * Rs
        expo = np.exp(np.minimum(Q * diode_voltage / (A * K * T), 700))
        f = Iph - Is * (expo - 1) - diode_voltage / Rp - current
        df = -Is * expo * Q * Rs / (A * K * T) - Rs / Rp - 1
        delta = -f / df
        current = np.where(active, current + delta, current)
        active &= np.abs(delta) >= 1e-6
        if not np.any(active):
            return float(current) if current.ndim == 0 else current
    raise RuntimeError("Newton-Raphson did not converge within 30 iterations.")


# Editable reference module and operating conditions (Shell Solar SQ150-PC).
NS = 72
VMP_MOD_REF = 34.0
IMP_MOD_REF = 4.4
VOC_MOD_REF = 43.4
ISC_MOD_REF = 4.8
TREF = 25 + 273.15
GREF = 1000.0
ALPHA = 0.0003  # relative Isc temperature coefficient [1/K]
BETA = -0.0037  # relative Voc temperature coefficient [1/K]
TEMPERATURE_C = 25.0
IRRADIANCE = 1000.0

FACTORS = np.array([0.00, 0.20, 0.40, 0.60, 0.70, 0.80, 0.85, 0.88, 0.90, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03, 1.04, 1.05, 1.08, 1.12], dtype=float)

# Active MATLAB measurement rows. Columns:
# V_test [V], I_test [A], V_emulator [V], I_emulator [A], R_load [ohm].
MEASUREMENTS = np.array([
    [45.001869, 9.165687, 23.428795, 4.771824, 4.909819599],
    [45.002346, 8.243192, 25.957676, 4.754732, 5.459335247],
    [45.00293, 7.68873, 27.723593, 4.736563, 5.8531034],
    [45.002434, 7.429503, 28.600252, 4.721648, 6.057260516],
    [45.006149, 7.207458, 29.371901, 4.70373, 6.244384988],
    [45.009201, 6.976861, 30.186039, 4.679128, 6.451210354],
    [45.00576, 6.623515, 31.42337, 4.624589, 6.794845985],
    [45.007038, 6.370541, 32.285755, 4.569901, 7.064869677],
    [45.006256, 6.092615, 33.184868, 4.492322, 7.387019007],
    [45.003723, 5.867533, 33.870052, 4.415938, 7.669956417],
    [45.006004, 5.55717, 34.74757, 4.290498, 8.098726535],
    [45.006401, 5.199488, 35.65852, 4.119548, 8.655930214],
    [45.007046, 4.942261, 36.230217, 3.97847, 9.106570365],
    [45.009201, 4.642657, 36.874702, 3.803591, 9.694707449],
    [45.006428, 4.315445, 37.490471, 3.594777, 10.42915068],
    [45.005337, 3.995132, 38.099937, 3.382138, 11.26504507],
    [45.006454, 3.076956, 39.200233, 2.680002, 14.6269417],
    [45.010605, 1.953921, 40.635674, 1.764004, 23.03604414],
    [45.032486, 0.160229, 43.159355, 0.153565, 281.0494253],
], dtype=float)

# Inactive measurement datasets preserved from the original MATLAB comments.
# To use one, also update the reference module and G/T settings above.
ALTERNATIVE_MEASUREMENTS = {
    'CanadianSolar CS6P-250P': np.array([
        [36.470619, 10.97875, 22.642405, 6.816042, 3.321928621],
        [36.469147, 10.580585, 23.414888, 6.793228, 3.446798488],
        [36.011593, 10.135595, 24.056635, 6.770828, 3.552982737],
        [36.008194, 9.758461, 24.833183, 6.729959, 3.689945659],
        [36.011219, 9.268014, 25.833637, 6.648664, 3.885538057],
        [36.012985, 8.819846, 26.70089, 6.539245, 4.083176269],
        [36.00943, 8.417071, 27.416254, 6.408448, 4.278142539],
        [36.009499, 8.022567, 28.057871, 6.251021, 4.488526114],
        [36.011665, 7.394004, 28.952494, 5.944597, 4.870388018],
        [36.014778, 6.823575, 29.632019, 5.614259, 5.277992875],
        [36.012032, 6.198761, 30.29425, 5.214557, 5.809553908],
        [36.01263, 5.744209, 30.687252, 4.894783, 6.26937946],
        [36.012844, 5.41966, 30.963587, 4.659785, 6.644853142],
        [36.013351, 4.885607, 31.42931, 4.263731, 7.371316342],
        [36.014511, 4.044705, 31.895897, 3.582153, 8.904113532],
        [36.019165, 2.940577, 32.502819, 2.653505, 12.24901366],
        [36.026188, 1.518282, 33.319305, 1.404204, 23.72825102],
        [36.04068, 0.271503, 34.06929, 0.256652, 132.7450789],
    ], dtype=float),
    'Kyocera KB260-6BPA': np.array([
        [38.970139, 10.476988, 24.569763, 6.605496, 3.719593956],
        [38.970062, 11.051579, 23.384464, 6.631636, 3.526198362],
        [38.9701, 10.006746, 25.593412, 6.571879, 3.894382718],
        [38.478722, 9.397497, 26.668896, 6.513233, 4.09457116],
        [38.481834, 8.949709, 27.661718, 6.433278, 4.299785895],
        [38.481712, 8.526176, 28.55562, 6.326908, 4.513361029],
        [38.97081, 8.168219, 29.45533, 6.17379, 4.771028817],
        [38.484024, 7.527627, 30.391579, 5.944713, 5.11237111],
        [38.485474, 6.842545, 31.405592, 5.583774, 5.624438238],
        [38.480385, 6.376169, 31.972301, 5.297785, 6.03503181],
        [38.480625, 6.013403, 32.411503, 5.064976, 6.399142464],
        [38.485367, 5.482672, 32.913673, 4.688922, 7.019454152],
        [38.481773, 5.186515, 33.195103, 4.473986, 7.419581331],
        [38.482933, 4.712598, 33.651192, 4.120906, 8.165969328],
        [38.483788, 3.500327, 34.301086, 3.119885, 10.99434306],
        [38.970074, 1.690633, 35.331486, 1.532781, 23.0505767],
        [38.513306, 0.255853, 36.170616, 0.24029, 150.5290108],
    ], dtype=float),
    'Kyocera KC85TS - R load': np.array([
        [22.491886, 11.192251, 4.290837, 2.135176, 2.009594057],
        [22.785145, 4.797872, 10.134062, 2.133931, 4.749011097],
        [22.49477, 3.785123, 12.669167, 2.1318, 5.94294352],
        [22.494724, 3.463998, 13.826797, 2.12921, 6.493862512],
        [22.785276, 3.220089, 15.001449, 2.120053, 7.075978289],
        [22.492458, 3.039437, 15.612347, 2.109718, 7.40020562],
        [22.494221, 2.910696, 16.172581, 2.092691, 7.728126608],
        [22.784504, 2.841723, 16.606222, 2.071157, 8.017847995],
        [22.493919, 2.67594, 17.095009, 2.03367, 8.405989664],
        [22.494114, 2.549509, 17.511776, 1.984805, 8.822920136],
        [22.494652, 2.356512, 18.038937, 1.889737, 9.545739434],
        [22.495331, 2.168992, 18.432518, 1.777257, 10.37132953],
        [22.497145, 1.950675, 18.831503, 1.632836, 11.53300331],
        [22.500511, 1.719844, 19.143658, 1.46326, 13.08288206],
        [22.504717, 1.484131, 19.473244, 1.284213, 15.16356243],
        [22.512615, 0.991903, 19.819204, 0.873231, 22.6964045],
        [22.518124, 0.453164, 20.202765, 0.406569, 49.69086428],
        [22.51823, 0.317759, 20.30147, 0.286478, 70.86572093],
    ], dtype=float),
    'Kyocera KC85TS - R+L load': np.array([
        [22.491726, 11.401453, 4.212106, 2.135191, 1.972706891],
        [22.784645, 4.763543, 10.206792, 2.133915, 4.7831296],
        [22.494654, 3.790882, 12.650108, 2.131843, 5.933883499],
        [22.494001, 3.46787, 13.811143, 2.129245, 6.486403866],
        [22.494892, 3.18538, 14.974236, 2.12042, 7.061919808],
        [22.493298, 3.035744, 15.629344, 2.10937, 7.409484348],
        [22.49297, 2.908536, 16.181227, 2.092373, 7.733433284],
        [22.495255, 2.81023, 16.588015, 2.072266, 8.004771106],
        [22.493977, 2.670704, 17.113245, 2.031851, 8.422490133],
        [22.493221, 2.52665, 17.581903, 1.974965, 8.902387131],
        [22.494623, 2.340781, 18.077917, 1.881181, 9.609876455],
        [22.493725, 2.189714, 18.393583, 1.790574, 10.27245062],
        [22.785147, 1.989317, 18.813585, 1.64257, 11.45374931],
        [22.501472, 1.710468, 19.156639, 1.456207, 13.15516201],
        [22.504978, 1.49112, 19.463346, 1.28959, 15.09266201],
        [22.512096, 1.000689, 19.813057, 0.880713, 22.49661013],
        [22.518711, 0.457511, 20.199619, 0.410394, 49.22006413],
        [22.517979, 0.306336, 20.309837, 0.276297, 73.50726573],
    ], dtype=float),
    'Kyocera KC200GT': np.array([
        [30.799538, 10.383743, 12.533917, 4.225679, 2.966130887],
        [30.407265, 8.72176, 14.705232, 4.217923, 3.486368054],
        [30.411446, 7.384953, 17.256174, 4.190397, 4.118028435],
        [30.412125, 7.033544, 18.053396, 4.175287, 4.323869473],
        [30.408638, 6.741723, 18.770275, 4.161449, 4.510514246],
        [30.799961, 6.556941, 19.421268, 4.134554, 4.69730665],
        [30.411461, 6.154503, 20.231714, 4.094382, 4.941335225],
        [30.410652, 5.776842, 21.174362, 4.022306, 5.264234496],
        [30.410862, 5.474993, 21.89076, 3.941084, 5.554502264],
        [30.411789, 5.199258, 22.50289, 3.847137, 5.849256213],
        [30.411108, 4.899352, 23.111118, 3.723294, 6.20716978],
        [30.409916, 4.590909, 23.679567, 3.574845, 6.623942297],
        [30.410358, 4.117438, 24.440063, 3.309085, 7.385746513],
        [30.410532, 3.747369, 24.939966, 3.073253, 8.11516852],
        [30.412104, 3.481879, 25.292833, 2.895774, 8.734394673],
        [30.413876, 3.004618, 25.808134, 2.549612, 10.12237705],
        [30.411911, 2.620044, 26.232281, 2.259961, 11.60740429],
        [30.799515, 1.95355, 26.693003, 1.693082, 15.76592451],
        [30.428387, 1.320891, 27.100122, 1.176412, 23.03625091],
        [30.434181, 0.844936, 27.427097, 0.761451, 36.01951669],
        [30.799423, 0.298323, 27.814808, 0.269414, 103.2418805],
    ], dtype=float),
    'ME Solar MESM-50W': np.array([
        [23.415052, 11.967575, 5.922083, 3.026813, 1.956540758],
        [23.112965, 5.560591, 12.551384, 3.019652, 4.156566386],
        [23.116062, 5.030003, 13.841137, 3.0118, 4.595636164],
        [23.115747, 4.57807, 15.125275, 2.995558, 5.049234567],
        [23.115692, 4.439414, 15.540461, 2.984576, 5.2069242],
        [23.114132, 4.291178, 15.992638, 2.96906, 5.386431396],
        [23.116524, 4.171403, 16.356771, 2.951598, 5.541666243],
        [23.114021, 4.05022, 16.715565, 2.929032, 5.706856395],
        [23.115362, 3.948944, 17.007088, 2.905429, 5.853554845],
        [23.114918, 3.798939, 17.418226, 2.862687, 6.084572292],
        [23.114498, 3.65326, 17.789234, 2.811599, 6.327087896],
        [23.114141, 3.503831, 18.138689, 2.749611, 6.596820059],
        [23.114452, 3.3509, 18.464882, 2.676852, 6.8979839],
        [23.11492, 3.15309, 18.834526, 2.569205, 7.330877061],
        [23.116348, 2.903193, 19.264559, 2.419445, 7.962387655],
        [23.115911, 2.630473, 19.65144, 2.236233, 8.787742601],
        [23.116068, 2.396255, 19.971132, 2.070245, 9.64674809],
        [23.123291, 1.748199, 20.552086, 1.553807, 13.22692329],
        [23.136288, 0.80102, 21.464041, 0.743124, 28.8835255],
        [23.140648, 0.309762, 21.969181, 0.29408, 74.70477761],
    ], dtype=float),
    'Renogy RNG-50DB-H': np.array([
        [23.424229, 11.58681, 5.894531, 2.915734, 2.021628516],
        [23.72953, 5.693411, 12.125871, 2.909352, 4.167894088],
        [23.426764, 4.922009, 13.816704, 2.902917, 4.759593195],
        [23.424603, 4.170353, 16.145973, 2.874516, 5.616936208],
        [23.426395, 4.002863, 16.715784, 2.856222, 5.852410632],
        [23.423584, 3.848856, 17.231913, 2.831469, 6.08585614],
        [23.424973, 3.715216, 17.664825, 2.801653, 6.305143785],
        [23.426899, 3.545583, 18.180428, 2.751547, 6.607347794],
        [23.427967, 3.344581, 18.727219, 2.673501, 7.004754814],
        [23.427265, 3.14929, 19.188766, 2.579515, 7.4389046],
        [23.427248, 2.939382, 19.602093, 2.459445, 7.970128627],
        [23.426811, 2.739423, 19.970377, 2.335244, 8.551730355],
        [23.426422, 2.418296, 20.420042, 2.107949, 9.687161312],
        [23.429234, 2.159837, 20.759869, 1.91376, 10.84768675],
        [23.437695, 1.558239, 21.242069, 1.412265, 15.04113534],
        [23.447426, 0.686716, 21.980999, 0.643768, 34.14428645],
        [23.450655, 0.307881, 22.318256, 0.293014, 76.16788276],
    ], dtype=float),
}


def piecewise_coefficients(points_V, points_I):
    points_V, points_I = np.asarray(points_V), np.asarray(points_I)
    if len(points_V) < 2 or len(points_V) != len(points_I) or np.any(np.diff(points_V) <= 0):
        raise ValueError("Piecewise interpolation requires increasing voltage points.")
    a = np.diff(points_I) / np.diff(points_V)
    b = points_I[:-1] - a * points_V[:-1]
    return np.column_stack((a, b))


def piecewise_pv_model(V_in, points_V, coeffs):
    """I = a*V + b; return zero outside the defined segments, as in MATLAB."""
    values = np.asarray(V_in, dtype=float)
    result = np.zeros_like(values)
    found = np.zeros(values.shape, dtype=bool)
    for index, (a, b) in enumerate(coeffs):
        mask = (~found) & (values >= points_V[index] - 1e-9) & (values <= points_V[index + 1] + 1e-9)
        result = np.where(mask, a * values + b, result)
        found |= mask
    return float(result) if result.ndim == 0 else result


def main(argv=None):
    parser = plot_arguments(__doc__)
    parser.add_argument("--temperature", type=float, default=TEMPERATURE_C, help="Temperature in degC.")
    parser.add_argument("--irradiance", type=float, default=IRRADIANCE, help="Irradiance in W/m^2.")
    args = parser.parse_args(argv)
    configure_plots()
    T, G = args.temperature + 273.15, args.irradiance
    if T <= 0 or G < 0:
        parser.error("Temperature must exceed absolute zero and irradiance must be nonnegative.")
    parameters = estimate_parameters(VOC_MOD_REF, ISC_MOD_REF, VMP_MOD_REF,
                                     IMP_MOD_REF, NS, TREF)
    print_parameters(parameters)
    Voc_estimation = VOC_MOD_REF * (1 + BETA * (T - TREF))
    Vmp_estimation = VMP_MOD_REF * Voc_estimation / VOC_MOD_REF
    if Voc_estimation <= 0:
        raise ValueError("The linear Voc estimate is nonpositive at this temperature.")
    V_mod = FACTORS * Vmp_estimation
    V_mod = np.unique(np.append(V_mod[(V_mod < Voc_estimation) & (V_mod >= 0)], Voc_estimation))
    I_mod = solve_I_V_2_11(V_mod, *parameters, G=G, Gref=GREF, alpha=ALPHA,
                          T=T, Tref=TREF, ns=NS)
    coeffs = piecewise_coefficients(V_mod, I_mod)
    measurements = np.asarray(MEASUREMENTS, dtype=float)
    if measurements.size:
        V_test, I_test, V_emu, I_emu, R_data = measurements.T
        resistance_values = np.unique(R_data)
        vmax, imax = V_test.max(), I_test.max()
    else:
        # The MATLAB empty-measurements branch uses max([]), which cannot
        # define valid plot limits. Fall back to model extents here.
        V_test = I_test = V_emu = I_emu = np.array([])
        resistance_values = np.sort(np.r_[.1, .5, np.arange(1, 16),
                                         np.arange(20, 51, 5), np.arange(60, 201, 10),
                                         300, 400, 500, 1000, np.inf])
        vmax, imax = Voc_estimation, max(I_mod.max(), 1e-9)
    fig, ax = plt.subplots(figsize=(12, 7), layout="constrained")
    V_line = np.linspace(0, vmax * 1.25, 100)
    load_handle = None
    for r_val in resistance_values:
        if np.isinf(r_val):
            x, y = V_line, np.zeros_like(V_line)
        elif r_val == 0:
            x, y = np.zeros_like(V_line), np.linspace(0, imax * 1.1, 100)
        else:
            x, y = V_line, V_line / r_val
        handle, = ax.plot(x, y, "--", color=(.5, .5, .5), linewidth=.8)
        if load_handle is None:
            load_handle = handle
    V_plot = np.linspace(0, V_mod.max(), 500)
    model_handle, = ax.plot(V_plot, piecewise_pv_model(V_plot, V_mod, coeffs),
                            "r-", linewidth=2)
    handles, labels = [model_handle, load_handle], ["Single-diode model I–V curve", "Load line"]
    if measurements.size:
        handles.append(ax.scatter(V_test, I_test, s=40, color="m"))
        labels.append("Test point")
        handles.append(ax.scatter(V_emu, I_emu, s=70, color="k", marker="s"))
        labels.append("Emulation Points")
    ax.set(xlabel=r"$V_{out}^{OT}\,[\mathrm{V}]$", ylabel=r"$I_{out}^{OT}\,[\mathrm{A}]$",
           xlim=(0, 1.1 * Voc_estimation), ylim=(0, 1.1 * imax))
    ax.legend(handles, labels, loc="upper left", fontsize=20, frameon=True)
    ax.grid(False)
    finish([("single_diode_emulator", fig)], args)
    return parameters, V_mod, I_mod, coeffs


if __name__ == "__main__":
    main()
