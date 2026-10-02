"""Estimate the five single-diode parameters (CanadianSolar reference data)."""
import argparse
import numpy as np

from scipy.optimize import least_squares

Q = 1.60217662e-19
K = 1.38064852e-23
E_G0 = 1.166
K1 = 4.73e-4
K2 = 636.0


def residuals_2_20(x, Voc, Isc, Vmp, Imp, Tref):
    """
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


# Editable module parameters, copied from compute_5_parameters.m.
NS = 60
VMP_MOD_REF = 30.1
IMP_MOD_REF = 8.30
VOC_MOD_REF = 37.2
ISC_MOD_REF = 8.87
TREF = 25 + 273.15
GREF = 1000.0


def main():
    argparse.ArgumentParser(description=__doc__).parse_args()
    parameters = estimate_parameters(VOC_MOD_REF, ISC_MOD_REF, VMP_MOD_REF,
                                     IMP_MOD_REF, NS, TREF)
    print_parameters(parameters)
    return parameters


if __name__ == "__main__":
    main()
