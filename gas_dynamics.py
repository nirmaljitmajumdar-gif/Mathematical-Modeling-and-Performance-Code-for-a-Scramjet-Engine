"""
gas_dynamics.py
----------------
Compressible-flow relations used to march the flow through the scramjet:

  * Isentropic flow relations (stagnation-to-static ratios, area-Mach relation)
    -> used for the inlet (compression) and the nozzle (expansion).
  * Rayleigh-flow relations (frictionless, constant-area duct with heat addition)
    -> used for the combustor, including the thermal-choking limit.

All inverse relations (Mach number from an area ratio or a stagnation-
temperature ratio) are solved with a bisection root-find on the physically
correct branch (M > 1 throughout a scramjet flow path).
"""

import math

GAMMA_DEFAULT = 1.4


# ---------------------------------------------------------------------------
# Isentropic relations
# ---------------------------------------------------------------------------

def T0_over_T(M, gamma=GAMMA_DEFAULT):
    """Stagnation-to-static temperature ratio."""
    return 1.0 + 0.5 * (gamma - 1.0) * M ** 2


def p0_over_p(M, gamma=GAMMA_DEFAULT):
    """Stagnation-to-static pressure ratio (isentropic)."""
    return T0_over_T(M, gamma) ** (gamma / (gamma - 1.0))


def A_over_Astar(M, gamma=GAMMA_DEFAULT):
    """Isentropic area ratio A/A* as a function of Mach number."""
    g = gamma
    term = (2.0 / (g + 1.0)) * T0_over_T(M, g)
    return (1.0 / M) * term ** ((g + 1.0) / (2.0 * (g - 1.0)))


def solve_M_from_area_ratio(AR, gamma=GAMMA_DEFAULT, supersonic=True, M_guess_bounds=None):
    """
    Invert A/A* = AR for Mach number, on the supersonic (M>1) branch by
    default (the physically relevant branch for a scramjet).
    """
    if AR < 1.0 - 1e-9:
        raise ValueError(f"A/A* = {AR:.4f} < 1 is not physically achievable (min is 1 at M=1).")

    lo, hi = (M_guess_bounds if M_guess_bounds else ((1.0 + 1e-6, 60.0) if supersonic else (1e-4, 1.0 - 1e-6)))

    def f(M):
        return A_over_Astar(M, gamma) - AR

    return _bisect(f, lo, hi)


# ---------------------------------------------------------------------------
# Rayleigh-flow relations (constant-area duct, frictionless, with heat addition)
# ---------------------------------------------------------------------------

def rayleigh_T0_over_T0star(M, gamma=GAMMA_DEFAULT):
    """Stagnation-temperature ratio T0/T0* for Rayleigh flow."""
    g = gamma
    num = (g + 1.0) * M ** 2 * (2.0 + (g - 1.0) * M ** 2)
    den = (1.0 + g * M ** 2) ** 2
    return num / den


def rayleigh_p_over_pstar(M, gamma=GAMMA_DEFAULT):
    """Static-pressure ratio p/p* for Rayleigh flow."""
    return (gamma + 1.0) / (1.0 + gamma * M ** 2)


def rayleigh_T_over_Tstar(M, gamma=GAMMA_DEFAULT):
    """Static-temperature ratio T/T* for Rayleigh flow."""
    return (M * (gamma + 1.0) / (1.0 + gamma * M ** 2)) ** 2


def solve_M_from_rayleigh_T0ratio(T0_ratio, gamma=GAMMA_DEFAULT, supersonic=True):
    """
    Invert T0/T0* = T0_ratio for Mach number on the supersonic branch
    (this is the branch that applies to a scramjet combustor, since the
    flow stays supersonic all the way through in true "scram" operation).

    T0/T0* increases monotonically from 0 to 1 as M -> 1 on the supersonic
    branch (M decreasing from infinity down to 1), so T0_ratio must be <= 1;
    T0_ratio == 1 is the thermal-choking limit (M -> 1).
    """
    if T0_ratio > 1.0 + 1e-9:
        raise ValueError("Requested stagnation-temperature ratio exceeds 1 -> thermal choking exceeded.")
    T0_ratio = min(T0_ratio, 1.0 - 1e-9)

    lo, hi = (1.0 + 1e-6, 60.0) if supersonic else (1e-4, 1.0 - 1e-6)

    def f(M):
        return rayleigh_T0_over_T0star(M, gamma) - T0_ratio

    return _bisect(f, lo, hi)


# ---------------------------------------------------------------------------
# Generic bisection root-finder (keeps this module dependency-free)
# ---------------------------------------------------------------------------

def _bisect(f, lo, hi, tol=1e-10, max_iter=200):
    f_lo, f_hi = f(lo), f(hi)
    if f_lo * f_hi > 0:
        raise ValueError(f"Bisection bracket does not straddle a root: f({lo})={f_lo}, f({hi})={f_hi}")
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        f_mid = f(mid)
        if abs(f_mid) < tol or (hi - lo) < tol:
            return mid
        if f_lo * f_mid <= 0:
            hi, f_hi = mid, f_mid
        else:
            lo, f_lo = mid, f_mid
    return 0.5 * (lo + hi)


if __name__ == "__main__":
    # quick self-check: an isentropic area ratio inverted should return the
    # same Mach number it was generated from, and likewise for Rayleigh flow.
    for M_test in [1.5, 2.0, 3.0, 5.0, 8.0]:
        AR = A_over_Astar(M_test)
        M_back = solve_M_from_area_ratio(AR)
        T0r = rayleigh_T0_over_T0star(M_test)
        M_back_ray = solve_M_from_rayleigh_T0ratio(T0r)
        print(f"M={M_test:4.1f}  A/A*={AR:7.4f} -> M_back={M_back:6.4f}   "
              f"T0/T0*={T0r:7.4f} -> M_back={M_back_ray:6.4f}")
