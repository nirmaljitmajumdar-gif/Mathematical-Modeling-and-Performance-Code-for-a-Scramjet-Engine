"""
scramjet_cycle.py
------------------
1-D thermodynamic + gas-dynamic performance model of a scramjet engine.

Flow path modelled:  freestream (0) -> inlet -> combustor entry (3)
                      -> combustor (heat addition) -> combustor exit (4)
                      -> nozzle -> exit (10)

Method
------
* Inlet:      isentropic area-Mach relation (A0/A3 = specified contraction
              ratio CR) gives the ideal Mach-number drop from M0 to M3;
              the MIL-E-5008B empirical correlation gives the stagnation-
              pressure recovery PR(M0), which is the standard first-order
              way to fold inlet shock/viscous losses into a cycle-level
              model without solving the full oblique-shock chain.
* Combustor:  constant-area duct, Rayleigh flow. Heat addition from the
              fuel is computed from the equivalence ratio and combustion
              efficiency; the model checks this against the thermal-
              choking limit (T0/T0* <= 1) automatically.
* Nozzle:     isentropic expansion to the ambient (freestream) static
              pressure (perfectly-expanded-nozzle assumption), corrected
              by a nozzle thermodynamic efficiency eta_e.

All of these are the standard first-order simplifications used at the
cycle-analysis stage of a scramjet design (see accompanying
Scramjet_Engine_Mathematical_Model.docx for the underlying derivations
and literature sources) -- they trade some physical fidelity (e.g. no
explicit shock train, calorically-perfect gas, no dissociation) for a
transparent, fast, hand-checkable model, which is what a cycle code is
for. Assumptions are called out in the docstrings/comments below so they
can be stated explicitly in a project report.
"""

import math
from dataclasses import dataclass, field

from atmosphere import standard_atmosphere, R_AIR, G0
import gas_dynamics as gd

# ---------------------------------------------------------------------------
# Fuel data (lower heating value, stoichiometric fuel/air mass ratio)
# ---------------------------------------------------------------------------
FUELS = {
    "H2":   {"Hf": 120.0e6, "f_stoic": 0.0291},   # hydrogen
    "JP7":  {"Hf": 43.5e6,  "f_stoic": 0.0687},   # kerosene-type hydrocarbon
}


@dataclass
class ScramjetInputs:
    M0: float                      # flight Mach number
    altitude_m: float               # flight altitude, m
    CR: float = 4.0                 # inlet contraction ratio A0/A3
    equivalence_ratio: float = 1.0  # fuel equivalence ratio, phi
    fuel: str = "H2"
    eta_b: float = 0.90             # combustion efficiency
    eta_e: float = 0.95             # nozzle expansion efficiency
    gamma: float = 1.4              # calorically-perfect-gas assumption throughout
    cp: float = field(default=None) # derived from gamma, R if not given


@dataclass
class ScramjetResults:
    inputs: ScramjetInputs
    # freestream
    T_inf: float = None
    p_inf: float = None
    rho_inf: float = None
    V0: float = None
    T0_0: float = None
    p0_0: float = None
    # inlet / station 3
    PR: float = None
    M3: float = None
    T3: float = None
    p3: float = None
    T0_3: float = None
    p0_3: float = None
    # combustor / station 4
    f: float = None
    q_added: float = None
    choked: bool = False
    M4: float = None
    T4: float = None
    p4: float = None
    T0_4: float = None
    # nozzle / station 10
    M10_ideal: float = None
    T10_ideal: float = None
    T10: float = None
    V10: float = None
    # performance
    thrust_specific: float = None   # N per (kg/s of captured air)
    Isp_fuel_s: float = None        # standard aerospace Isp, based on fuel flow
    Isp_engine_s: float = None      # (V10-V0)/g -- ties back to cycle Eq. (38a)
    eta_p: float = None
    eta_th: float = None
    eta_o: float = None


def run_cycle(inp: ScramjetInputs) -> ScramjetResults:
    g = inp.gamma
    R = R_AIR
    cp = inp.cp if inp.cp else g * R / (g - 1.0)

    res = ScramjetResults(inputs=inp)

    # ---------------- 1. Freestream (station 0) ----------------
    T_inf, p_inf, rho_inf = standard_atmosphere(inp.altitude_m)
    V0 = inp.M0 * math.sqrt(g * R * T_inf)
    T0_0 = T_inf * gd.T0_over_T(inp.M0, g)
    p0_0 = p_inf * gd.p0_over_p(inp.M0, g)
    res.T_inf, res.p_inf, res.rho_inf, res.V0 = T_inf, p_inf, rho_inf, V0
    res.T0_0, res.p0_0 = T0_0, p0_0

    # ---------------- 2. Inlet: freestream -> station 3 ----------------
    # (a) Ideal isentropic Mach drop from the specified area (contraction) ratio.
    AR0 = gd.A_over_Astar(inp.M0, g)
    AR3 = AR0 / inp.CR
    M3 = gd.solve_M_from_area_ratio(AR3, g, supersonic=True)

    # (b) MIL-E-5008B empirical stagnation-pressure-recovery correlation,
    #     standard textbook benchmark for supersonic inlet losses (Hill &
    #     Peterson, "Mechanics and Thermodynamics of Propulsion"):
    if inp.M0 <= 1.0:
        PR = 1.0
    elif inp.M0 <= 5.0:
        PR = 1.0 - 0.075 * (inp.M0 - 1.0) ** 1.35
    else:
        PR = 800.0 / (inp.M0 ** 4 + 935.0)

    T0_3 = T0_0                      # adiabatic inlet: stagnation temperature conserved
    p0_3 = PR * p0_0
    T3 = T0_3 / gd.T0_over_T(M3, g)
    p3 = p0_3 / gd.p0_over_p(M3, g)

    res.PR, res.M3, res.T3, res.p3, res.T0_3, res.p0_3 = PR, M3, T3, p3, T0_3, p0_3

    # ---------------- 3. Combustor: station 3 -> station 4 (Rayleigh flow) ----------------
    fuel = FUELS[inp.fuel]
    f = inp.equivalence_ratio * fuel["f_stoic"]
    q_added = f * fuel["Hf"] * inp.eta_b          # heat released per unit AIR mass flow
    # energy balance on the (1+f) kg of combustion products per kg of air:
    T0_4_request = (T0_3 + q_added / cp) / 1.0     # simplified: heat raises T0 of the (~) same flow
    # (approximation: fuel mass addition to the flow's heat capacity neglected -- standard
    #  first-order simplification since f << 1 for realistic equivalence ratios)

    T0r_3 = gd.rayleigh_T0_over_T0star(M3, g)
    T0star = T0_3 / T0r_3
    T0r_4_request = T0_4_request / T0star

    choked = False
    if T0r_4_request >= 1.0:
        choked = True
        T0r_4 = 1.0
        M4 = 1.0
        T0_4 = T0star
    else:
        T0r_4 = T0r_4_request
        M4 = gd.solve_M_from_rayleigh_T0ratio(T0r_4, g, supersonic=True)
        T0_4 = T0_4_request

    p4 = p3 * gd.rayleigh_p_over_pstar(M4, g) / gd.rayleigh_p_over_pstar(M3, g)
    T4 = T0_4 / gd.T0_over_T(M4, g)

    res.f, res.q_added, res.choked = f, q_added, choked
    res.M4, res.T4, res.p4, res.T0_4 = M4, T4, p4, T0_4

    # ---------------- 4. Nozzle: station 4 -> station 10 (isentropic expansion) ----------------
    # perfectly-expanded-nozzle assumption: p10 = p_inf
    p0_4 = p4 * gd.p0_over_p(M4, g)
    p0_over_p_10 = p0_4 / p_inf
    M10_ideal = gd.solve_M_from_area_ratio  # placeholder to avoid confusion; not used directly

    # invert isentropic p0/p relation for the ideal exit Mach number
    def f_p(M):
        return gd.p0_over_p(M, g) - p0_over_p_10
    M10_ideal_val = gd._bisect(f_p, 1.0 + 1e-6, 60.0)
    T10_ideal = T0_4 / gd.T0_over_T(M10_ideal_val, g)

    # nozzle efficiency applied to the static-enthalpy drop (h4-h10)=eta_e*(h4-h10'):
    T10 = T4 - inp.eta_e * (T4 - T10_ideal)
    V4 = M4 * math.sqrt(g * R * T4)
    V10 = math.sqrt(max(V4 ** 2 + 2.0 * cp * (T4 - T10), 0.0))

    res.M10_ideal, res.T10_ideal, res.T10, res.V10 = M10_ideal_val, T10_ideal, T10, V10

    # ---------------- 5. Overall performance ----------------
    thrust_specific = (1.0 + f) * V10 - V0          # N per kg/s of captured air
    Isp_fuel = thrust_specific / (f * G0) if f > 0 else float("nan")
    Isp_engine = (V10 - V0) / G0

    eta_p = 2.0 * V0 / (V10 + V0)
    eta_th = ((1.0 + f) * V10 ** 2 - V0 ** 2) / (2.0 * f * fuel["Hf"]) if f > 0 else float("nan")
    eta_o = eta_p * eta_th

    res.thrust_specific = thrust_specific
    res.Isp_fuel_s = Isp_fuel
    res.Isp_engine_s = Isp_engine
    res.eta_p, res.eta_th, res.eta_o = eta_p, eta_th, eta_o

    return res


def print_report(res: ScramjetResults):
    inp = res.inputs
    print("=" * 72)
    print(f" SCRAMJET CYCLE PERFORMANCE — M0 = {inp.M0}, altitude = {inp.altitude_m/1000:.1f} km")
    print("=" * 72)
    print(f" Fuel: {inp.fuel}   phi = {inp.equivalence_ratio}   eta_b = {inp.eta_b}   eta_e = {inp.eta_e}   CR = {inp.CR}")
    print("-" * 72)
    print(f" Freestream:  T_inf={res.T_inf:8.2f} K   p_inf={res.p_inf:9.1f} Pa   V0={res.V0:8.1f} m/s")
    print(f" Station 3 (combustor entry):  PR={res.PR:.4f}   M3={res.M3:.3f}   T3={res.T3:8.2f} K   p3={res.p3:9.1f} Pa")
    choke_flag = "  *** THERMAL CHOKING LIMIT REACHED — heat input capped ***" if res.choked else ""
    print(f" Station 4 (combustor exit):   f={res.f:.4f}   M4={res.M4:.3f}   T4={res.T4:8.2f} K   p4={res.p4:9.1f} Pa{choke_flag}")
    print(f" Station 10 (nozzle exit):     M10_ideal={res.M10_ideal:.3f}   T10={res.T10:8.2f} K   V10={res.V10:8.1f} m/s")
    print("-" * 72)
    print(f" Specific thrust        : {res.thrust_specific:9.1f}  N per kg/s air")
    print(f" Specific impulse (fuel): {res.Isp_fuel_s:9.1f}  s")
    print(f" Specific impulse (eng.): {res.Isp_engine_s:9.1f}  s   (cycle Eq. 38a check)")
    print(f" Propulsive efficiency  : {res.eta_p*100:9.2f}  %")
    print(f" Thermal efficiency     : {res.eta_th*100:9.2f}  %")
    print(f" Overall efficiency     : {res.eta_o*100:9.2f}  %")
    print("=" * 72)


if __name__ == "__main__":
    inputs = ScramjetInputs(M0=8.0, altitude_m=30000.0, CR=4.0,
                             equivalence_ratio=1.0, fuel="H2",
                             eta_b=0.90, eta_e=0.95)
    results = run_cycle(inputs)
    print_report(results)
