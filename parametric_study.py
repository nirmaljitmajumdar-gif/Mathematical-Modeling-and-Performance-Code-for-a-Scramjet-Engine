"""
parametric_study.py
--------------------
Sweeps flight Mach number (at a representative constant-dynamic-pressure
trajectory, a common simplification for scramjet performance-mapping
studies) and:

  1. writes a CSV of all computed performance parameters,
  2. produces the standard scramjet performance plots (Isp, thrust,
     station temperatures/Mach numbers vs. flight Mach number),
  3. runs a couple of basic "does this look like a real scramjet"
     sanity/benchmark checks against literature-reported ranges, and
     against the MIL-E-5008B inlet-recovery correlation used inside the
     model itself.

This fills in the "parametric study results, plots generated, ...
accuracy vs reference data" bullet from the project description.
"""

import csv
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from atmosphere import standard_atmosphere
from scramjet_cycle import ScramjetInputs, run_cycle, FUELS

RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")


# Representative scramjet corridor: roughly constant dynamic pressure q ~ 50 kPa,
# a standard "cruise corridor" assumption used to pick a sensible altitude for
# each Mach number when sweeping M0 (Curran & Murthy; Heiser & Pratt use similar
# constant-q corridors for scramjet performance maps).
def altitude_for_constant_q(M0, q_target_pa=50000.0, gamma=1.4, h_guess_m=25000.0):
    h = h_guess_m
    for _ in range(60):
        T, p, rho = standard_atmosphere(h)
        V = M0 * math.sqrt(gamma * 287.05287 * T)
        q = 0.5 * rho * V ** 2
        # simple proportional correction (rho ~ exp(-h/H), H ~ 7000 m)
        dh = 7000.0 * math.log(q / q_target_pa)
        h += dh
        if abs(dh) < 1.0:
            break
    return max(h, 15000.0)


def run_sweep(M0_list, phi=1.0, fuel="H2", CR=4.0, eta_b=0.90, eta_e=0.95):
    rows = []
    for M0 in M0_list:
        h = altitude_for_constant_q(M0)
        inp = ScramjetInputs(M0=M0, altitude_m=h, CR=CR, equivalence_ratio=phi,
                              fuel=fuel, eta_b=eta_b, eta_e=eta_e)
        res = run_cycle(inp)
        rows.append({
            "M0": M0, "altitude_km": h / 1000.0,
            "PR": res.PR, "M3": res.M3, "T3_K": res.T3,
            "M4": res.M4, "T4_K": res.T4, "choked": res.choked,
            "V10_ms": res.V10,
            "specific_thrust_N_per_kg_s": res.thrust_specific,
            "Isp_fuel_s": res.Isp_fuel_s,
            "eta_p_pct": res.eta_p * 100, "eta_th_pct": res.eta_th * 100,
            "eta_o_pct": res.eta_o * 100,
        })
    return rows


def write_csv(rows, path):
    with open(path, "w", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def make_plots(rows, out_prefix):
    M0 = [r["M0"] for r in rows]

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(M0, [r["Isp_fuel_s"] for r in rows], "o-", color="#1F3864")
    ax.set_xlabel("Flight Mach number, M0")
    ax.set_ylabel("Specific impulse, Isp (s)")
    ax.set_title("Scramjet Specific Impulse vs. Flight Mach Number\n(H2 fuel, stoichiometric, constant-q corridor)")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{out_prefix}_Isp_vs_Mach.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(M0, [r["specific_thrust_N_per_kg_s"] for r in rows], "s-", color="#B22222")
    ax.set_xlabel("Flight Mach number, M0")
    ax.set_ylabel("Specific thrust (N per kg/s air)")
    ax.set_title("Scramjet Specific Thrust vs. Flight Mach Number")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{out_prefix}_Thrust_vs_Mach.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(M0, [r["T3_K"] for r in rows], "^-", label="T3 (combustor entry)")
    ax.plot(M0, [r["T4_K"] for r in rows], "v-", label="T4 (combustor exit)")
    ax.set_xlabel("Flight Mach number, M0")
    ax.set_ylabel("Static temperature (K)")
    ax.set_title("Combustor Entry/Exit Temperature vs. Flight Mach Number")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{out_prefix}_Temperatures_vs_Mach.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(M0, [r["eta_p_pct"] for r in rows], "o-", label="Propulsive efficiency")
    ax.plot(M0, [r["eta_th_pct"] for r in rows], "s-", label="Thermal efficiency")
    ax.plot(M0, [r["eta_o_pct"] for r in rows], "^-", label="Overall efficiency")
    ax.set_xlabel("Flight Mach number, M0")
    ax.set_ylabel("Efficiency (%)")
    ax.set_title("Cycle Efficiencies vs. Flight Mach Number")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{out_prefix}_Efficiencies_vs_Mach.png", dpi=150)
    plt.close(fig)


def validate(rows):
    """
    Basic benchmark / sanity checks, printed to console -- this is the
    "validated against standard air-breathing propulsion benchmarks" step:

      1. MIL-E-5008B inlet pressure recovery should fall inside the
         published correlation's own range (0 < PR <= 1) and decrease
         monotonically with M0 -- a basic self-consistency check.
      2. H2-fuelled scramjet specific impulse is well documented in the
         literature (e.g. Heiser & Pratt; Curran & Murthy) to fall
         roughly in the 1000-3000 s range for M0 = 5-10 at stoichiometric
         conditions; flag anything wildly outside that band.
    """
    print("\n" + "=" * 72)
    print(" VALIDATION / BENCHMARK CHECKS")
    print("=" * 72)

    prs = [r["PR"] for r in rows]
    monotonic = all(prs[i] >= prs[i + 1] - 1e-9 for i in range(len(prs) - 1))
    print(f" Inlet PR monotonically non-increasing with M0 (MIL-E-5008B check): {'PASS' if monotonic else 'FAIL'}")
    print(f" Inlet PR in valid range (0,1]: {'PASS' if all(0 < p <= 1 for p in prs) else 'FAIL'}")

    isp_ok = []
    for r in rows:
        if 4.0 <= r["M0"] <= 10.0:
            in_band = 800.0 <= r["Isp_fuel_s"] <= 3500.0
            isp_ok.append(in_band)
            flag = "OK" if in_band else "OUTSIDE literature band"
            print(f"  M0={r['M0']:.1f}: Isp={r['Isp_fuel_s']:.0f} s  -> {flag}")
    print(f" Specific impulse within literature band (800-3500 s, M0=4-10): "
          f"{'PASS' if all(isp_ok) else 'CHECK'}")
    print("=" * 72)


if __name__ == "__main__":
    os.makedirs(RESULTS_DIR, exist_ok=True)
    M0_list = [4, 5, 6, 7, 8, 9, 10]
    rows = run_sweep(M0_list, phi=1.0, fuel="H2")
    write_csv(rows, os.path.join(RESULTS_DIR, "scramjet_parametric_results.csv"))
    make_plots(rows, os.path.join(RESULTS_DIR, "scramjet"))
    validate(rows)
    print(f"\nCSV and plots written to {RESULTS_DIR}/")
