"""
main.py
-------
Command-line entry point: compute scramjet performance for a given
flight Mach number and altitude.

Examples
--------
    python3 main.py --M0 8 --altitude_km 30
    python3 main.py --M0 6 --altitude_km 25 --phi 0.6 --fuel H2
    python3 main.py --M0 5 --altitude_km 22 --fuel JP7 --phi 1.0
"""

import argparse
from scramjet_cycle import ScramjetInputs, run_cycle, print_report


def main():
    ap = argparse.ArgumentParser(description="1-D scramjet cycle performance model")
    ap.add_argument("--M0", type=float, default=8.0, help="flight Mach number")
    ap.add_argument("--altitude_km", type=float, default=30.0, help="flight altitude, km")
    ap.add_argument("--CR", type=float, default=4.0, help="inlet contraction ratio A0/A3")
    ap.add_argument("--phi", type=float, default=1.0, help="fuel equivalence ratio")
    ap.add_argument("--fuel", type=str, default="H2", choices=["H2", "JP7"], help="fuel type")
    ap.add_argument("--eta_b", type=float, default=0.90, help="combustion efficiency")
    ap.add_argument("--eta_e", type=float, default=0.95, help="nozzle expansion efficiency")
    args = ap.parse_args()

    inp = ScramjetInputs(
        M0=args.M0,
        altitude_m=args.altitude_km * 1000.0,
        CR=args.CR,
        equivalence_ratio=args.phi,
        fuel=args.fuel,
        eta_b=args.eta_b,
        eta_e=args.eta_e,
    )
    res = run_cycle(inp)
    print_report(res)


if __name__ == "__main__":
    main()
