"""
atmosphere.py
-------------
International Standard Atmosphere (ISA), 1976 model, valid 0-47 km.
Used to get freestream static temperature, pressure and density at a
given flight altitude for the scramjet cycle code.
"""

import math

G0 = 9.80665        # m/s^2, standard gravity
R_AIR = 287.05287    # J/(kg.K), specific gas constant for air

# Layer base data: (base altitude [m], base temperature [K], lapse rate [K/m], base pressure [Pa])
T0 = 288.15
P0 = 101325.0

_LAYERS = [
    (0.0,     288.15, -0.0065),
    (11000.0, 216.65,  0.0),
    (20000.0, 216.65,  0.0010),
    (32000.0, 228.65,  0.0028),
    (47000.0, 270.65,  0.0),
]


def _layer_pressures():
    """Pre-compute base pressure at the start of every layer by marching up from sea level."""
    pressures = [P0]
    for i in range(1, len(_LAYERS)):
        h0, T_b0, L0 = _LAYERS[i - 1]
        h1, T_b1, _ = _LAYERS[i]
        p0 = pressures[-1]
        if abs(L0) < 1e-12:
            p1 = p0 * math.exp(-G0 * (h1 - h0) / (R_AIR * T_b0))
        else:
            p1 = p0 * (T_b1 / T_b0) ** (-G0 / (L0 * R_AIR))
        pressures.append(p1)
    return pressures


_BASE_PRESSURES = _layer_pressures()


def standard_atmosphere(h_m):
    """
    Return (T [K], p [Pa], rho [kg/m^3]) at geopotential altitude h_m (metres).
    Valid 0 - 47000 m (covers the whole scramjet flight corridor, ~20-40 km).
    """
    if h_m < 0:
        h_m = 0.0
    if h_m > 47000.0:
        h_m = 47000.0  # clamp; scramjet cruise is well inside this range

    # find layer
    idx = 0
    for i in range(len(_LAYERS) - 1):
        if h_m >= _LAYERS[i][0]:
            idx = i
    h_b, T_b, L = _LAYERS[idx]
    p_b = _BASE_PRESSURES[idx]

    T = T_b + L * (h_m - h_b)
    if abs(L) < 1e-12:
        p = p_b * math.exp(-G0 * (h_m - h_b) / (R_AIR * T_b))
    else:
        p = p_b * (T / T_b) ** (-G0 / (L * R_AIR))

    rho = p / (R_AIR * T)
    return T, p, rho


if __name__ == "__main__":
    for h_km in [0, 11, 20, 25, 30, 35, 40]:
        T, p, rho = standard_atmosphere(h_km * 1000)
        print(f"h={h_km:>3} km   T={T:7.2f} K   p={p:10.1f} Pa   rho={rho:.5f} kg/m^3")
