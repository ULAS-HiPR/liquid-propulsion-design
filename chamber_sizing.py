"""
Combustion chamber sizing 
Sizes the chamber from the injector face to the throat, plus the nozzle
expansion ratio (exit area / throat area).

Inputs you need to have ;Isp, Cstar, Gamma (all from CEA run), Lstar (empirical) and Contraction ratio (empirical, H&H)

Method:
    mdot = F / (Isp * g0)                 mass flow from delivered Isp
    At   = mdot * c* / Pc                 throat area (choked-flow / c* definition)
    Vc   = L* * At                        chamber volume (incl. convergent)
    Ac   = (Ac/At) * At                   contraction ratio is an input
    Lcyl = (Vc - Vconv) / Ac              cylinder length after convergent frustum
    Me   = from isentropic Pe/Pc          exit Mach number
    Ae/At = (1/Me) * [2/(g+1) * (1 + (g-1)/2 * Me^2)]^((g+1)/(2(g-1)))

-Remember the assumptions you are making. Real flow is not isentropic nor perfect combustion efficiency. Apply correction factors to the ideal numbers you get
-Know the limitations of the 1D isentropic relations & area-mach relations.
"""
import math

G0 = 9.80665  # m/s^2

# ---- Inputs ---------------------------------------------------------------
F     = 5000.0    # N      thrust
ISP   = 185.0     # s      delivered Isp
PC    = 20e5      # Pa     chamber pressure
GAMMA = 1.245      # -      ratio of specific heats (chamber)
CSTAR = 1465   # m/s    characteristic velocity (use delivered c* if you have it)
LSTAR = 0.8       # m      characteristic length
EPS_C = 4.2       # -      contraction ratio Ac/At
PE    = 1.01325e5 # Pa     nozzle exit pressure (design expansion)

BETA_DEG = 30.0   # deg    convergent half-angle
# ---------------------------------------------------------------------------


def area_ratio(M, g):
    """A/At for a given Mach number (isentropic)."""
    return (1 / M) * ((2 / (g + 1)) * (1 + (g - 1) / 2 * M**2)) ** ((g + 1) / (2 * (g - 1)))


def subsonic_mach(eps, g):
    """Solve A/At = eps for the subsonic Mach number (bisection)."""
    lo, hi = 1e-6, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if area_ratio(mid, g) > eps:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main():
    # 1. Mass flow
    mdot = F / (ISP * G0)

    # 2. Throat
    At = mdot * CSTAR / PC
    Dt = math.sqrt(4 * At / math.pi)
    Rt = Dt / 2

    # 3. Chamber volume
    Vc = LSTAR * At

    # 4. Chamber diameter from contraction ratio
    eps_c = EPS_C
    Ac = eps_c * At
    Dc = math.sqrt(4 * Ac / math.pi)
    Rc = Dc / 2
    Mc = subsonic_mach(eps_c, GAMMA)

    # 5. Convergent section (frustum)
    beta = math.radians(BETA_DEG)
    L_conv = (Rc - Rt) / math.tan(beta)
    V_conv = math.pi / 3 * L_conv * (Rc**2 + Rc * Rt + Rt**2)

    # 6. Cylindrical section
    L_cyl = (Vc - V_conv) / Ac
    L_tot = L_cyl + L_conv

    # 7. Nozzle expansion ratio from Pe/Pc (isentropic)
    g = GAMMA
    Me = math.sqrt(2 / (g - 1) * ((PC / PE) ** ((g - 1) / g) - 1))
    eps_e = area_ratio(Me, g)
    Ae = eps_e * At
    De = math.sqrt(4 * Ae / math.pi)

    print("Inputs")
    print(f"  F = {F:.0f} N, Isp = {ISP:.1f} s, Pc = {PC / 1e5:.1f} bar, "
          f"Pe = {PE / 1e5:.3f} bar")
    print(f"  gamma = {GAMMA:.3f}, c* = {CSTAR:.0f} m/s, L* = {LSTAR:.2f} m, "
          f"Ac/At = {EPS_C:.2f}, beta = {BETA_DEG:.0f} deg\n")

    print("Results")
    print(f"  Mass flow          mdot   = {mdot:8.3f} kg/s")
    print(f"  Throat area        At     = {At * 1e6:8.1f} mm^2")
    print(f"  Throat diameter    Dt     = {Dt * 1e3:8.2f} mm")
    print(f"  Chamber volume     Vc     = {Vc * 1e6:8.0f} cm^3")
    print(f"  Contraction ratio  Ac/At  = {eps_c:8.2f}")
    print(f"  Chamber diameter   Dc     = {Dc * 1e3:8.2f} mm")
    print(f"  Chamber Mach       Mc     = {Mc:8.3f}")
    print(f"  Convergent length  Lconv  = {L_conv * 1e3:8.1f} mm")
    print(f"  Convergent volume  Vconv  = {V_conv * 1e6:8.0f} cm^3")
    print(f"  Cylinder length    Lcyl   = {L_cyl * 1e3:8.1f} mm")
    print(f"  Face-to-throat     Ltot   = {L_tot * 1e3:8.1f} mm")
    print(f"  Lcyl / Dc                 = {L_cyl / Dc:8.2f}")
    print(f"  Exit Mach          Me     = {Me:8.3f}")
    print(f"  Expansion ratio    Ae/At  = {eps_e:8.2f}")
    print(f"  Exit area          Ae     = {Ae * 1e6:8.1f} mm^2")
    print(f"  Exit diameter      De     = {De * 1e3:8.2f} mm")

    if L_cyl <= 0:
        print("\n  WARNING: convergent volume exceeds Vc - increase L* or beta.")
    if Mc > 0.3:
        print("\n  WARNING: chamber Mach > 0.3 - expect Rayleigh pressure loss.")


if __name__ == "__main__":
    main()
