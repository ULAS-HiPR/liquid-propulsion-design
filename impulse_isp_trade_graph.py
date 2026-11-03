"""
Isp trade study for Lúin V2 (Dagda, EuRoC 2027).

For a fixed total impulse, how much more propellant is needed if Isp drops?

    I_total = Isp * g0 * m_prop   ->   m_prop = I_total / (Isp * g0)

Edit the inputs below and run:  python impulse_isp_trade.py
"""

G0 = 9.80665  # m/s^2

# ---- Inputs ---------------------------------------------------------------
TOTAL_IMPULSE = 40_000.0   # N*s   (5 kN x 7.5 s baseline)
THRUST        = 5_000.0    # N     (held constant -> burn time fixed, mdot changes)
ISP_BASELINE  = 190.0      # s     effective (delivered) Isp
ISP_STEP      = 5.0        # s     step size for the sweep
N_STEPS       = 6          # number of steps below baseline to show

OF_RATIO      = 3.0        # oxidiser/fuel mass ratio (N2O/IPA) - update from RPA
RHO_N2O       = 786.0      # kg/m^3 liquid N2O (~20 C, saturated) - adjust for your tank temp
RHO_IPA       = 786.0      # kg/m^3 isopropanol
# ---------------------------------------------------------------------------


def propellant_mass(impulse, isp):
    """Total propellant mass [kg] for a given total impulse and Isp."""
    return impulse / (isp * G0)


def split(m_prop, of):
    """Split total propellant into oxidiser and fuel masses."""
    m_fuel = m_prop / (1.0 + of)
    m_ox = m_prop - m_fuel
    return m_ox, m_fuel


def main():
    m_base = propellant_mass(TOTAL_IMPULSE, ISP_BASELINE)
    ox_base, fuel_base = split(m_base, OF_RATIO)
    burn_time = TOTAL_IMPULSE / THRUST

    print(f"Total impulse : {TOTAL_IMPULSE:,.0f} N*s")
    print(f"Thrust        : {THRUST:,.0f} N  ->  burn time {burn_time:.2f} s")
    print(f"O/F           : {OF_RATIO:.2f}")
    print(f"Baseline Isp  : {ISP_BASELINE:.1f} s  ->  {m_base:.2f} kg propellant "
          f"({ox_base:.2f} kg N2O + {fuel_base:.2f} kg IPA)\n")

    header = (f"{'Isp [s]':>8} {'m_prop [kg]':>12} {'Δm [kg]':>9} {'Δm [%]':>8} "
              f"{'N2O [kg]':>9} {'IPA [kg]':>9} {'ΔN2O [L]':>9} {'ΔIPA [L]':>9} "
              f"{'mdot [kg/s]':>12}")
    print(header)
    print("-" * len(header))

    for i in range(N_STEPS + 1):
        isp = ISP_BASELINE - i * ISP_STEP
        m = propellant_mass(TOTAL_IMPULSE, isp)
        ox, fuel = split(m, OF_RATIO)
        dm = m - m_base
        d_ox_L = (ox - ox_base) / RHO_N2O * 1000
        d_fuel_L = (fuel - fuel_base) / RHO_IPA * 1000
        mdot = m / burn_time
        print(f"{isp:8.1f} {m:12.2f} {dm:9.2f} {100 * dm / m_base:8.2f} "
              f"{ox:9.2f} {fuel:9.2f} {d_ox_L:9.2f} {d_fuel_L:9.2f} {mdot:12.3f}")

    # Handy rule of thumb: sensitivity of propellant mass to Isp
    dm_per_s = m_base / ISP_BASELINE
    print(f"\nRule of thumb: each 1 s of Isp lost costs ~{dm_per_s:.3f} kg "
          f"({100 / ISP_BASELINE:.2f} %) of propellant at this Isp.")

    plot()


def plot():
    """Propellant mass vs Isp for the fixed total impulse."""
    import numpy as np
    import matplotlib.pyplot as plt

    isp_min = ISP_BASELINE - N_STEPS * ISP_STEP
    isp = np.linspace(isp_min, ISP_BASELINE + ISP_STEP, 200)
    m = propellant_mass(TOTAL_IMPULSE, isp)
    ox, fuel = split(m, OF_RATIO)

    series = [("Total", m, "#2a78d6"),
              ("N2O", ox, "#eb6834"),
              ("IPA", fuel, "#1baf7a")]

    fig, ax = plt.subplots(figsize=(8, 5))
    for name, y, colour in series:
        ax.plot(isp, y, color=colour, lw=2, label=name)
        ax.text(isp[0] - 0.5, y[0], name, color="#52514e",
                ha="right", va="center", fontsize=10)

    # Mark baseline and baseline - 5 s on the total line
    m_base = propellant_mass(TOTAL_IMPULSE, ISP_BASELINE)
    # baseline label below the line, the -5 s label above it so they don't clash
    for x, dy, va in ((ISP_BASELINE, -14, "top"), (ISP_BASELINE - 5, 14, "bottom")):
        y = propellant_mass(TOTAL_IMPULSE, x)
        ax.plot(x, y, "o", ms=8, color="#2a78d6",
                markeredgecolor="white", markeredgewidth=2, zorder=3)
        ax.annotate(f"{x:.0f} s\n{y:.2f} kg (+{y - m_base:.2f})",
                    (x, y), xytext=(0, dy), textcoords="offset points",
                    ha="center", va=va, fontsize=9, color="#0b0b0b")

    ax.set_xlabel("Effective Isp [s]")
    ax.set_ylabel("Propellant mass [kg]")
    ax.set_title(f"Propellant needed for {TOTAL_IMPULSE / 1000:.1f} kN·s "
                 f"(O/F = {OF_RATIO:g})", loc="left")
    ax.set_ylim(0, None)
    ax.set_xlim(isp_min - 6, isp[-1])
    ax.grid(True, color="#e6e5e0", lw=0.8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    fig.tight_layout()
    fig.savefig("isp_trade.png", dpi=150)
    plt.show()


if __name__ == "__main__":
    main()
