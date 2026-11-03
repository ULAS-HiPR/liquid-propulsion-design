#!/usr/bin/env python3
# Pintle injection geometry calculator using SPI, HEM & NHNE models
# Debugged version with fixes and improvements

import math
import CoolProp.CoolProp as CP

# --- Basic check that CoolProp works ---
print("If the number below is ~373.1243, CoolProp database is connected:")
try:
    print(CP.PropsSI('T', 'P', 101325, 'Q', 0, 'Water'))
except Exception as e:
    print("CoolProp test failed:", e)
    raise

# ----------------------
# INPUTS (user / case)
# ----------------------
T_C = 10.0              # inlet temperature in °C
p_inj = 5e6           # upstream N2O pressure [Pa]
p_cham = 2e6            # chamber/backpressure [Pa]
mdot_N2O = 0.736340877  # N2O mass flow [kg/s]
Re_guess = 1e6
Cd_N2O = 0.6
num_hole = 30

# ----------------------
# N2O properties at inlet (liquid, Q=0)
# ----------------------
T_K = T_C + 273.15
s_in = CP.PropsSI('S', 'T', T_K, 'Q', 0, 'NitrousOxide')
rho_inl = CP.PropsSI('D', 'Q', 0, 'T', T_K, 'NitrousOxide')
print("\nLiquid N2O density at inlet (kg/m^3):", rho_inl)

# ----------------------
# SPI (incompressible) model mass flux (kg / m^2 / s)
# ----------------------
deltaP = p_inj - p_cham
if deltaP <= 0:
    raise ValueError("Upstream pressure must exceed chamber pressure (deltaP <= 0).")

G_SPI = Cd_N2O * math.sqrt(2.0 * rho_inl * deltaP)
print("\nIncompressible model mass flux G_SPI (kg/m^2/s):", G_SPI)

# SPI orifice sizing
A_SPI_m2 = mdot_N2O / G_SPI
A_SPI_mm2 = A_SPI_m2 * 1e6
hole_area_SPI_m2 = A_SPI_m2 / num_hole
hole_radius_SPI_m = math.sqrt(hole_area_SPI_m2 / math.pi)
hole_radius_SPI_mm = hole_radius_SPI_m * 1e3
hole_diameter_SPI_mm = 2 * hole_radius_SPI_mm
print("Area of injector using SPI model: {:.6f} m^2  ({:.3f} mm^2)".format(A_SPI_m2, A_SPI_mm2))
print("SPI - Radius: {:.4f} mm, Diameter: {:.4f} mm (suggested diameter: {:.1f} mm)".format(
    hole_radius_SPI_mm, hole_diameter_SPI_mm, round(hole_diameter_SPI_mm, 1)))

# ----------------------
# HEM (homogeneous equilibrium) model
# ----------------------
h_in = CP.PropsSI('H', 'T', T_K, 'Q', 0, 'NitrousOxide')
h_out = CP.PropsSI('H', 'P', p_cham, 'S', s_in, 'NitrousOxide')
rho_out = CP.PropsSI('D', 'P', p_cham, 'S', s_in, 'NitrousOxide')

if (h_in - h_out) <= 0 or rho_out <= 0:
    raise RuntimeError("Non-physical HEM result: check states/inputs.")

G_HEM = Cd_N2O * rho_out * math.sqrt(2.0 * (h_in - h_out))
A_HEM_m2 = mdot_N2O / G_HEM
A_HEM_mm2 = A_HEM_m2 * 1e6
print("\nArea of injector using HEM model: {:.6f} m^2  ({:.3f} mm^2) at {} °C".format(
    A_HEM_m2, A_HEM_mm2, T_C))

# HEM orifice sizing
hole_area_HEM_m2 = A_HEM_m2 / num_hole
hole_radius_HEM_m = math.sqrt(hole_area_HEM_m2 / math.pi)
hole_radius_HEM_mm = hole_radius_HEM_m * 1e3
hole_diameter_HEM_mm = 2 * hole_radius_HEM_mm
print("HEM - Radius: {:.4f} mm, Diameter: {:.4f} mm (suggested diameter: {:.1f} mm)".format(
    hole_radius_HEM_mm, hole_diameter_HEM_mm, round(hole_diameter_HEM_mm, 1)))

# ----------------------
# NHNE (non-homogeneous non-equilibrium) estimate
# ----------------------
L = 0.0015  # characteristic length [m]
p_sat = CP.PropsSI('P', 'T', T_K, 'Q', 0, 'NitrousOxide')

tau_b = math.sqrt((3.0 / 2.0) * rho_inl / deltaP)
tau_r = L * math.sqrt(rho_inl / (2.0 * deltaP))
kappa = 1.0 / (1.0 + (tau_b / tau_r))

G_NHNEest = 0.5 * (G_SPI + G_HEM)
inj_totalarea = mdot_N2O / G_NHNEest
mdot_NHNE = 0.0
tolerance = 1e-8
step_size = max(inj_totalarea * 1e-5, 1e-12)
max_iters = 200000

areas = []

for it in range(max_iters):
    G_NHNE = (kappa / (1.0 + kappa)) * G_SPI + (1.0 / (1.0 + kappa)) * G_HEM
    mdot_NHNE = G_NHNE * inj_totalarea
    
    if abs(mdot_N2O - mdot_NHNE) < tolerance:
        break
    
    if mdot_NHNE < mdot_N2O:
        inj_totalarea += step_size
    else:
        inj_totalarea = max(inj_totalarea - step_size, 1e-12)

    areas.append(inj_totalarea)
else:
    print("Warning: NHNE iteration did not converge within max iterations.")

inj_totalarea_mm2 = inj_totalarea * 1e6
print("\nFINAL injector total area (NHNE estimate): {:.6f} m^2  ({:.3f} mm^2)".format(
    inj_totalarea, inj_totalarea_mm2))

hole_area_m2 = inj_totalarea / num_hole
hole_radius_m = math.sqrt(hole_area_m2 / math.pi)
hole_radius_mm = hole_radius_m * 1e3
hole_diameter_mm = 2 * hole_radius_mm
print("Number of holes:", num_hole)
print("Single hole area: {:.6e} m^2  ({:.3f} mm^2)".format(hole_area_m2, hole_area_m2 * 1e6))
print("NHNE - Radius: {:.4f} mm, Diameter: {:.4f} mm (suggested diameter: {:.1f} mm)".format(
    hole_radius_mm, hole_diameter_mm, round(hole_diameter_mm, 1)))

# ----------------------
# IPA (isopropanol) annulus - using incompressible SPI model
# ----------------------
mdot_IPA = 0.245446959  # kg/s
p_in_IPA = 2.5e6
T_IPA_C = T_C + 10  # BUG FIX: IPA is 10°C warmer than N2O
T_IPA_K = T_IPA_C + 273.15

fluid_IPA = 'Isopropanol'
try:
    # BUG FIX: Use T and P (not P,T order) for consistency, specify liquid phase
    rho_in_IPA = CP.PropsSI('D', 'T', T_IPA_K, 'P', p_in_IPA, fluid_IPA)
except Exception:
    fluid_IPA = 'Ethanol'
    rho_in_IPA = CP.PropsSI('D', 'T', T_IPA_K, 'P', p_in_IPA, fluid_IPA)

print("\nLiquid {} density at inlet (kg/m^3): {:.3f}".format(fluid_IPA, rho_in_IPA))

Cd_IPA = 0.9
Annulus_deltaP = p_in_IPA - p_cham
if Annulus_deltaP <= 0:
    raise ValueError("IPA upstream pressure must exceed chamber pressure (Annulus_deltaP <= 0).")

G_SPI_IPA = Cd_IPA * math.sqrt(2.0 * rho_in_IPA * Annulus_deltaP)
A_IPA_m2 = mdot_IPA / G_SPI_IPA
print("Area of annulus: {:.6f} m^2  ({:.3f} mm^2)".format(A_IPA_m2, A_IPA_m2 * 1e6))

# ----------------------
# Annular geometry calculation
# ----------------------
radius_N2O = 0.0075  # 7.5 mm radius for pintle body
annular_radius = math.sqrt((A_IPA_m2 + (math.pi * radius_N2O ** 2)) / math.pi)
annular_width = annular_radius - radius_N2O
pintle_diameter = 2 * radius_N2O
annulus_diameter = 2 * annular_radius

print("\nAnnular radius (from center): {:.3f} mm".format(annular_radius * 1e3))
print("Annular width: {:.3f} mm".format(annular_width * 1e3))
print("Radius of pintle (to hole centroid): {:.3f} mm".format(radius_N2O * 1e3))
print("Diameter of pintle body: {:.3f} mm".format(pintle_diameter * 1e3))
print("Diameter of annulus (outer): {:.3f} mm".format(annulus_diameter * 1e3))

# ----------------------
# Film cooling holes (IPA)
# ----------------------
film_cooling_percent = 10.0
A_film_IPA = (mdot_IPA * (film_cooling_percent / 100.0)) / G_SPI_IPA
print("\n--- Film Cooling ---")
print("Film cooling percentage: {:.0f}%".format(film_cooling_percent))
print("Area of film holes (m^2 / mm^2): {:.6e} m^2  ({:.3f} mm^2)".format(A_film_IPA, A_film_IPA * 1e6))
film_num_hole = 20
film_hole_area = A_film_IPA / film_num_hole
film_hole_radius_m = math.sqrt(film_hole_area / math.pi)
film_hole_diameter_m = 2 * film_hole_radius_m
print("IPA film cooling holes - Radius: {:.4f} mm, Diameter: {:.4f} mm (suggested diameter: {:.1f} mm)".format(
    film_hole_radius_m * 1e3, film_hole_diameter_m * 1e3, round(film_hole_diameter_m * 1e3, 1)))

print("\nN2O outlet density rho_out (HEM, kg/m^3): {:.3f}".format(rho_out))
# ----------------------
# Spray angle estimate
# ----------------------
# BUG FIX: Use rho_in_IPA (inlet density) for velocity calculation
v_IPA = Cd_IPA * math.sqrt((2.0 * Annulus_deltaP) / rho_in_IPA)
v_N2O = Cd_N2O * math.sqrt((2.0 * deltaP) / rho_out)

if mdot_IPA <= 0 or v_IPA <= 0:
    print("\nCannot compute TMR/spray angle: zero IPA mass flow or velocity.")
else:
    TMR = (mdot_N2O * v_N2O) / (mdot_IPA * v_IPA)
    acos_arg = 1.0 / (1.0 + TMR) if (1.0 + TMR) != 0 else 1.0
    acos_arg = max(min(acos_arg, 1.0), -1.0)
    rad = math.acos(acos_arg)
    theta_deg = rad * 180.0 / math.pi
    print("\nSpray angle of atomised mixture (degrees): {:.2f}".format(theta_deg))

print("\n-- END OF PROGRAM --")
