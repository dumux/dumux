#!/usr/bin/env python3
"""
Analytical solution for pipe-flow BHE heat transport — two formulations compared.

Formulation A (OGS original): U referred to r_po, characteristic length X uses r_pi (inconsistent)
Formulation B (consistent):   U and X both referred to r_pi — Ramey (1962)
"""
import glob
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ── Input ──────────────────────────────────────────────────────────────────────
dumux_output_name = "test_wellbore_heat_transport_1d"  # DuMux output base name (csv/vtp prefix)

# ── Physical parameters (OGS benchmark) ───────────────────────────────────────
T_d = 55            # undisturbed formation temperature, °C
T_i = 20            # inlet fluid temperature, °C
Z = np.insert(np.linspace(0.3, 30, 100), 0, 1e-10)  # depth, m
q = 0.0002          # flow rate, m³/s
rho_f = 1000        # fluid density, kg/m³
c_p_f = 4190        # fluid specific heat, J/kg/K
mu_f = 1.14e-3      # dynamic viscosity, Pa·s
lambda_f = 0.59     # fluid thermal conductivity, W/m/K
rho_re = 1800       # formation density, kg/m³
c_p_re = 1778       # formation specific heat, J/kg/K
lambda_re = 2.78018 # formation thermal conductivity, W/m/K
lambda_g = 0.73     # grout thermal conductivity, W/m/K
lambda_pi = 1.3     # pipe wall thermal conductivity, W/m/K
r_pi = 0.12913      # pipe inner radius, m
r_b = 0.14          # borehole radius, m
t_pi = 0.00587      # pipe wall thickness, m
t = 86400 * 5       # simulation time, s

# ── Helper functions ───────────────────────────────────────────────────────────

def outlet_temp(T_d, T_i, z, X):
    return T_d + (T_i - T_d) * np.exp(-z / X)

def coefficient_x(q, rho_f, c_p_f, lambda_re, r_ref, U, f_t):
    return (q * rho_f * c_p_f) * (lambda_re + r_ref * U * f_t) / (2 * np.pi * r_ref * U * lambda_re)

def dimensionless_time(lambda_re, delta_t, rho_re, c_p_re, r_b):
    return lambda_re * delta_t / (rho_re * c_p_re * r_b**2)

def time_function(t_D):
    if t_D > 1.5:
        return (0.4063 + 0.5 * np.log(t_D)) * (1 + 0.6 / t_D)
    return 1.1281 * np.sqrt(t_D) * (1 - 0.3 * np.sqrt(t_D))

# ── Nusselt number and convective heat transfer coefficient ────────────────────
v = q / (np.pi * r_pi**2)
Pr = mu_f * c_p_f / lambda_f
Re = rho_f * v * (2 * r_pi) / mu_f

if Re < 2300:
    Nu_p = 4.364
elif Re > 10000:
    xi = (1.8 * np.log(Re) - 1.5) ** -2
    Nu_p = (xi / 8 * Re * Pr) / (1 + 12.7 * np.sqrt(xi / 8) * (Pr**(2/3) - 1)) * (1 + (r_pi / Z)**(2/3))
else:
    gamma = (Re - 2300) / (10000 - 2300)
    xi_t = 0.0308  # (1.8*log10(1e4) - 1.5)^-2
    Nu_t = (xi_t / 8 * 1e4 * Pr) / (1 + 12.7 * np.sqrt(xi_t / 8) * (Pr**(2/3) - 1)) * (1 + (r_pi / Z)**(2/3))
    Nu_p = (1 - gamma) * 4.364 + gamma * Nu_t

h = lambda_f * Nu_p / (2 * r_pi)

# ── Overall heat transfer coefficients ────────────────────────────────────────
r_po = r_pi + t_pi

# Formulation A: U referred to r_po, X uses r_pi (inconsistent)
U_A = 1 / (r_po / (r_pi * h) + r_po * (np.log(r_po / r_pi) / lambda_pi + np.log(r_b / r_po) / lambda_g))

# Formulation B: U and X both referred to r_pi (consistent)
U_B = 1 / (1 / h + r_pi * (np.log(r_po / r_pi) / lambda_pi + np.log(r_b / r_po) / lambda_g))

print(f"Re={Re:.1f}, Pr={Pr:.3f}, Nu={Nu_p:.4f}, h={h:.4f} W/m²K")
print(f"  U_A (ref r_po={r_po:.5f} m) = {U_A:.4f} W/m²K")
print(f"  U_B (ref r_pi={r_pi:.5f} m) = {U_B:.4f} W/m²K")

# ── Compute analytical solutions ───────────────────────────────────────────────
time_range = range(7200, t + 1, 7200)

data_A = pd.DataFrame(index=Z, columns=list(time_range))
data_B = pd.DataFrame(index=Z, columns=list(time_range))

for delta_z in Z:
    T_A, T_B = [], []
    for delta_t in time_range:
        t_D = dimensionless_time(lambda_re, delta_t, rho_re, c_p_re, r_b)
        f_t = time_function(t_D)
        T_A.append(outlet_temp(T_d, T_i, delta_z, coefficient_x(q, rho_f, c_p_f, lambda_re, r_pi, U_A, f_t)))
        T_B.append(outlet_temp(T_d, T_i, delta_z, coefficient_x(q, rho_f, c_p_f, lambda_re, r_pi, U_B, f_t)))
    data_A.loc[delta_z] = T_A
    data_B.loc[delta_z] = T_B

# ── Load DuMux output ──────────────────────────────────────────────────────────
script_dir = os.path.dirname(os.path.abspath(__file__))
csv_name = dumux_output_name + ".csv"

dumux_time_d = dumux_T = None
for path in [csv_name, os.path.join(script_dir, csv_name)]:
    if os.path.exists(path):
        raw = np.genfromtxt(path, delimiter=",", skip_header=1)
        if raw.ndim == 1:
            raw = raw[np.newaxis, :]
        dumux_time_d = raw[:, 0]
        dumux_T = raw[:, 1]
        print(f"Loaded DuMux output: {path}  ({len(dumux_time_d)} points)")
        break
else:
    print(f"'{csv_name}' not found — run the simulation first.")

# ── Load DuMux temperature profile along depth (final-timestep 1d VTP) ─────────
def load_dumux_profile(prefix, search_dirs):
    """Return (x, T[°C]) along the pipe at the final timestep, from the 1d VTP output."""
    try:
        import vtk
        from vtk.util.numpy_support import vtk_to_numpy
    except ImportError:
        print("python3-vtk not available — skipping temperature-distribution plot.")
        return None, None

    vtp = None
    for d in search_dirs:
        files = glob.glob(os.path.join(d, prefix + "-*.vtp"))
        if files:
            vtp = max(files, key=lambda f: int(re.search(r"-(\d+)\.vtp$", f).group(1)))
            break
    if vtp is None:
        print(f"'{prefix}-*.vtp' not found — skipping temperature-distribution plot.")
        return None, None

    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(vtp)
    reader.Update()
    poly = reader.GetOutput()
    x = vtk_to_numpy(poly.GetPoints().GetData())[:, 0]
    T = vtk_to_numpy(poly.GetPointData().GetArray("T")) - 273.15  # K -> °C
    order = np.argsort(x)
    print(f"Loaded DuMux profile: {vtp}  ({len(x)} points)")
    return x[order], T[order]

dumux_x, dumux_T_profile = load_dumux_profile(dumux_output_name, [script_dir, "."])

# ── OGS reference data ─────────────────────────────────────────────────────────
ogs_time_s = np.linspace(0, 5 * 86400, 31)
ogs_T = np.array([
    20.        , 25.67074983, 26.67948325, 26.54800871, 26.27755879, 26.03968084,
    25.84599915, 25.68828698, 25.55846907, 25.45030474, 25.35907264, 25.28119084,
    25.21392656, 25.15518214, 25.10333715, 25.05713156, 25.01557889, 24.97790140,
    24.94348150, 24.91182514, 24.88253415, 24.85528518, 24.82981365, 24.80590137,
    24.78336707, 24.76205898, 24.74184908, 24.72262861, 24.70430449, 24.68679654,
    24.67003519,
])
ogs_error_A = np.array([
    -7.74064046, -0.56492394,  0.71049789,  0.76404746,  0.67145754,  0.56154154,
     0.46794475,  0.39242137,  0.33228395,  0.28454449,  0.24660294,  0.21634100,
     0.19208005,  0.17250760,  0.15660424,  0.14358037,  0.13282477,  0.12386385,
     0.11632997,  0.10993699,  0.10446154,  0.09972859,  0.09560050,  0.09196854,
     0.08874648,  0.08586559,  0.08327079,  0.08091778,  0.07877067,  0.07680025,
     0.07498260,
])

analytical_B_at_ogs = np.array([
    outlet_temp(T_d, T_i, Z[-1],
                coefficient_x(q, rho_f, c_p_f, lambda_re, r_pi, U_B,
                              time_function(dimensionless_time(lambda_re, dt, rho_re, c_p_re, r_b))))
    for dt in ogs_time_s
])
ogs_error_B = ogs_T - analytical_B_at_ogs

# OGS distributed fluid temperature along depth at the final time (t = 5 d),
# sampled at 101 points from the inlet (0 m) to the bottom (30 m)
ogs_x_profile = np.linspace(0, 30, 101)
ogs_T_profile = np.array([
    20.        , 20.0499748 , 20.0999496 , 20.1499244 , 20.19959088, 20.24910318,
    20.29861549, 20.34812923, 20.39764582, 20.44716241, 20.49667901, 20.5457412 ,
    20.59480339, 20.64386559, 20.69293331, 20.74200381, 20.7910743 , 20.83999301,
    20.88860817, 20.93722332, 20.98583847, 21.03446673, 21.08309499, 21.13172325,
    21.18004725, 21.22821912, 21.276391  , 21.32456885, 21.37275866, 21.42094847,
    21.46913828, 21.5168706 , 21.56460292, 21.61233525, 21.66008277, 21.70783788,
    21.755593  , 21.80319523, 21.85049171, 21.89778819, 21.94508467, 21.99240882,
    22.03973297, 22.08705711, 22.1340747 , 22.18093901, 22.22780331, 22.27467848,
    22.32157535, 22.36847222, 22.4153691 , 22.46180488, 22.50824065, 22.55467643,
    22.6011372 , 22.64761047, 22.69408374, 22.74040287, 22.78641374, 22.8324246 ,
    22.87843546, 22.92448876, 22.97054207, 23.01659537, 23.06233949, 23.10792902,
    23.15351855, 23.19912389, 23.24476083, 23.29039778, 23.33603473, 23.38120648,
    23.42637824, 23.47154999, 23.51675669, 23.56198086, 23.60720503, 23.65227365,
    23.69703115, 23.74178866, 23.78654617, 23.83136112, 23.87617606, 23.92099101,
    23.96549383, 24.00984059, 24.05418734, 24.09855493, 24.14296419, 24.18737345,
    24.2317827 , 24.27572213, 24.31966155, 24.36360098, 24.40758592, 24.45159361,
    24.4956013 , 24.5394493 , 24.58297793, 24.62650656, 24.67003519,
])

# ── DuMux error vs Formulation B ──────────────────────────────────────────────
dumux_time_s = dumux_error = None
if dumux_time_d is not None:
    dumux_time_s = dumux_time_d * 86400
    analytical_at_dumux = np.array([
        outlet_temp(T_d, T_i, Z[-1],
                    coefficient_x(q, rho_f, c_p_f, lambda_re, r_pi, U_B,
                                  time_function(dimensionless_time(lambda_re, dt, rho_re, c_p_re, r_b))))
        for dt in dumux_time_s
    ])
    dumux_error = dumux_T - analytical_at_dumux

# ── Plot ───────────────────────────────────────────────────────────────────────
fig, ax1 = plt.subplots(figsize=(10, 8))

ax1.plot(list(time_range), data_A.iloc[-1, :], "k.", markersize=10, markerfacecolor="none",
         label=f"Ramey A — OGS, different r in U and X")
ax1.plot(list(time_range), data_B.iloc[-1, :], "b.", markersize=10, markerfacecolor="none",
         label=f"Ramey B — same r in U and X")

# skip t=0
mask_ogs = ogs_time_s > 0
ax1.plot(ogs_time_s[mask_ogs], ogs_T[mask_ogs], "red", label="OGS")

if dumux_time_d is not None:
    mask_d = dumux_time_s > 0
    ax1.plot(dumux_time_s[mask_d], dumux_T[mask_d], "green", label="DuMux")

ax1.set_ylabel("Outlet temperature (°C)", fontsize=20)
ax1.set_xlabel("Time (d)", fontsize=20)
ax1.set_xticks(np.arange(0, 432000, 86400), np.arange(0, 5))
ax1.spines["left"].set_linewidth(2)
ax1.spines["bottom"].set_linewidth(2)
ax1.tick_params(axis="both", labelsize=16)

ax2 = ax1.twinx()
ax2.plot(ogs_time_s[mask_ogs], ogs_error_A[mask_ogs], color="red", linestyle="--",
         label="OGS error vs Ramey A")
ax2.plot(ogs_time_s[mask_ogs], ogs_error_B[mask_ogs], color="red", linestyle="-.",
         label="OGS error vs Ramey B")
if dumux_error is not None:
    ax2.plot(dumux_time_s[mask_d], dumux_error[mask_d], color="green", linestyle="--",
             label="DuMux error vs Ramey B")
ax2.set_ylabel("Absolute error (°C)", fontsize=20)
ax2.spines["right"].set_linewidth(2)
ax2.tick_params(axis="y", labelsize=16)
ax2.grid(axis="y", color="gray", linewidth=0.5)

print("Outlet temperature at t=5 days:")
print(f"Analytical B: {data_B.iloc[-1, -1]:.4f}")
if dumux_T is not None:
    print(f"DuMux:        {dumux_T[-1]:.4f}")
    print(f"Relative error: {(dumux_T[-1] - data_B.iloc[-1, -1]) / (data_B.iloc[-1, -1] - 20.0) * 100:.2f} %")

plt.figlegend(fontsize=14, loc="center right", bbox_to_anchor=(0.88, 0.5), frameon=False)
fig.tight_layout()
fig.savefig("wellbore_outlet_temperature.png", dpi=150)
print("Saved: wellbore_outlet_temperature.png")

# ── Second plot: temperature distribution along the borehole at final time ─────
# (mirrors OGS benchmark Fig. 3: distributed fluid temperature + absolute error
#  at the final timestep, t = 5 days)
f_t_final = time_function(dimensionless_time(lambda_re, t, rho_re, c_p_re, r_b))
X_final_A = coefficient_x(q, rho_f, c_p_f, lambda_re, r_pi, U_A, f_t_final)
X_final_B = coefficient_x(q, rho_f, c_p_f, lambda_re, r_pi, U_B, f_t_final)

ogs_error_A_profile = ogs_T_profile - outlet_temp(T_d, T_i, ogs_x_profile, X_final_A)
ogs_error_B_profile = ogs_T_profile - outlet_temp(T_d, T_i, ogs_x_profile, X_final_B)

fig2, bx1 = plt.subplots(figsize=(10, 8))
bx1.plot(Z, data_A.iloc[:, -1], "k.", markersize=10, markerfacecolor="none",
         label="Ramey A — OGS, different r in U and X")
bx1.plot(Z, data_B.iloc[:, -1], "b.", markersize=10, markerfacecolor="none",
         label="Ramey B — same r in U and X")
bx1.plot(ogs_x_profile, ogs_T_profile, "red", label="OGS")
if dumux_x is not None:
    bx1.plot(dumux_x, dumux_T_profile, "green", label="DuMux")
bx1.set_xlabel("Distance along borehole (m)", fontsize=20)
bx1.set_ylabel("Fluid temperature (°C)", fontsize=20)
bx1.set_xlim(0, 30)
bx1.spines["left"].set_linewidth(2)
bx1.spines["bottom"].set_linewidth(2)
bx1.tick_params(axis="both", labelsize=16)

bx2 = bx1.twinx()
bx2.plot(ogs_x_profile, ogs_error_A_profile, color="red", linestyle="--",
         label="OGS error vs Ramey A")
bx2.plot(ogs_x_profile, ogs_error_B_profile, color="red", linestyle="-.",
         label="OGS error vs Ramey B")
if dumux_x is not None:
    dumux_error_profile = dumux_T_profile - outlet_temp(T_d, T_i, dumux_x, X_final_B)
    bx2.plot(dumux_x, dumux_error_profile, color="green", linestyle="--",
             label="DuMux error vs Ramey B")
bx2.set_ylabel("Absolute error (°C)", fontsize=20)
bx2.spines["right"].set_linewidth(2)
bx2.tick_params(axis="y", labelsize=16)
bx2.grid(axis="y", color="gray", linewidth=0.5)

fig2.legend(fontsize=14, loc="center right", bbox_to_anchor=(0.88, 0.5), frameon=False)
fig2.tight_layout()
fig2.savefig("wellbore_temperature_distribution.png", dpi=150)
print("Saved: wellbore_temperature_distribution.png")

print("Temperature distribution at t=5 days:")
print(f"  Outlet (z=30 m) Ramey B:      {outlet_temp(T_d, T_i, 30.0, X_final_B):.4f}")
print(f"  Outlet (z=30 m) OGS:          {ogs_T_profile[-1]:.4f}")
print(f"  OGS max |error| vs Ramey B:   {np.abs(ogs_error_B_profile).max():.4f}")
if dumux_x is not None:
    print(f"  Outlet (z=30 m) DuMux:        {dumux_T_profile[-1]:.4f}")
    print(f"  DuMux max |error| vs Ramey B: {np.abs(dumux_error_profile).max():.4f}")

plt.show()
