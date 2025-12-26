import os
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

FLUID = "CO2"
OUTDIR = "out"
os.makedirs(OUTDIR, exist_ok=True)

Tc = PropsSI("Tcrit", FLUID)
Pc = PropsSI("pcrit", FLUID)

# Temperaturas: acima do crítico (principal) + uma abaixo (ilustrativa)
T_list_above = [Tc + dT for dT in [0.1, 1.0, 2.0, 5.0, 10.0]]
T_below = Tc - 1.0

def safe_props(prop, in1, v1, in2, v2):
    try:
        return PropsSI(prop, in1, v1, in2, v2, FLUID)
    except Exception:
        return np.nan

def compute_vs_P(T, P_vals):
    Z = np.array([safe_props("Z", "T", T, "P", P) for P in P_vals], dtype=float)
    rho = np.array([safe_props("Dmolar", "T", T, "P", P) for P in P_vals], dtype=float)  # mol/m3
    ok = np.isfinite(Z) & np.isfinite(rho)
    return P_vals[ok], Z[ok], rho[ok]

def compute_vs_rho(T, rho_vals):
    Z = np.array([safe_props("Z", "T", T, "Dmolar", r) for r in rho_vals], dtype=float)
    P = np.array([safe_props("P", "T", T, "Dmolar", r) for r in rho_vals], dtype=float)
    ok = np.isfinite(Z) & np.isfinite(P)
    return rho_vals[ok], Z[ok], P[ok]

def add_top_axis_rho_from_P(ax, P_curve, rho_curve, nticks=6):
    """
    Add a top x-axis where tick labels show rho corresponding to selected P ticks.
    Assumes monotonic rho(P) over the plotted range (true above critical for chosen ranges).
    """
    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())

    # Choose P tick positions from bottom axis ticks (or evenly spaced)
    P_ticks = np.linspace(P_curve.min(), P_curve.max(), nticks)
    rho_ticks = np.interp(P_ticks, P_curve, rho_curve)

    ax_top.set_xticks(P_ticks / 1e6)  # but top axis uses same transform as bottom; we'll relabel
    # Trick: set ticks in MPa coordinates on bottom, then relabel with rho.
    # Better: keep same scale as bottom (MPa), and label rho.
    ax_top.set_xlim(ax.get_xlim())
    ax_top.set_xticks(np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], nticks))

    # Map those MPa positions back to Pa to interpolate rho
    P_ticks_Pa = np.linspace(P_curve.min(), P_curve.max(), nticks)
    rho_ticks = np.interp(P_ticks_Pa, P_curve, rho_curve)

    ax_top.set_xticklabels([f"{rt:,.0f}" for rt in rho_ticks])
    ax_top.set_xlabel("Densidade molar ρ (mol/m³) correspondente aos pontos da curva")
    return ax_top

def add_top_axis_P_from_rho(ax, rho_curve, P_curve, nticks=6):
    """
    Add a top x-axis where tick labels show P corresponding to selected rho ticks.
    Assumes monotonic P(rho) over the plotted range.
    """
    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())

    rho_ticks = np.linspace(rho_curve.min(), rho_curve.max(), nticks)
    P_ticks = np.interp(rho_ticks, rho_curve, P_curve)

    ax_top.set_xticks(np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], nticks))
    ax_top.set_xticklabels([f"{pt/1e6:,.2f}" for pt in P_ticks])
    ax_top.set_xlabel("Pressão P (MPa) correspondente aos pontos da curva")
    return ax_top

# -----------------------------
# FIG 1: Z vs P (bottom), with rho on top (same curve)
# -----------------------------
P_vals = np.linspace(0.2 * Pc, 3.0 * Pc, 600)

plt.figure(figsize=(9, 5))
ax = plt.gca()

# Plot above-critical temperatures
for T in T_list_above:
    P_ok, Z_ok, rho_ok = compute_vs_P(T, P_vals)
    ax.plot(P_ok / 1e6, Z_ok, label=f"T = Tc + {T - Tc:.1f} K")

# Plot below-critical (illustrative)
P_ok_b, Z_ok_b, rho_ok_b = compute_vs_P(T_below, P_vals)
ax.plot(P_ok_b / 1e6, Z_ok_b, linestyle="--", label=f"T = Tc - {Tc - T_below:.1f} K (ilustrativo)")

ax.set_xlabel("Pressão P (MPa)")
ax.set_ylabel("Fator de compressibilidade Z (-)")
ax.set_title("CO₂ – Z vs P (eixo superior: ρ correspondente ao longo da curva)")
ax.grid(True)
ax.legend(loc="best", fontsize=8)

# Top axis uses rho(P) for ONE representative curve (pick Tc+1 K for clarity)
T_ref = Tc + 1.0
P_ref, Z_ref, rho_ref = compute_vs_P(T_ref, P_vals)
# Set top axis labels based on this reference mapping
ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
# choose tick positions in MPa and map to Pa
tick_pos = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 6)  # MPa positions
P_ticks_Pa = tick_pos * 1e6
rho_ticks = np.interp(P_ticks_Pa, P_ref, rho_ref)
ax_top.set_xticks(tick_pos)
ax_top.set_xticklabels([f"{rt:,.0f}" for rt in rho_ticks])
ax_top.set_xlabel("Densidade molar ρ (mol/m³) (mapeada usando T = Tc + 1 K)")

plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "co2_Z_vs_P_with_top_rho.png"), dpi=200)

# -----------------------------
# FIG 2: Z vs rho (bottom), with P on top (same curve)
# -----------------------------
rho_vals = np.linspace(500.0, 25000.0, 700)

plt.figure(figsize=(9, 5))
ax = plt.gca()

for T in T_list_above:
    rho_ok, Z_ok, P_ok = compute_vs_rho(T, rho_vals)
    ax.plot(rho_ok, Z_ok, label=f"T = Tc + {T - Tc:.1f} K")

rho_ok_b, Z_ok_b, P_ok_b = compute_vs_rho(T_below, rho_vals)
ax.plot(rho_ok_b, Z_ok_b, linestyle="--", label=f"T = Tc - {Tc - T_below:.1f} K (ilustrativo)")

ax.set_xlabel("Densidade molar ρ (mol/m³)")
ax.set_ylabel("Fator de compressibilidade Z (-)")
ax.set_title("CO₂ – Z vs ρ (eixo superior: P correspondente ao longo da curva)")
ax.grid(True)
ax.legend(loc="best", fontsize=8)

# Top axis labels based on reference mapping at Tc+1 K
rho_ref, Z_ref, P_ref = compute_vs_rho(T_ref, rho_vals)
ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
tick_pos = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 6)  # rho positions
P_ticks = np.interp(tick_pos, rho_ref, P_ref)
ax_top.set_xticks(tick_pos)
ax_top.set_xticklabels([f"{pt/1e6:,.2f}" for pt in P_ticks])
ax_top.set_xlabel("Pressão P (MPa) (mapeada usando T = Tc + 1 K)")

plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "co2_Z_vs_rho_with_top_P.png"), dpi=200)

print("Gerado:")
print(" - out/co2_Z_vs_P_with_top_rho.png")
print(" - out/co2_Z_vs_rho_with_top_P.png")
print(f"Tc = {Tc:.6f} K, Pc = {Pc/1e6:.6f} MPa")
