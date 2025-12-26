import os
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

FLUID = "CO2"

# Constantes críticas do próprio CoolProp (pra não hardcodar)
Tc = PropsSI("Tcrit", FLUID)       # K
Pc = PropsSI("pcrit", FLUID)       # Pa

# Temperatura "perto do crítico" mas acima para ficar monofásico
T = Tc + 1.0  # K

OUTDIR = "out"
os.makedirs(OUTDIR, exist_ok=True)

def safe_props(prop, in1, v1, in2, v2, fluid=FLUID):
    try:
        return PropsSI(prop, in1, v1, in2, v2, fluid)
    except Exception:
        return np.nan

# -----------------------
# 1) Z vs P (T fixo)
# -----------------------
# Faixa de pressão em torno do crítico: 0.2 Pc a 3 Pc
P_vals = np.linspace(0.2 * Pc, 3.0 * Pc, 400)
Z_vs_P = np.array([safe_props("Z", "T", T, "P", P, FLUID) for P in P_vals])

plt.figure()
plt.plot(P_vals / 1e6, Z_vs_P)
plt.xlabel("Pressão P (MPa)")
plt.ylabel("Fator de compressibilidade Z (-)")
plt.title(f"CO₂: Z vs P (T = Tc + 1 K = {T:.3f} K)")
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "co2_Z_vs_P_near_critical.png"), dpi=200)

# -----------------------
# 2) Z vs rho_molar (T fixo) e P vs rho_molar
# -----------------------
# Faixa típica de densidade molar perto do crítico do CO2:
# crítico ~ 10.6 kmol/m3 (≈ 10600 mol/m3). Vamos varrer bem ao redor.
rho_vals = np.linspace(500.0, 25000.0, 500)  # mol/m3

Z_vs_rho = np.array([safe_props("Z", "T", T, "Dmolar", rho, FLUID) for rho in rho_vals])
P_vs_rho = np.array([safe_props("P", "T", T, "Dmolar", rho, FLUID) for rho in rho_vals])

plt.figure()
plt.plot(rho_vals, Z_vs_rho)
plt.xlabel("Densidade molar ρ (mol/m³)")
plt.ylabel("Fator de compressibilidade Z (-)")
plt.title(f"CO₂: Z vs ρ (T = Tc + 1 K = {T:.3f} K)")
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "co2_Z_vs_rho_near_critical.png"), dpi=200)

plt.figure()
plt.plot(rho_vals, P_vs_rho / 1e6)
plt.xlabel("Densidade molar ρ (mol/m³)")
plt.ylabel("Pressão P (MPa)")
plt.title(f"CO₂: P vs ρ (T = Tc + 1 K = {T:.3f} K)")
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "co2_P_vs_rho_near_critical.png"), dpi=200)

print("Gerado:")
print(" - out/co2_Z_vs_P_near_critical.png")
print(" - out/co2_Z_vs_rho_near_critical.png")
print(" - out/co2_P_vs_rho_near_critical.png")
print(f"Tc = {Tc:.6f} K, Pc = {Pc/1e6:.6f} MPa, T usado = {T:.6f} K")
