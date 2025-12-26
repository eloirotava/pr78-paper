import os
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

FLUID = "CO2"
OUTDIR = "out"
os.makedirs(OUTDIR, exist_ok=True)

Tc = PropsSI("Tcrit", FLUID)
Pc = PropsSI("pcrit", FLUID)

# Temperaturas desejadas (inclui abaixo e acima)
dT_list = [-5.0, -2.0, -1.0, -0.1, 0.1, 1.0, 2.0, 5.0, 10.0, 50]
T_list = [Tc + dT for dT in dT_list]

def safe_props(prop, in1, v1, in2, v2):
    try:
        return PropsSI(prop, in1, v1, in2, v2, FLUID)
    except Exception:
        return np.nan

def series_from_P(T, P_min, P_max, n=220):
    P = np.linspace(P_min, P_max, n)
    Z = np.array([safe_props("Z", "T", T, "P", p) for p in P], dtype=float)
    ok = np.isfinite(Z)
    P = P[ok]; Z = Z[ok]
    s = np.linspace(0.0, 1.0, len(Z))  # parâmetro normalizado
    return s, Z

def series_from_rho(T, rho_min, rho_max, n=220):
    rho = np.linspace(rho_min, rho_max, n)
    Z = np.array([safe_props("Z", "T", T, "Dmolar", r) for r in rho], dtype=float)
    ok = np.isfinite(Z)
    rho = rho[ok]; Z = Z[ok]
    s = np.linspace(0.0, 1.0, len(Z))  # parâmetro normalizado
    return s, Z

# Faixas (ajuste se quiser)
P_min, P_max = 0.2 * Pc, 3.0 * Pc
rho_min, rho_max = 500.0, 25000.0

plt.figure(figsize=(10, 6))
ax = plt.gca()

# Para não lotar a legenda com 18 entradas, vamos:
# - legendas por temperatura (linhas) e marcador indica tipo (P vs rho)
# - e incluir um "proxy" para os marcadores no final
for dT, T in zip(dT_list, T_list):
    # P-varredura (bolinha)
    sP, ZP = series_from_P(T, P_min, P_max)
    ax.plot(sP, ZP, marker="o", markersize=2.5, linewidth=0.8, label=f"T = Tc {dT:+g} K (varre P)")

    # rho-varredura (quadrado)
    sR, ZR = series_from_rho(T, rho_min, rho_max)
    ax.plot(sR, ZR, marker="s", markersize=2.5, linewidth=0.8, label=f"T = Tc {dT:+g} K (varre ρ)")

ax.set_xlabel("Parâmetro de varredura normalizado s (0→1)")
ax.set_ylabel("Fator de compressibilidade Z (-)")
ax.set_title("CO₂ – Comparação do comportamento numérico: varredura por P (●) vs varredura por ρ (■)")
ax.grid(True)

# Legenda: pode ficar grande; alternativa é salvar sem legenda e fazer no Word.
ax.legend(loc="best", fontsize=7, ncol=2)

plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "co2_overlay_Z_P_vs_rho_s.png"), dpi=200)

print("Gerado: out/co2_overlay_Z_P_vs_rho_s.png")
print(f"Tc={Tc:.6f} K, Pc={Pc/1e6:.3f} MPa")
