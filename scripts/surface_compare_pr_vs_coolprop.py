import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

FLUID = "CO2"
INFILE = "out/pr78_co2_grid.csv"
OUTDIR = "out"
os.makedirs(OUTDIR, exist_ok=True)

df = pd.read_csv(INFILE)

# Remove linhas inválidas do PR (NaN)
df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["Z_PR", "P_Pa"])

# Calcula referência CoolProp nos mesmos pontos
T = df["T_K"].to_numpy()
rho = df["rho_mol_m3"].to_numpy()

Z_cp = np.array([PropsSI("Z", "T", t, "Dmolar", r, FLUID) for t, r in zip(T, rho)], dtype=float)

df["Z_CP"] = Z_cp
df["dZ"] = df["Z_PR"] - df["Z_CP"]
df["rel_dZ"] = df["dZ"] / df["Z_CP"]

# Transformar em grade (T, rho) -> valores
Ts = np.sort(df["T_K"].unique())
rhos = np.sort(df["rho_mol_m3"].unique())

def gridify(col):
    pivot = df.pivot(index="T_K", columns="rho_mol_m3", values=col)
    pivot = pivot.reindex(index=Ts, columns=rhos)
    return pivot.to_numpy()

ZPR = gridify("Z_PR")
ZCP = gridify("Z_CP")
dZ  = gridify("dZ")
rel = gridify("rel_dZ")

TT, RR = np.meshgrid(rhos, Ts)  # RR no x? aqui x=rho, y=T (pela pivot)

def surface_plot(Z, title, fname):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(RR, TT, Z, linewidth=0, antialiased=True)
    ax.set_ylabel("ρ (mol/m³)")
    ax.set_xlabel("T (K)")
    ax.set_zlabel("Z (-)")
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, fname), dpi=200)

surface_plot(ZPR, "CO₂ – Superfície Z(ρ,T) – PR78", "surf_Z_PR78.png")
surface_plot(ZCP, "CO₂ – Superfície Z(ρ,T) – CoolProp", "surf_Z_CoolProp.png")

# Erro (ΔZ) e erro relativo
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.plot_surface(RR, TT, dZ, linewidth=0, antialiased=True)
ax.set_ylabel("ρ (mol/m³)")
ax.set_xlabel("T (K)")
ax.set_zlabel("ΔZ = Z_PR − Z_ref")
ax.set_title("CO₂ – Erro absoluto ΔZ(ρ,T) – PR78 vs CoolProp")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "surf_dZ.png"), dpi=200)

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.plot_surface(RR, TT, rel, linewidth=0, antialiased=True)
ax.set_ylabel("ρ (mol/m³)")
ax.set_xlabel("T (K)")
ax.set_zlabel("Erro relativo (ΔZ/Z_ref)")
ax.set_title("CO₂ – Erro relativo em Z(ρ,T) – PR78 vs CoolProp")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "surf_rel_dZ.png"), dpi=200)

print("Gerado:")
print(" - out/surf_Z_PR78.png")
print(" - out/surf_Z_CoolProp.png")
print(" - out/surf_dZ.png")
print(" - out/surf_rel_dZ.png")
