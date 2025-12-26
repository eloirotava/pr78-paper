import time
import numpy as np
from CoolProp.CoolProp import PropsSI

FLUID = "CO2"
Tc = PropsSI("Tcrit", FLUID)

# use exatamente os mesmos ranges do Rust
t_min = Tc + 0.1
t_max = Tc + 10.0
n_t = 81

rho_min = 2000.0
rho_max = 16000.0
n_rho = 120

Ts = np.linspace(t_min, t_max, n_t)
rhos = np.linspace(rho_min, rho_max, n_rho)

# warm-up (importante)
for _ in range(500):
    PropsSI("Z", "T", float(Ts[0]), "Dmolar", float(rhos[0]), FLUID)

t0 = time.perf_counter()
acc = 0.0
count = 0

for T in Ts:
    for rho in rhos:
        z = PropsSI("Z", "T", float(T), "Dmolar", float(rho), FLUID)
        acc += z
        count += 1

dt = time.perf_counter() - t0
print(f"CoolProp: pontos={count}, tempo={dt:.6f} s, pontos/s={count/dt:.3e}, checksum={acc:.6f}")
