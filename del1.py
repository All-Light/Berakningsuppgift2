import numpy as np
import matplotlib.pyplot as plt
import pyfluids as pf
from scipy.optimize import brentq
WATER = pf.Fluid(pf.FluidsList.Water)

# Ett system utan matarvattenförvärmning
# 4 olika tillstånd
#   1 - Precis efter kondensator - Mättad vätska. P1 = P4, S1 = S2
#   2 - Precis efter pump - Trycksatt vätska. P2 = P3, S1 = S2
#   3 - Precis efter ångpanna. # T3, P2 = P3, S3 = S4
#   4 - Precis efter turbin. # P1 = P4, S3 = S4


P2 = P3 = 13*1e6 # Pa
#print(P3)
T3 = 600 # C 
WATER.update(pf.Input.pressure(P3), pf.Input.temperature(T3))

H3 = WATER.enthalpy # J/(kg)
S3 = WATER.entropy # J/(kg*K)
#print("H3: ",H3)
#print("S3: ",S3)

def f(P):
    WATER.update(pf.Input.pressure(P), pf.Input.quality(0))
    s_f = WATER.entropy
    WATER.update(pf.Input.pressure(P), pf.Input.quality(100))
    s_g = WATER.entropy
    return s_f + 0.85 * (s_g - s_f) - S3  # s3 from state 3

def show_graph():
    x = np.linspace(1000, 1e7, 1001)
    y = np.zeros_like(x)
    for index, val in enumerate(x):
        y[index] = f(val)
    plt.plot(x,y)
    plt.show()

# Vi löser ut trycket eftersom vi inte kan använda entropi och ångkvalité för att uppdatera vattnet
P4 = brentq(f,1000, 1e7)
#print("P4: ",P4)

# Tillstånd 4
# Ångkvalitén minst 85% i tillstånd 4. Entropin är samma som i tillstånd 3
# Trycket i tillstånd 4 är samma som tillstånd 1
WATER.update(pf.Input.pressure(P4),  pf.Input.quality(85))
H4 = WATER.enthalpy
T4 = WATER.temperature
#print("H4: ", H4)
#print("T4: ", T4)

# Tillstånd 1 - Efter kondensorn
WATER.update(pf.Input.pressure(P4), pf.Input.quality(0))
H1 = WATER.enthalpy
S2 = S1 = WATER.entropy
#print("H1: ", H1)
#print("S1: ", S1)

# Tillstånd 2 - Efter pumpen
WATER.update(pf.Input.pressure(P2), pf.Input.entropy(S2))
H2 = WATER.enthalpy

#print("H2: ", H2)

W_pump = (H2 - H1)/1000
W_turb =(H3 - H4)/1000
Q_in = (H3 - H2)/1000
print("Pump energy usage (KJ/kg): ", W_pump)
print("Turbine energy generation (KJ/kg): ", W_turb)
print("Heater energy input (KJ/kg): ", Q_in)

W_net = W_turb - W_pump 
verkningsgrad = W_net/Q_in
print(f"Verkningsgrad (%): {verkningsgrad*100}")
print(f"Kondensortryck {P4/1000} kPa")
