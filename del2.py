import numpy as np
import matplotlib.pyplot as plt
import pyfluids as pf
from scipy.optimize import brentq

WATER = pf.Fluid(pf.FluidsList.Water)

# Ett system MED matarvattenförvärmning
# 7 olika tillstånd
#   1 - Precis efter kondensator - Mättad vätska. P1 = P7, S1 = S7
#   2 - Precis efter pump 1 - Trycksatt vätska. P2 = P3 = P6, S1 = S2
#   3 - Precis efter matarvattenförvärmning. # P2 = P3 = P6
#   4 - Precis efter pump 2. # P4 = P5, S3 = S4
#   5 - Precis efter ångpanna. # P4 = P5, T5
#   6 - Precis efter turbin fraktion y. # S5 = S6 = S7
#   7 - Precis efter turbin fraktion 1-y. # S5 = S6 = S7

P5 = 13*1e6 # Pa
T5 = 600 # C 
WATER.update(pf.Input.pressure(P5), pf.Input.temperature(T5))

H5 = WATER.enthalpy # J/(kg)
S5 = WATER.entropy # J/(kg*K)
#print("H5: ", H5)
#print("S5: ", S5)

# Givet från förra uppgiften:
#   Kondensatortryck:
P1 = P7 = 28030.389667251104 # Pa

WATER.update(pf.Input.pressure(P7),  pf.Input.quality(85))
H7 = WATER.enthalpy
S7 = WATER.entropy
T7 = WATER.temperature
#print("H7: ", H7)
#print("S7: ", S7)
#print("T7: ", T7)

# Tillstånd 1 har samma tryck som tillstånd 7 och ångkvalité 0% (mättad vätska)
WATER.update(pf.Input.pressure(P7), pf.Input.quality(0))
H1 = WATER.enthalpy
S1 = WATER.entropy
T1 = WATER.temperature
#print("H1: ", H1)
#print("S1: ", S1)
#print("T1: ", T1)

# SÖKES: tryck i tillstånd 2, 3, 6
#   Vi vet att tillstånd 3 är mättad vätska (ångkvalité 0%) och samma tryck som tillstånd 6.
#   Vi vet att tillstånd 6 är överhettad ånga
def f(P, y): 
    # Beräkna trycket P som krävs vid fraktion y från turbinen (Newtons metod)
    #   för att energiprincipen ska vara uppfylld 
    
    # isentrop pump mellan 1-2 ger H2
    WATER.update(pf.Input.pressure(P), pf.Input.entropy(S1)) 
    H2 = WATER.enthalpy

    # givet mättad vätska vid samma tryck som 2 och 6 ger H3
    WATER.update(pf.Input.pressure(P), pf.Input.quality(0))
    H3 = WATER.enthalpy

    # isentrop turbin ger att entropin är samma i 5-6-7 och ger då H6 för trycket
    WATER.update(pf.Input.pressure(P), pf.Input.entropy(S7))
    H6 = WATER.enthalpy

    # Energiprincipen => delta E_sys = 0 sök P givet y.
    return ((1-y)*H2 + y*H6) - H3


# SÖKES: tryck i tillstånd 2, 3, 6
#   Vi vet att tillstånd 3 är mättad vätska (ångkvalité 0%) och samma tryck som tillstånd 6.
#   Vi vet att tillstånd 6 är överhettad ånga

'''
Visualisera tryck-energy
def show_graph():
    percentage = np.linspace(0, 1, 101)
    x = np.linspace(1000, 1e7, 101)
    y = np.zeros(shape=[len(percentage), len(x)])
    for index1, percent in enumerate(percentage):
        #print("new percentage: ", percent)
        for index2, pressure in enumerate(x):
            y[index1, index2] = f(pressure, percent)
    #print(y)
    plt.plot(x.T, y.T)

    plt.xlabel("Pressure")
    plt.ylabel("Energy ")
    plt.show()
show_graph()
'''

percentages = np.linspace(0, 1, 1001) # olika fraktioner y

# Skapa matris av trycket P vid de olika fraktionerna y. 
#   sol har struktur [[y_0, P_0], [y_1, P_1]] då vi behöver
#   vet vilken fraktion som motsvarar vilket tryck.

sol = np.zeros([len(percentages),2])
for index, y in enumerate(percentages):
    try:
        min_pressure = 1000 #Pa
        max_pressure = P5 # impossible to have higher pressure than P5
        sol[index,:] = (y, brentq(f, min_pressure, max_pressure, args=y))
    except:
        pass

# bevara endast där trycket är större än 0 (dvs vi har en lösning)
sol = sol[sol[:, 1] > 0]

# då vi har en matris av tryck P får vi även matriser för tillstånd 2,3 och 6
S3 = np.zeros_like(sol[:,1])
H3 = np.zeros_like(sol[:,1])

S2 = np.zeros_like(sol[:,1])
H2 = np.zeros_like(sol[:,1])

S6 = np.zeros_like(sol[:,1])
H6 = np.zeros_like(sol[:,1])

for index, row in enumerate(sol):
    P = row[1]
    if P == 0:
        continue
    # 3 - mättad ånga vid tryck P
    WATER.update(pf.Input.pressure(P), pf.Input.quality(0))
    S3[index] = WATER.entropy
    H3[index] = WATER.enthalpy
    # 2 - samma entropi som 1 (isentrop pump) vid tryck P
    WATER.update(pf.Input.pressure(P), pf.Input.entropy(S1))
    S2[index] = WATER.entropy
    H2[index] = WATER.enthalpy
    # 6 - samma entropi som 5 (isentrop turbin) vid tryck P
    WATER.update(pf.Input.pressure(P), pf.Input.entropy(S5))
    S6[index] = WATER.entropy
    H6[index] = WATER.enthalpy

# Tillstånd 4 har samma entropi som tillstånd 3 och samma tryck som tillstånd 5
#   pga isobar ångpanna
H4 = np.empty_like(S3)
for index, entropy in enumerate(S3):
    WATER.update(pf.Input.entropy(entropy), pf.Input.pressure(P5))
    H4[index] = WATER.enthalpy 


# vi får matriser för pumparbetet och värmetillförsel
W_pump1 = H2 - H1 # arbete i pump 1 är entalpiskillnaden
W_pump2 = H4 - H3 # arbete i pump 2 är entalpiskillnaden
Q_in = H5 - H4 # värmetillförsel i ångpannan är entalpiskillnaden

W_turb = np.empty_like(sol)
for index, row in enumerate(sol):
    y = row[0]
    # systemets arbete på turbinen på formen (fraktion, energiutvinning)
    W_turb[index] = (y, H5 - H6[index] +(1-y)*(H6[index] - H7))

W_net = W_turb[:,1] - W_pump1 - W_pump2
verkningsgrad = W_net / Q_in
max_index = np.argmax(verkningsgrad) # maximala verkningsgradens index

#print("W_net: ", W_net)
#print("Q_in: ", Q_in)
#print("verkningsgrad: ",verkningsgrad)
print(f"Maximala verkningsgraden fås vid y={sol[max_index,0]*100}% då verkningsgraden är {np.max(verkningsgrad)*100}%. ")
print(f"Trycket uppgår till ca {sol[max_index,1]/1e6} MPa")