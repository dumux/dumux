import numpy as np
import matplotlib.pyplot as plt

defaultPoreRadius = 2e-4

# default is water air system, 0.2mm pore radius, and cubic pores
def Pc(sw, sigma=0.0725, poreRadius=defaultPoreRadius, expFactor=6.83):
    Pc = 2.0*sigma / (poreRadius*(1.0 - np.exp(-expFactor*sw)))
    return Pc

def Sw(pc, sigma=0.0725, poreRadius=defaultPoreRadius, expFactor=6.83):
    Sw =-1.0/expFactor* np.log(1.0 - 2.0*sigma/(poreRadius*pc))
    return Sw

def dPc_dSw(sw, sigma=0.0725, poreRadius=defaultPoreRadius, expFactor=6.83):
    e = np.exp(expFactor*sw)
    dPc_dSw = -(2.0*expFactor*sigma*e) / (poreRadius*(1.0-e)*(1.0-e))
    return dPc_dSw


def PcEntry(theta=0, sigma=0.0725, radius=defaultPoreRadius*0.5):
    cosTheta = np.cos(theta)
    sinTheta = np.sin(theta)
    shapeFactor = 1.0/16.0
    factorD = np.pi - 3.0*theta + 3*sinTheta*cosTheta - (cosTheta*cosTheta) / (4.0*shapeFactor)
    factorF = (1 + np.sqrt(1 + 4*shapeFactor*factorD / (cosTheta*cosTheta))) / (1.0 + 2.0*np.sqrt(np.pi*shapeFactor))
    PcEntry = sigma / radius * cosTheta * (1 + 2*np.sqrt(np.pi*shapeFactor)) * factorF
    return PcEntry

Sw_array = np.linspace(0.01, 1, 1000)
Pc_array = Pc(Sw_array)
dPc_dSw_array = dPc_dSw(Sw_array)
PcValue = PcEntry()
SwValue = Sw(PcValue)
dPc_dSwValue = dPc_dSw(SwValue)
plt.plot(Sw_array, dPc_dSw_array)
plt.axhline(y = dPc_dSwValue, color = 'r', linestyle = '-')
plt.xlabel("Sw [-]")
plt.ylabel("dPc_dSw [Pa]")
plt.show()
