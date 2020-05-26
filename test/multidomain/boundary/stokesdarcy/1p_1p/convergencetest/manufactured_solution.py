#!/usr/bin/env python3

"""
This script determines the source terms needed for analytical solutions for coupled Stokes/Darcy problems.
Given an analytical solution, it evaluates the momentum and mass balance equations and outputs source terms
such that the balance equations are fulfilled. It furthermore calculates the Darcy velocity from a given pressure solution.
"""

import sympy as sp
import numpy as np
from sympy import *

# case = "Rybak"
# case = "Shiue_1"
case = "Shiue_2"

# divergence of a matrix
def div(M):
    return np.array([sp.diff(M[0,0],x) + sp.diff(M[0,1],y) , sp.diff(M[1,0],x) + sp.diff(M[1,1],y) ])

# gradient of a vector
def grad(v):
    return np.array([[sp.diff(v[0],x),  sp.diff(v[0],y)], [sp.diff(v[1],x),  sp.diff(v[1],y)]])

# coordinates
x, y = sp.symbols("x y")

# analytical solutions for v and p in free-flow domain (see reference given in problems)
def analyticalSolutionStokes(case):
    if case == "Rybak":
        vFF = np.array([-sp.cos(sp.pi*x)*sp.sin(sp.pi*y), sp.sin(sp.pi*x)*sp.cos(sp.pi*y)])
        pFF = 0.5*y*sp.sin(sp.pi*x)
    elif case == "Shiue_1":
        vFF = np.array([-1.0/sp.pi * sp.exp(y) * sp.sin(sp.pi*x), (sp.exp(y) - sp.exp(1)) * sp.cos(sp.pi*x)])
        pFF = 2.0*sp.exp(y) * sp.cos(sp.pi*x)
    else: # Shiue_2
        vFF = np.array([(y-1.0)*(y-1.0) + x*(y-1.0) + 3.0*x - 1.0, x*(x-1.0) - 0.5*(y-1.0)*(y-1.0) - 3.0*y + 1.0])
        pFF = 2.0*x + y - 1.0
    return [vFF, pFF]

vFF, pFF = analyticalSolutionStokes(case)

# individual terms of the Navier-Stokes eq.
vvT = np.outer(vFF, vFF)
gradV = grad(vFF)
gradVT = grad(vFF).T
pI = np.array([[pFF,  0], [0,  pFF]])

# complete momentum flux and its divergence
#momentumFlux = vvT - (gradV + gradVT) +pI
momentumFlux = - (gradV + gradVT) +pI # only Stokes
divMomentumFlux = div(momentumFlux)

print("Source terms for case", case)

# solution for source term
print(" \nStokes:")
print("Source term mass:", sp.diff(vFF[0],x) + sp.diff(vFF[1], y))
print("Source term momentum x:", divMomentumFlux[0])
print("Source term momentum y:", divMomentumFlux[1])

# analytical solutions for p in Darcy domain (see reference given in problems)
def analyticalSolutionDarcy(case):
    if case == "Rybak":
        pD = 0.5*y*y*sp.sin(sp.pi*x)
    elif case == "Shiue_1":
        pD = (sp.exp(y) - y*sp.exp(1)) * sp.cos(sp.pi*x)
    else: # Shiue_2
        pD = x*(1.0-x)*(y-1.0) + pow(y-1.0, 3)/3.0 + 2.0*x + 2.0*y + 4.0
    return pD

pD = analyticalSolutionDarcy(case)

gradPdK = np.array([sp.diff(pD,x), sp.diff(pD,y)])
vD = -gradPdK
divVd = sp.diff(vD[0],x) + sp.diff(vD[1],y)

print("\nDarcy:")
print("Source term mass:", simplify(divVd))
print("v x:", simplify(vD[0]))
print("v_y", simplify(vD[1]))
