#!/usr/bin/env python3

import sympy as sp
import numpy as np
from sympy import *

# divergence of a matrix
def div(M):
    return np.array([sp.diff(M[0,0],x) + sp.diff(M[0,1],y) , sp.diff(M[1,0],x) + sp.diff(M[1,1],y) ])

# gradient of a vector
def grad(v):
    return np.array([[sp.diff(v[0],x),  sp.diff(v[0],y)], [sp.diff(v[1],x),  sp.diff(v[1],y)]])

# coordinates
x, y = sp.symbols("x y")

kinematicViscosity = 0.1
t = 1.0/(10*kinematicViscosity)

# analytical solutions for v and p in free-flow domain (see reference given in problems)
def analyticalSolutionStokes():

    vFF = np.array([- 2.0 * sp.pi * sp.exp(- 5.0 * kinematicViscosity * sp.pi * sp.pi * t) * sp.cos(sp.pi * x) * sp.sin(2.0 * sp.pi * y), sp.pi * sp.exp(- 5.0 * kinematicViscosity * sp.pi*sp.pi * t) * sp.sin(sp.pi * x) * sp.cos(2.0 * sp.pi * y)])
    pFF = - 0.25 * sp.exp(-10.0 * kinematicViscosity * sp.pi * sp.pi * t) * sp.pi * sp.pi * (4.0 * sp.cos(2.0 * sp.pi * x) + sp.cos(4.0 * sp.pi * y))
    return [vFF, pFF]

vFF, pFF = analyticalSolutionStokes()

# individual terms of the Navier-Stokes eq.
vvT = np.outer(vFF, vFF)
gradV = grad(vFF)
gradVT = grad(vFF).T
pI = np.array([[pFF,  0], [0,  pFF]])

# complete momentum flux and its divergence
momentumFlux = vvT - (gradV + gradVT) +pI
#momentumFlux = - (gradV + gradVT) +pI # only Stokes
divMomentumFlux = div(momentumFlux)

# solution for source term
print("Source term mass:", sp.diff(vFF[0],x) + sp.diff(vFF[1], y))
print("Source term momentum x:", divMomentumFlux[0])
print("Source term momentum y:", divMomentumFlux[1])
