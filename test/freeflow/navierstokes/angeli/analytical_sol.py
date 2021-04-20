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

# coordinates and time
x, y, t = sp.symbols("x y t")

kinematicViscosity = 0.1

# analytical solutions for v and p in free-flow domain (see reference given in problems)
def analyticalSolutionStokes():
    psi = sp.exp(-5*kinematicViscosity*sp.pi**2*t)*sp.cos(sp.pi*x)*sp.cos(2*sp.pi*y)
    ux = sp.diff(psi,y)
    uy = -sp.diff(psi,x)

    vFF = np.array([- 2.0 * sp.pi * sp.exp(- 5.0 * kinematicViscosity * sp.pi * sp.pi * t) * sp.cos(sp.pi * x) * sp.sin(2.0 * sp.pi * y), sp.pi * sp.exp(- 5.0 * kinematicViscosity * sp.pi*sp.pi * t) * sp.sin(sp.pi * x) * sp.cos(2.0 * sp.pi * y)])

    assert(vFF[0]-ux == 0)
    assert(vFF[1]-uy == 0)

    pFF = - 0.25 * sp.exp(-10.0 * kinematicViscosity * sp.pi**2 * t) * sp.pi**2 * (4.0 * sp.cos(2.0 * sp.pi * x) + sp.cos(4.0 * sp.pi * y))
    return [vFF, pFF]

vFF, pFF = analyticalSolutionStokes()

# individual terms of the Navier-Stokes eq.
#vvT = np.outer(vFF, vFF)
vvT = np.tensordot(vFF, vFF, axes=0)
gradV = grad(vFF)
gradVT = grad(vFF).T
pI = np.array([[pFF,  0], [0,  pFF]])
storageTerm = np.array([sp.diff(vFF[0], t), sp.diff(vFF[1], t)])

# complete momentum flux and its divergence
momentumFlux = vvT - kinematicViscosity*(gradV + gradVT) +pI
divMomentumFlux = div(momentumFlux)

# solution for source term
print("Source term mass:", sp.diff(vFF[0],x) + sp.diff(vFF[1], y))
print("Source term momentum x:", sp.simplify(divMomentumFlux[0] + storageTerm[0]))
print("Source term momentum y:", sp.simplify(divMomentumFlux[1] + storageTerm[1]))
