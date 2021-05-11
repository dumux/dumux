#!/usr/bin/env python3

import sympy as sp
import numpy as np
from sympy import *

# divergence of a matrix
def div(M):
    return np.array([sp.diff(M[0,0],x) + sp.diff(M[0,1],y) + sp.diff(M[0,2],z),
                     sp.diff(M[1,0],x) + sp.diff(M[1,1],y) + sp.diff(M[1,2],z),
                     sp.diff(M[2,0],x) + sp.diff(M[2,1],y) + sp.diff(M[2,2],z)])

# gradient of a vector
def grad(v):
    return np.array([[sp.diff(v[0],x),  sp.diff(v[0],y), sp.diff(v[0],z)],
                     [sp.diff(v[1],x),  sp.diff(v[1],y), sp.diff(v[1],z)],
                     [sp.diff(v[2],x),  sp.diff(v[2],y), sp.diff(v[2],z)]])

# coordinates and time
x, y, z, t = sp.symbols("x y z t")

kinematicViscosity = 1.0

# analytical solutions for v and p in free-flow domain (see reference given in problems)
def analyticalSolutionStokes():
    a = sp.pi / 4.0
    d = sp.pi / 2.0
    nu = kinematicViscosity
    ax = a*x
    ay = a*y
    az = a*z
    dx = d*x
    dy = d*y
    dz = d*z
    f = sp.exp(-(d**2)*nu*t)

    eax = sp.exp(ax)
    eay = sp.exp(ay)
    eaz = sp.exp(az)
    ux = -a * (eax*sp.sin(ay + dz) + eaz*sp.cos(ax + dy)) * f
    uy = -a * (eay*sp.sin(az + dx) + eax*sp.cos(ay + dz)) * f
    uz = -a * (eaz*sp.sin(ax + dy) + eay*sp.cos(az + dx)) * f

    vFF = np.array([ux, uy, uz])

    # pFF = -0.5*a**2 *(sp.exp(2.0*ax) + sp.exp(2.0*ay) + sp.exp(2.0*az)
    #       + 2.0*sp.sin(ax + dy)*sp.cos(az + dx)
    #       + 2.0*sp.sin(ay + dz)*sp.cos(ax + dy)
    #       + 2.0*sp.sin(az + dx)*sp.cos(ay + dz)) * sp.exp(-2.0*(d**2)*nu*t)
    pFF = -0.5*a**2 *(sp.exp(2.0*ax) + sp.exp(2.0*ay) + sp.exp(2.0*az)
          + 2.0*sp.sin(ax + dy)*sp.cos(az + dx)*sp.exp(a*(y+z))
          + 2.0*sp.sin(ay + dz)*sp.cos(ax + dy)*sp.exp(a*(z+x))
          + 2.0*sp.sin(az + dx)*sp.cos(ay + dz)*sp.exp(a*(x+y))) * sp.exp(-2.0*(d**2)*nu*t)

    return [vFF, pFF]

vFF, pFF = analyticalSolutionStokes()

# individual terms of the Navier-Stokes eq.
vvT = np.outer(vFF, vFF)
# vvT = np.tensordot(vFF, vFF, axes=0)
gradV = grad(vFF)
gradVT = grad(vFF).T
pI = np.array([[pFF, 0, 0],
               [0,  pFF, 0],
               [0,  0, pFF]])

storageTerm = np.array([sp.diff(vFF[0], t),
                        sp.diff(vFF[1], t),
                        sp.diff(vFF[2], t)])

# complete momentum flux and its divergence
momentumFlux = vvT - kinematicViscosity*(gradV + gradVT) + pI
divMomentumFlux = div(momentumFlux)

divV = sp.simplify(sp.diff(vFF[0],x) + sp.diff(vFF[1], y) + sp.diff(vFF[2], z))

# solution for source term
print("Source term mass:", divV)
print("Source term momentum x:", sp.simplify(divMomentumFlux[0] + storageTerm[0]))
# print("Source term momentum y:", sp.simplify(divMomentumFlux[1] ))#+ storageTerm[1]))
#
print("unsteady terms balance viscous terms in the momentum equations")
print("storage - viscous terms 0", sp.simplify(storageTerm[0] - kinematicViscosity*div(gradV + gradVT)[0]))
print("storage - viscous terms 1", sp.simplify(storageTerm[1] - kinematicViscosity*div(gradV + gradVT)[1]))
print("storage - viscous terms 2", sp.simplify(storageTerm[2] - kinematicViscosity*div(gradV + gradVT)[2]))

print("pressure term balances advection term")
print("advection term + pressure term 0", sp.simplify(div(vvT) + div(pI))[0])
