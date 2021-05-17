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
x, y, t, mu = sp.symbols("x y t mu")

kinematicViscosity = 1.0

# analytical solutions for v and p in free-flow domain (see reference given in problems)
def analyticalSolutionStokes():
    psi = sp.exp(-5*mu*sp.pi**2*t)*sp.cos(sp.pi*x)*sp.cos(2*sp.pi*y)
    ux = sp.diff(psi,y)
    uy = -sp.diff(psi,x)

    vFF = np.array([- 2.0 * sp.pi * sp.exp(- 5.0 * mu * sp.pi * sp.pi * t) * sp.cos(sp.pi * x) * sp.sin(2.0 * sp.pi * y), sp.pi * sp.exp(- 5.0 * mu * sp.pi*sp.pi * t) * sp.sin(sp.pi * x) * sp.cos(2.0 * sp.pi * y)])

    #assert(vFF[0]-ux == 0)
    #assert(vFF[1]-uy == 0)

    pFF = - 0.25 * sp.exp(-10.0 * mu * sp.pi**2 * t) * sp.pi**2 * (4.0 * sp.cos(2.0 * sp.pi * x) + sp.cos(4.0 * sp.pi * y))
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
momentumFlux = vvT - mu*(gradV + gradVT) +pI
divMomentumFlux = div(momentumFlux)

# solution for source term
print("Source term mass:", sp.diff(vFF[0],x) + sp.diff(vFF[1], y))
print("Source term momentum x:", sp.simplify(divMomentumFlux[0] + storageTerm[0]))
print("Source term momentum y:", sp.simplify(divMomentumFlux[1] + storageTerm[1]))




print('\n##### Integration over control volume')
fluxBottom = sp.integrate(vFF.dot(np.array([0, -1])).subs(y, 0), (x, 0, 0.1))
fluxTop= sp.integrate(vFF.dot(np.array([0, 1])).subs(y, 0.1), (x, 0, 0.1))
print("fluxBottom", sp.simplify(fluxBottom))
print("fluxTop", sp.simplify(fluxTop))

print("Sum vert", sp.simplify(fluxBottom+fluxTop))

fluxLeft= sp.integrate(vFF.dot(np.array([-1, 0])).subs(x, 0), (y, 0, 0.1))
fluxRight= sp.integrate(vFF.dot(np.array([1, 0])).subs(x, 0.1), (y, 0, 0.1))

print("fluxLeft", sp.simplify(fluxLeft))
print("fluxRight", sp.simplify(fluxRight))

print("Total Sum", sp.simplify(fluxBottom+fluxTop+fluxLeft+fluxRight))



print("Storage cancels diffusion: storage-diff:", sp.simplify(storageTerm-div(mu*(gradV + gradVT))))
print("Advection cancels pressure: Advection+press:", sp.simplify(div(vvT)+div(pI)))

print('v interior at 5e-7', vFF[1].subs(mu, 0.1).subs(x, 0.01).subs(y,0.02).subs(t,0.5e-6).evalf())
print('v boundary at 5e-7', vFF[1].subs(mu, 0.1).subs(x, 0.01).subs(y,0.00).subs(t,0.5e-6).evalf())
print('v boundary at 1e-6', vFF[1].subs(mu, 0.1).subs(x, 0.01).subs(y,0.00).subs(t,1e-6).evalf())
print('v boundary left at 1e-6', vFF[0].subs(mu, 0.1).subs(x, 0.00).subs(y,0.01).subs(t,1e-6).evalf())
print('p at 0.5e-6', pI[0][0].subs(mu, 0.1).subs(x, 0.01).subs(y,0.01).subs(t,0.5e-6).evalf())
print('p at 1e-6', pI[0][0].subs(mu, 0.1).subs(x, 0.01).subs(y,0.01).subs(t,1e-6).evalf())


volume = 0.02**2

print("Eval at x = 0.02, y= 0.77, t = 1e-6")
print("Storage term", storageTerm[0].subs(mu, 0.1).subs(x, 0.02).subs(y,0.77).subs(t,1e-6).evalf()*volume)

v0 = vFF[0].subs(mu, 0.1).subs(x, 0.02).subs(y,0.77).subs(t,0.5e-6).evalf()
v1 = vFF[0].subs(mu, 0.1).subs(x, 0.02).subs(y,0.77).subs(t,1e-6).evalf()

print("v1", v1)
print("v0", v0)
print("v(0)", vFF[0].subs(mu, 0.1).subs(x, 0.02).subs(y,0.77).subs(t,0).evalf())
#print("v(1e-6)", vFF[0].subs(mu, 0.1).subs(x, 0.02).subs(y,0.77).subs(t,1e-6).evalf()


print("Storage term, numerical", (v1-v0)/0.5e-6*volume)
print("diffusion term", div(mu*(gradV + gradVT))[0].subs(mu, 0.1).subs(x, 0.02).subs(y,0.77).subs(t,1e-6).evalf()*volume)

#print("\nNeumann fluxes ######################")
##print("Shape momentumFlux", momentumFlux.shape)
##momFluxNormal = momentumFlux.dot(np.array([1,0]))
##print("Shape momFluxNormal", momFluxNormal.shape)
#neumann = sp.simplify(momentumFlux)
#print(" momFlux 0 0", neumann[0][0])
#print(" momFlux 0 1", neumann[0][1])
#print(" momFlux 1 0", neumann[1][0])
#print(" momFlux 1 1", neumann[1][1])
##print("momentumFlux * normal", sp.simplify(momFluxNormal))
#print("######################\n")





#print("Sol at ccc", vFF[0].subs(t,0))

#intStor = sp.simplify(sp.integrate(storageTerm[1].subs(x,0.01).subs(y,0.04), (t, 0, 1e-8)))

#print("integralStorageTerm 1e-8:", intStor.evalf() * 0.0002)
#print("test:", sp.simplify(sp.integrate(vFF[1].subs(x,0.01).subs(y,0.04), (t, 0, 1e-8))).evalf() * 0.0002)


#numericStor0 = sp.simplify(vFF[1].subs(x,0.03).subs(y,0.02).subs(t,0))
#numericStor1 = sp.simplify(vFF[1].subs(x,0.03).subs(y,0.02).subs(t,1e-8))

#numericStor = (numericStor1 - numericStor0).evalf() / 1e-8 * 0.0002

#print("numeric storage", numericStor)



##intStor = sp.simplify(sp.integrate(storageTerm[1].subs(x,0.01).subs(y,0.04), (t, 0, 1e-3)))

##print("integralStorageTerm 1e-8:", intStor.evalf() * 0.0002)
##print("test:", sp.simplify(sp.integrate(vFF[1].subs(x,0.01).subs(y,0.04), (t, 0, 1e-3))).evalf() * 0.0002)


##numericStor0 = sp.simplify(vFF[1].subs(x,0.01).subs(y,0.04).subs(t,0))
##numericStor1 = sp.simplify(vFF[1].subs(x,0.01).subs(y,0.04).subs(t,1e-3))

##numericStor = (numericStor1 - numericStor0).evalf() / 1e-3 * 0.0002

##print("numeric storage", numericStor)

#print("time deriv at 1e-8", storageTerm[1].subs(x,0.03).subs(y,0.02).subs(t, 1e-8).evalf() * 0.0002)
