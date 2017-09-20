#!usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

# function to calculate rho dependent on pressure
rho_min = 1440;
rho_max = 1480;
k = 5e-7;

def rho(p):
    return rho_min + (rho_max - rho_min)/(1 + rho_min*np.exp(-1.0*k*(rho_max - rho_min)*p));

# sample pressure in range (1e4, 1e7) and compute corresponding densities
p = np.logspace(4, 7, 100)
r = rho(p)

# plot density vs. pressure
plt.semilogx(p, r)
plt.show()
