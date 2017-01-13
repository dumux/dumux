#! /usr/bin/python

import subprocess
from math import log, floor, fabs
import numpy

radius = 0.05
# # maximum allowed refinement is numCells = 2/radius then the average operator has to be implemented
# # only even numbers allowed otherwise log(0) will be evaluated
# levels = [int(i) for i in numpy.linspace(7, 21, 5)]
# levels = [(i+1 if i%2 == 0 else i) for i in levels]

# #initialize log file
# file = open("log.txt", "w+")
# file.write("#l2-norms\n")
# file.close()

# for numCells in levels:
#     subprocess.call(['./dangelo', '-VesselGrid.Cells', str(numCells),
#                                   '-TissueGrid.Cells', " ".join([str(numCells)]*3),
#                                   '-SpatialParams.Radius', str(radius)])

# calculate the rates
hmax3 = []
hmax1 = []
p3 = []
p1 = []
for line in open("log.txt"):
    line = line.strip()
    if line.startswith('#'):
        continue
    list = line.split(' ')
    hmax3.append(float(list[0]))
    hmax1.append(float(list[1]))
    p3.append(float(list[2]))
    p1.append(float(list[3]))

rates3d = []
rates1d = []
for i in range(1, len(p1)):
    try:
        rate1 = (log(p1[i-1]) - log(p1[i]))/(log(hmax1[i-1]) - log(hmax1[i]))
        rate3 = (log(p3[i-1]) - log(p3[i]))/(log(hmax3[i-1]) - log(hmax3[i]))
    except ValueError:
        continue
    print('1d: {0:.5f}  3d: {1:.5f}'.format(rate1, rate3))
    rates3d.append(rate3)
    rates1d.append(rate1)


# plot the results
import matplotlib.pyplot as plt
plt.grid("on")
plt.plot(hmax1, p1, label='|pv-pv^h|L2')
plt.plot(hmax3, p3, label='|pt-pt^h|L2')
plt.plot(hmax1, [h*fabs(log(h)) for h in hmax1], label='h') #label='h*|log h|')
plt.plot(hmax3, [h*fabs(log(h)) for h in hmax3], label='h') #label='h*|log h|')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()
