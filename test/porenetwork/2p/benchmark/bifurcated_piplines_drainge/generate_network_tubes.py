"""
This is a script used to generate column-like pore-networks in .dgf file.
Column-like means for each column the pores have the same y-coordinates.
"""

import numpy as np
import random as rd


def poreCoordinates(poreNumberX, poreNumberY, fixThroatLength_Y, fixPoreRadius_Y, distanceX):
    poreX_unit = []
    poreY_unit = []
    for i in range(poreNumberX):
        poreX_unit.append(i * distanceX)
    poreX = np.repeat(poreX_unit, poreNumberY)
    distance_Y = 2 * fixPoreRadius_Y + fixThroatLength_Y
    for i in range(poreNumberY):
        poreY_unit.append((i*distance_Y))
    poreY = np.array(list(poreY_unit) * poreNumberX)
    poreZ = np.zeros( (poreNumberX* poreNumberY))
    poreX = np.append(poreX, -distanceX)
    poreY = np.append(poreY, 0.0)
    poreZ = np.append(poreZ, 0.0)
    pore_Coordinates = [poreX, poreY, poreZ]
    return pore_Coordinates


def writeDGF(filename, vData, eData):
    # write DGF file
    with open(filename, "w") as outputfile:
        outputfile.write("DGF\n")
        outputfile.write(
            "% Vertex parameters: PoreInscribedRadius PoreLabel"
            + "\n"
        )
        outputfile.write(
            "% Element parameters: ThroatInscribedRadius ThroatLength ThroatLabel"
            + "\n"
        )
        outputfile.write("Vertex\n")
        outputfile.write("parameters " + str(vData.shape[1] - 3) + "\n")  # vertex data, -3 because first 3 parameters are coordinates, shape[1] is column number
        for i in range(len(vData)):
            outputfile.write(" ".join([str(v) for v in vData[i]]) + "\n")
        outputfile.write("\n#\n")
        outputfile.write("SIMPLEX\n")
        outputfile.write("parameters " + str(3) + "\n")  # cell data, first 2 parameters are pb to pth links, shape[0] and shape[1] are row and column number
        for i in range(len(eData)):
            outputfile.write(
                " ".join(
                    [str(v if idx > 1 else int(v)) for idx, v in enumerate(eData[i])]
                )
                + "\n"
            )
        outputfile.write("\n#\n")
        outputfile.write("BOUNDARYDOMAIN\ndefault 1\n#")

poreNumberX = 4
poreNumberY = 4

pore_Coordinates = poreCoordinates(poreNumberX, poreNumberY, 0.1*1e-3, 0.2*1e-3, 0.5*1e-3)

print(pore_Coordinates[0])
print(pore_Coordinates[1])
print(pore_Coordinates[2])

numberOfPores = poreNumberX * poreNumberY
poreLabel = np.ones(numberOfPores)*(-1)
poreRadius = np.ones(numberOfPores + 1)*(0.2*1e-3)
# 3 outlet
for i in range(poreNumberX):
    poreLabel[(i + 1) * poreNumberY - 1] = 3
poreLabel = np.append(poreLabel, 2)
vData = np.array( [pore_Coordinates[0], pore_Coordinates[1], pore_Coordinates[2], poreRadius, poreLabel] ).T


throat1 = []
throat2 = []

for i in range(poreNumberX):
    for j in range(poreNumberY - 1):
        throat1.append(poreNumberY * i + j)
        throat2.append(poreNumberY * i + j + 1)

for i in range(poreNumberX - 1):
    throat1.append(poreNumberY * i)
    throat2.append(poreNumberY * (i +1))


numberOfThroats = len(throat1)
throatLength = np.ones(numberOfThroats) * (0.1*1e-3)
throatLabel = np.ones(numberOfThroats) * (-1)

for i in range(poreNumberX):
    throatLabel[(poreNumberY - 2) + (poreNumberY - 1)* i] = 3

throatRadius = []
for i in range(numberOfThroats):
    throatRadius.append(1e-3 * rd.uniform(0.1, 0.15))
throat1.append(numberOfPores)
throat2.append(0)
throatLength = np.append(throatLength, 0.1*1e-3)
throatRadius.append(0.18*1e-3)
throatLabel = np.append(throatLabel, 2)
eData = np.array([throat1, throat2, throatRadius, throatLength, throatLabel]).T

writeDGF("capillary_multiple_tubes.dgf", vData, eData)
