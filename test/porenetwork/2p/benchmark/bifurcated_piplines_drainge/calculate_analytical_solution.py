"""
This is a script used to generate column-like pore-networks in .dgf file,
and then to calculate the analytical solution of the invading path.
"""

import sys
sys.path.append('../1d_pipeline_drainage/')
from plot_PcS_and_derivatives import *
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

def find_path(startingThroat, swEntry):
    neighbor = [None]*16
    neighbor[0] = [15, 1, 12]
    neighbor[1] = [0, 2]
    neighbor[2] = [1]
    neighbor[3] = [12, 3, 14]
    neighbor[4] = [3, 5]
    neighbor[5] = [4]
    neighbor[6] = [13, 7, 14]
    neighbor[7] = [6, 8]
    neighbor[8] = [7]
    neighbor[9] = [14, 10]
    neighbor[10] = [9, 11]
    neighbor[11] = [10]
    neighbor[12] = [15, 0, 3, 13]
    neighbor[13] = [12, 3, 6, 14]
    neighbor[14] = [13, 6, 9]
    neighbor[15] = [0, 12]

    allthroats = np.linspace(0, 15, 16, dtype=int)
    endThroats = [2, 5, 8, 11]
    throatsInvaded = []
    throatNotInvaded = np.setdiff1d(allthroats, throatsInvaded)
    # print(throatNotInvaded)

    ongoingthroat = startingThroat
    throatsInvaded.append(ongoingthroat)
    while (ongoingthroat not in endThroats):
        allneighbors = [neighbor[i][j] for i in throatsInvaded for j in range(len(neighbor[i]))]
        neighborNotInvaded = np.setdiff1d(allneighbors, throatsInvaded)
        # print("neighbors not invaded yet: ", neighborNotInvaded)
        neighborNotInvadedSw = [ swEntry[i] for i in neighborNotInvaded ]
        ongoingthroat =  list(swEntry).index(max(neighborNotInvadedSw))
        # print("throat currently invaded: ", ongoingthroat)
        throatsInvaded.append(ongoingthroat)
        # print("throat has been invaded: ", throatsInvaded)
    return throatsInvaded

def calculate_invasion_time(path, swEntry, poreVolume, volumeFlux):
    t = []
    t_inv = []
    for i, throat in enumerate(path):
        t.append( (1.0 - swEntry[throat]) * poreVolume / volumeFlux)
        t_inv.append( sum(t) )
        print(swEntry[throat])
    # print("time used for each inv: ", t)
    print("invasion time point: ", t_inv)

if __name__ == "__main__":

    poreNumberX = 4
    poreNumberY = 4

    pore_Coordinates = poreCoordinates(poreNumberX, poreNumberY, 0.1*1e-3, 0.2*1e-3, 0.5*1e-3)

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

    pcEntry = PcEntry(radius=np.array(throatRadius))
    poreVolume = 8*defaultPoreRadius*defaultPoreRadius*defaultPoreRadius
    volumeFlux = 5e-10/1000
    swEntry = Sw(pcEntry)
    # print(swEntry)
    path = find_path(15, swEntry)
    print(path)
    time = calculate_invasion_time(path, swEntry, poreVolume, volumeFlux)
