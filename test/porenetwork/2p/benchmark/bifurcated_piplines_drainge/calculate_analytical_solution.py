import sys
sys.path.append('../1d_pipeline_drainage/')
from plot_PcS_and_derivatives import *
import numpy as np


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
    pcEntry = np.genfromtxt("pcEntry.csv", skip_header=True).T
    poreVolume = 8*defaultPoreRadius*defaultPoreRadius*defaultPoreRadius
    volumeFlux = 5e-10/1000
    swEntry = Sw(pcEntry)
    # print(swEntry)
    path = find_path(15, swEntry)
    print(path)
    time = calculate_invasion_time(path, swEntry, poreVolume, volumeFlux)
