import numpy as np
from scipy.spatial.distance import cdist  #for eulidian distance between pairs of coordinates
import itertools as it

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
        outputfile.write("parameters " + str(3) + "\n")  # cell data, -2 because first 2 parameters are pb to pth links, shape[0] and shape[1] are row and column number
        for i in range(len(eData)):
            outputfile.write(
                " ".join(
                    [str(v if idx > 1 else int(v)) for idx, v in enumerate(eData[i])]
                )
                + "\n"
            )
        outputfile.write("\n#\n")
        outputfile.write("BOUNDARYDOMAIN\ndefault 1\n#")


###############################################
####### Parameter to generate the grid ########
###############################################
rowNumber = 5               # 5 row
colNumber = 10              # 10 column
distance = 5e-4             # 0.5mm between center of two pore bodies
pR = 2e-4                   # 0.2mm poreRadius
tL = distance - 2.0 * pR    # throat length



######################################################
poreX = []
poreY = []
poreLabel = []
for i in range(rowNumber):
    for j in range(colNumber):
        poreX.append(j * distance)
        poreY.append(i * distance)


pore_Number = len(poreX)
assert pore_Number == rowNumber*colNumber, "Pore number is not consistent! Please check the network geometry!"
poreRadius = np.ones(pore_Number) * pR
poreX = np.array(poreX)*1e-3
poreY = np.array(poreY)*1e-3
poreZ = np.zeros(pore_Number)
poreRadius = np.array(poreRadius)

for i in range(pore_Number):
    if poreX[i] == 0:
        poreLabel_i = 2
    elif poreX[i] == (colNumber-1) * distance * 1e-3:
        poreLabel_i = 3
    else:
        poreLabel_i = -1
    poreLabel.append(poreLabel_i)

vData = np.array([poreX, poreY, poreZ, poreRadius, poreLabel]).T

coordinates = np.array([poreX, poreY, poreZ]).T

maxDistance = distance + 1e-8

throat_indices_0 = []
throat_indices_1 = []

for i in range(rowNumber):
    for j in range(colNumber - 1):
        throat_indices_0.append(i*colNumber + j)
        throat_indices_1.append(i*colNumber + j + 1)

for i in range(rowNumber - 1):
    for j in range(colNumber -2):
        throat_indices_0.append(i * colNumber + j + 1)
        throat_indices_1.append( (i+1) * colNumber + j + 1)

throat_number = len(throat_indices_0)
throatLength = np.ones(throat_number) * tL
throatRadius   = np.random.normal(pR*0.5, pR*0.1, throat_number)
#throatRadius = np.ones(throat_number) * pR * 0.5
throatLabel = np.ones(throat_number) * -1
eData = np.array([throat_indices_0, throat_indices_1, throatRadius, throatLength, throatLabel]).T

writeDGF("lattic_network.dgf", vData, eData)
