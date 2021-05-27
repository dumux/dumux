import os
import math
import subprocess
import numpy as np
# import matplotlib.pyplot as plt

errorFileStandard = "errors_standard.txt"

# remove error norm files from previous runs
if os.path.exists(errorFileStandard):
    subprocess.call(["rm", errorFileStandard])

cells1 = [10, 10]
cells2 = [11, 10]
cellsMortar = 8
numRefinements = 2

hInverse = []
for refIdx in range(1, numRefinements+2):
    factor = math.pow(2, refIdx-1)
    cells1X = int(cells1[0]*factor)
    cells1Y = int(cells1[1]*factor)
    cells2X = int(cells2[0]*factor)
    cells2Y = int(cells2[1]*factor)

    strCells1 = str(cells1X) + " " + str(cells1Y)
    strCells2 = str(cells2X) + " " + str(cells2Y)
    strCellsMortar = str(int(cellsMortar*factor))

    # run schemes
    print("\n\nCalling standard solver for refinement " + str(refIdx-1))
    subprocess.call(["./test_md_mortar_darcy1p_stokes1p_convergence",
                     "-L2Error.OutputFile", errorFileStandard,
                     "-Domain1.Grid.Cells", strCells1,
                     "-Domain1.Problem.Name", "domain1",
                     "-Domain2.Grid.Cells", strCells2,
                     "-Domain2.Problem.Name", "domain2",
                     "-Mortar.VariableType", "Pressure",
                     "-Mortar.Grid.Cells", strCellsMortar])

    # add numCells for this level
    hInverse.append( (0.5*(cells1X + cells2X)) )

errorsStandard = np.loadtxt(errorFileStandard, delimiter=',')

# compute rates
ratesStandardPressure = []
ratesStandardFlux = []
ratesStandardMortar = []
ratesStandardIFFlux = []

for i in range(0, len(errorsStandard)-1):
    deltaEP = np.log(errorsStandard[i+1][0]) - np.log(errorsStandard[i][0])
    deltaEF = np.log(errorsStandard[i+1][1]) - np.log(errorsStandard[i][1])
    deltaEM = np.log(errorsStandard[i+1][2]) - np.log(errorsStandard[i][2])
    deltaEIF = np.log(errorsStandard[i+1][3]) - np.log(errorsStandard[i][3])
    deltaH = np.log(hInverse[i+1]) - np.log(hInverse[i])
    ratesStandardPressure.append( deltaEP / deltaH )
    ratesStandardFlux.append( deltaEF / deltaH )
    ratesStandardMortar.append( deltaEM / deltaH )
    ratesStandardIFFlux.append( deltaEIF / deltaH )

print("hInverse: ", hInverse)
print("Flat pressure: ", ratesStandardPressure)
print("Flat flux: ", ratesStandardFlux)
print("Flat mortar: ", ratesStandardMortar)
print("Flat ifFlux: ", ratesStandardIFFlux)
#
# # plot pressure error norms
# plt.figure(1)
# plt.loglog(hInverse, errorsStandard[:, 0], label='standard')
# plt.loglog(hInverse, errorsSharp[:, 0], label='sharp')
# plt.legend()
# plt.xlabel("1/h")
# plt.ylabel("absolute error")
# plt.savefig("pressure.pdf", bbox_inches='tight')
#
# # plot flux error norms
# plt.figure(2)
# plt.loglog(hInverse, errorsStandard[:, 1], label='standard')
# plt.loglog(hInverse, errorsSharp[:, 1], label='sharp')
# plt.legend()
# plt.xlabel("1/h")
# plt.ylabel("absolute error")
# plt.savefig("flux.pdf", bbox_inches='tight')
#
# # plot flux error norms
# plt.figure(3)
# plt.loglog(hInverse, errorsStandard[:, 2], label='standard')
# plt.loglog(hInverse, errorsSharp[:, 2], label='sharp')
# plt.legend()
# plt.xlabel("1/h")
# plt.ylabel("absolute error")
# plt.savefig("fluxmortar.pdf", bbox_inches='tight')
