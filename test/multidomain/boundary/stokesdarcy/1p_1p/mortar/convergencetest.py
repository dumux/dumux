import os
import math
import subprocess
import numpy as np
import matplotlib.pyplot as plt

errorFileStandard = "errors_standard.txt"
errorFileSharp = "errors_sharp.txt"

# remove error norm files from previous runs
if os.path.exists(errorFileStandard):
    subprocess.call(["rm", errorFileStandard])
if os.path.exists(errorFileSharp):
    subprocess.call(["rm", errorFileSharp])

cells1 = [11, 10]
cells2 = [10, 10]
cellsMortar = 4
numRefinements = 4

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
    subprocess.call(["./test_md_mortar_darcy1p_stokes1p",
                     "-L2Error.OutputFile", errorFileStandard,
                     "-Domain1.Grid.Cells", strCells1,
                     "-Domain1.Problem.Name", "domain1",
                     "-Domain2.Grid.Cells", strCells2,
                     "-Domain2.Problem.Name", "domain2",
                     "-Mortar.VariableType", "Flux",
                     "-Mortar.Grid.Cells", strCellsMortar])

    print("\n\nCalling sharp solver for refinement " + str(refIdx-1))
    subprocess.call(["./sharp/test_md_mortar_sharp_darcy1p_stokes1p", "params.input",
                     "-L2Error.OutputFile", errorFileSharp,
                     "-Domain1.Grid.Cells", strCells1,
                     "-Domain1.Problem.Name", "domain1_sharp",
                     "-Domain2.Grid.Cells", strCells2,
                     "-Domain2.Problem.Name", "domain2_sharp",
                     "-Mortar.Grid.Cells", strCellsMortar])

    # add numCells for this level
    hInverse.append( (0.5*(cells1X + cells2X)) )

errorsStandard = np.loadtxt(errorFileStandard, delimiter=',')
errorsSharp = np.loadtxt(errorFileSharp, delimiter=',')

# compute rates
ratesStandardPressure = []
ratesStandardFlux = []
ratesStandardMortar = []
ratesStandardIFFlux = []
ratesSharpPressure = []
ratesSharpFlux = []
ratesSharpMortar = []
ratesSharpIFFlux = []

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
    deltaEP = np.log(errorsSharp[i+1][0]) - np.log(errorsSharp[i][0])
    deltaEF = np.log(errorsSharp[i+1][1]) - np.log(errorsSharp[i][1])
    deltaEM = np.log(errorsSharp[i+1][2]) - np.log(errorsSharp[i][2])
    deltaEIF = np.log(errorsSharp[i+1][3]) - np.log(errorsSharp[i][3])
    ratesSharpPressure.append( deltaEP / deltaH )
    ratesSharpFlux.append( deltaEF / deltaH )
    ratesSharpMortar.append( deltaEM / deltaH )
    ratesSharpIFFlux.append( deltaEIF / deltaH )

print("hInverse: ", hInverse)
print("Flat pressure: ", ratesStandardPressure)
print("Flat flux: ", ratesStandardFlux)
print("Flat mortar: ", ratesStandardMortar)
print("Flat ifFlux: ", ratesStandardIFFlux)
print("Sharp pressure: ", ratesSharpPressure)
print("Sharp flux: ", ratesSharpFlux)
print("Sharp mortar: ", ratesSharpMortar)
print("Sharp ifFlux: ", ratesSharpIFFlux)

# plot pressure error norms
plt.figure(1)
plt.loglog(hInverse, errorsStandard[:, 0], label='standard')
plt.loglog(hInverse, errorsSharp[:, 0], label='sharp')
plt.legend()
plt.xlabel("1/h")
plt.ylabel("absolute error")
plt.savefig("pressure.pdf", bbox_inches='tight')

# plot flux error norms
plt.figure(2)
plt.loglog(hInverse, errorsStandard[:, 1], label='standard')
plt.loglog(hInverse, errorsSharp[:, 1], label='sharp')
plt.legend()
plt.xlabel("1/h")
plt.ylabel("absolute error")
plt.savefig("flux.pdf", bbox_inches='tight')

# plot flux error norms
plt.figure(3)
plt.loglog(hInverse, errorsStandard[:, 2], label='standard')
plt.loglog(hInverse, errorsSharp[:, 2], label='sharp')
plt.legend()
plt.xlabel("1/h")
plt.ylabel("absolute error")
plt.savefig("fluxmortar.pdf", bbox_inches='tight')
