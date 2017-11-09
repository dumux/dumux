#!/bin/sh
make test_spe10_tpfa

#inputArgs="-Newton.UseLineSearch 1"
#outputName="test_spe10_tpfa_superlu"

inputArgs="-Newton.UseLineSearch 1 -LinearSolver.Verbosity 1 -LinearSolver.MaxIterations 2000"
outputName="test_spe10_tpfa_ilu0bicgstab"

echo $outputName

for i in {0,10,20,30,40,50,60,70,80}
do
   outName=$outputName"_layer"$i
   echo $outName
   rm $outName".out"
  ./test_spe10_tpfa -Problem.Name $outName $inputArgs -SpatialParams.LayerIdx $i > $outName".out"
done
