#!/bin/sh
make test_spe10

inputArgs="-Newton.UseLineSearch 0
           -Newton.EnableChop 1 -TimeManager.DtInitial 5e6  -TimeManager.MaxTimeStepSize 5e6
           -Newton.MaxSteps 100 -Newton.TargetSteps 100"
outputName="test_spe10_mixedmimetic_superlu_chop_constDt"

#inputArgs="-Newton.UseLineSearch 0
#           -Newton.EnableChop 1 -LinearSolver.Verbosity 1 -LinearSolver.MaxIterations 2000 -TimeManager.DtInitial 5e6  -TimeManager.MaxTimeStepSize 5e6
#           -Newton.MaxSteps 100 -Newton.TargetSteps 100"
#outputName="test_spe10_mimetic_ilu0bicgstab_chop_constDt"

echo $outputName

for i in {0,10,20,30,40,50,60,70,80}
do
   outName=$outputName"_layer"$i
   echo $outName
   rm $outName".out"
  ./test_spe10 -Problem.Name $outName $inputArgs -SpatialParams.LayerIdx $i > $outName".out"
done
