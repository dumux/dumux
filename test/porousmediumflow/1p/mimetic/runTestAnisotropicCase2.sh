#!/bin/sh
make test_1pmimeticanisotropic
for i in {0..6}
do
   outName=$1
   gridName=$2$i
   gridFile="$gridName.dgf"
   inputArgs=""
  ./test_1pmimeticanisotropic  -ParameterFile ./test_1pmimeticanisotropic_case2.input -Grid.File ../implicit/grids/$gridFile -Problem.Name $outName $inputArgs -Problem.TestCase 2 >> $1".out"
done