#!/bin/sh
make test_1pmimeticanisotropic
for i in {0..6}
do
   outName=$1
   gridName=$2$i
   gridFile="$gridName.dgf"
   inputArgs=""
  ./test_1pmimeticanisotropic  -Grid.File ../implicit/grids/$gridFile -Problem.Name $outName $inputArgs -Problem.TestCase 1 >> $1".out"
done
