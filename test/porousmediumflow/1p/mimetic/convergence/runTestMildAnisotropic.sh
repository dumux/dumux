#!/bin/sh
make test_1panisotropic
for i in {0..6}
do
   outName=$1
   gridName=$2$i
   gridFile="$gridName.dgf"
   inputArgs="-Mpfa.UseNlTpfa $3 -Mpfa.UseMFDAnsatz 1 -Mpfa.ObjectiveFunction 1 -Mpfa.UseTpfa 0 -Mpfa.NewBoundaryApprox 1 -Mpfa.ComplexBoundary 0 -SpatialParams.k1 1 -SpatialParams.k2 1 -SpatialParams.k12 0.5"
  ./test_1panisotropic  -Grid.File ./grids/$gridFile -Problem.Name $outName $inputArgs -Problem.TestCase 1 >> $1".out"
done