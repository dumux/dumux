#!/bin/sh
make test_1pweakconv
for i in {3..7}
do
   outName=$1
   inputArgs="-Grid.Refinement $i"
  ./test_1pweakconv -Problem.Name $outName $inputArgs -Problem.TestCase $2 >> $1".out"
done
