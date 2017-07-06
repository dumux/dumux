#!/bin/sh
make test_1pmimeticanisotropic

outName=$1
inputArgs=""
./test_1pmimeticanisotropic  -ParameterFile ./test_1pmimeticanisotropic_case3.input -Problem.Name $outName $inputArgs -Problem.TestCase 3 >> $1".out"
