#!/bin/sh
outName=solutionFor_anisotropic_case3_mimetic
rm $outName".out"

make test_1panisotropic_mimetic

./test_1panisotropic_mimetic  -ParameterFile ./test_1panisotropic_mimetic_case3.input -Problem.Name $outName -Problem.TestCase 3 >> $outName".out"

#######################################################################################

outName=solutionFor_anisotropic_case3_mpfao
rm $outName".out"

make test_1panisotropic_mpfao

./test_1panisotropic_mpfao  -ParameterFile ./test_1panisotropic_mpfao_case3.input -Problem.Name $outName -Problem.TestCase 3 >> $outName".out"

#######################################################################################
