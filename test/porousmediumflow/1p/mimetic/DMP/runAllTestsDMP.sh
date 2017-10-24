#!/bin/sh

outName=solutionFor_DMP_mimetic
rm $outName".out"

make test_1pdmp_mimetic

./test_1pdmp_mimetic  -Problem.Name $outName -Problem.TestCase 1 >> $outName".out"

#######################################################################################

outName=solutionFor_DMP_mpfao
rm $outName".out"

make test_1pdmp_mpfao

./test_1pdmp_mpfao  -Problem.Name $outName -Problem.TestCase 1 >> $outName".out"

#######################################################################################
