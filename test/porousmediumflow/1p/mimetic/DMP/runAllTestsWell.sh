#!/bin/sh

outName=solutionFor_well_mimetic
rm $outName".out"

make test_1pdmpwell_mimetic

./test_1pdmpwell_mimetic  -Problem.Name $outName -Problem.TestCase 2 >> $outName".out"

#######################################################################################

outName=solutionFor_well_mpfao
rm $outName".out"

make test_1pdmpwell_mpfao

./test_1pdmpwell_mpfao  -Problem.Name $outName -Problem.TestCase 2 >> $outName".out"

#######################################################################################
