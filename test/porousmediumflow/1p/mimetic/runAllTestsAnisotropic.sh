#!/bin/sh
name=solutionFor_2dNonmatchingCubic_anisotropic_case1_mimetic
rm $name"_rates.txt"
rm $name".out"
sh ./runTestAnisotropic.sh $name 2dNonmatchingCubic_  1


name=solutionFor_2dRandCubic_anisotropic_case1_mimetic
rm $name"_rates.txt"
rm $name".out"
sh ./runTestAnisotropic.sh $name 2dRandCubic_  1

name=solutionFor_2dTwistedCubic_anisotropic_case1_mimetic
rm $name"_rates.txt"
rm $name".out"
sh ./runTestAnisotropic.sh $name 2dTwistedCubic_  1
