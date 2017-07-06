#!/bin/sh
name=solutionFor_2dNonmatchingCubic_anisotropic_case2_mimetic
rm $name"_rates.txt"
rm $name".out"
sh ./runTestAnisotropicCase2.sh $name 2dNonmatchingCubic_

name=solutionFor_2dRandCubic_anisotropic_case2_mimetic
rm $name"_rates.txt"
rm $name".out"
sh ./runTestAnisotropicCase2.sh $name 2dRandCubic_

name=solutionFor_2dTwistedCubic_anisotropic_case2_mimetic
rm $name"_rates.txt"
rm $name".out"
sh ./runTestAnisotropicCase2.sh $name 2dTwistedCubic_