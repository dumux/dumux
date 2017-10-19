#!/bin/sh
name=solutionFor_1pweakconv_harmPoints_case1_mimetic
rm $name"_rates.txt"
rm $name".out"
sh ./runTestWeakConv.sh $name 1

name=solutionFor_1pweakconv_harmPoints_case1_mpfal
rm $name"_rates.txt"
rm $name".out"
sh ./runTestWeakConvMPFAL.sh $name 1

############################################################
############################################################

name=solutionFor_1pweakconv_harmPoints_case2_mimetic
rm $name"_rates.txt"
rm $name".out"
sh ./runTestWeakConv.sh $name 2


name=solutionFor_1pweakconv_harmPoints_case2_mpfal
rm $name"_rates.txt"
rm $name".out"
sh ./runTestWeakConvMPFAL.sh $name 2

############################################################
############################################################

name=solutionFor_1pweakconv_harmPoints_case3_mimetic
rm $name"_rates.txt"
rm $name".out"
sh ./runTestWeakConv.sh $name 3

name=solutionFor_1pweakconv_harmPoints_case3_mpfal
rm $name"_rates.txt"
rm $name".out"
sh ./runTestWeakConvMPFAL.sh $name 3
