#!/bin/sh

make test1_conv_mimetic

inputArgs="-Newton.UseLineSearch 0"

outputName="test1_fvca6_mimetic_nonconvex"
echo $outputName
rm $outputName"_rates.txt"
rm $outputName".out"
./test1_conv_mimetic -Grid.File ../grids/nonconvex_4.dgf -Problem.Name $outputName $inputArgs > $outputName".out"
./test1_conv_mimetic -Grid.File ../grids/nonconvex_8.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test1_conv_mimetic -Grid.File ../grids/nonconvex_16.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test1_conv_mimetic -Grid.File ../grids/nonconvex_32.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test1_conv_mimetic -Grid.File ../grids/nonconvex_64.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"

outputName="test1_fvca6_mimetic_cossurf"
echo $outputName
rm $outputName"_rates.txt"
rm $outputName".out"
./test1_conv_mimetic -Grid.File ../grids/nonconvex_cossurf_4.dgf -Problem.Name $outputName $inputArgs > $outputName".out"
./test1_conv_mimetic -Grid.File ../grids/nonconvex_cossurf_8.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test1_conv_mimetic -Grid.File ../grids/nonconvex_cossurf_16.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test1_conv_mimetic -Grid.File ../grids/nonconvex_cossurf_32.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test1_conv_mimetic -Grid.File ../grids/nonconvex_cossurf_64.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"

outputName="test1_fvca6_mimetic_rand"
echo $outputName
rm $outputName"_rates.txt"
rm $outputName".out"
./test1_conv_mimetic -Grid.File ../grids/random_4.dgf -Problem.Name $outputName $inputArgs > $outputName".out"
./test1_conv_mimetic -Grid.File ../grids/random_8.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test1_conv_mimetic -Grid.File ../grids/random_16.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test1_conv_mimetic -Grid.File ../grids/random_32.dgf -Problem.Name $outputName  $inputArgs >> $outputName".out"
./test1_conv_mimetic -Grid.File ../grids/random_64.dgf -Problem.Name $outputName  $inputArgs >> $outputName".out"

outputName="test1_fvca6_mimetic_checkerboard"
echo $outputName
rm $outputName"_rates.txt"
rm $outputName".out"
./test1_conv_mimetic -Grid.File ../grids/test_1p_3d.dgf -Problem.Name $outputName -Grid.UseCheckerboard 1 $inputArgs > $outputName".out"
./test1_conv_mimetic -Grid.File ../grids/test_1p_3d.dgf -Problem.Name $outputName -Grid.Refinement 1 -Grid.UseCheckerboard 1  $inputArgs >> $outputName".out"
./test1_conv_mimetic -Grid.File ../grids/test_1p_3d.dgf -Problem.Name $outputName -Grid.Refinement 2 -Grid.UseCheckerboard 1 $inputArgs >> $outputName".out"
./test1_conv_mimetic -Grid.File ../grids/test_1p_3d.dgf -Problem.Name $outputName -Grid.Refinement 3 -Grid.UseCheckerboard 1 $inputArgs >> $outputName".out"
./test1_conv_mimetic -Grid.File ../grids/test_1p_3d.dgf -Problem.Name $outputName -Grid.Refinement 4 -Grid.UseCheckerboard 1 $inputArgs >> $outputName".out"
