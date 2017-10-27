#!/bin/sh

make test5_conv_mimetic

inputArgs="-Newton.UseLineSearch 0"

outputName="test5_conv_mimetic_nonconvex"
echo $outputName
rm $outputName"_rates.txt"
rm $outputName".out"
./test5_conv_mimetic -Grid.File ../grids/test_nonconvex_4.dgf -Problem.Name $outputName $inputArgs > $outputName".out"
./test5_conv_mimetic -Grid.File ../grids/test_nonconvex_8.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test5_conv_mimetic -Grid.File ../grids/test_nonconvex_16.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test5_conv_mimetic -Grid.File ../grids/test_nonconvex_32.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"

outputName="test5_conv_mimetic_cossurf"
echo $outputName
rm $outputName"_rates.txt"
rm $outputName".out"
./test5_conv_mimetic -Grid.File ../grids/test_nonconvex_cossurf_4.dgf -Problem.Name $outputName $inputArgs > $outputName".out"
./test5_conv_mimetic -Grid.File ../grids/test_nonconvex_cossurf_8.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test5_conv_mimetic -Grid.File ../grids/test_nonconvex_cossurf_16.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test5_conv_mimetic -Grid.File ../grids/test_nonconvex_cossurf_32.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"

outputName="test5_conv_mimetic_rand"
echo $outputName
rm $outputName"_rates.txt"
rm $outputName".out"
./test5_conv_mimetic -Grid.File ../grids/RandMesh4.msh.dgf -Problem.Name $outputName $inputArgs > $outputName".out"
./test5_conv_mimetic -Grid.File ../grids/RandMesh8.msh.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test5_conv_mimetic -Grid.File ../grids/RandMesh16.msh.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test5_conv_mimetic -Grid.File ../grids/RandMesh32.msh.dgf -Problem.Name $outputName  $inputArgs >> $outputName".out"

outputName="test5_conv_mimetic_checkerboard"
echo $outputName
rm $outputName"_rates.txt"
rm $outputName".out"
./test5_conv_mimetic -Grid.File ../grids/test_1p_3d.dgf -Problem.Name $outputName -Grid.UseCheckerboard 1 $inputArgs > $outputName".out"
./test5_conv_mimetic -Grid.File ../grids/test_1p_3d.dgf -Problem.Name $outputName -Grid.Refinement 1 -Grid.UseCheckerboard 1  $inputArgs >> $outputName".out"
./test5_conv_mimetic -Grid.File ../grids/test_1p_3d.dgf -Problem.Name $outputName -Grid.Refinement 2 -Grid.UseCheckerboard 1 $inputArgs >> $outputName".out"
./test5_conv_mimetic -Grid.File ../grids/test_1p_3d.dgf -Problem.Name $outputName -Grid.Refinement 3 -Grid.UseCheckerboard 1 $inputArgs >> $outputName".out"
