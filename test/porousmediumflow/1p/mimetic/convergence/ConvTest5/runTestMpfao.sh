#!/bin/sh

make test5_conv_mpfao

inputArgs="-Newton.UseLineSearch 0"

outputName="test5_conv_mpfao_nonconvex"
echo $outputName
rm $outputName"_rates.txt"
rm $outputName".out"
./test5_conv_mpfao -Grid.File ../grids/test_nonconvex_4.dgf -Problem.Name $outputName $inputArgs > $outputName".out"
./test5_conv_mpfao -Grid.File ../grids/test_nonconvex_8.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test5_conv_mpfao -Grid.File ../grids/test_nonconvex_16.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test5_conv_mpfao -Grid.File ../grids/test_nonconvex_32.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"

outputName="test5_conv_mpfao_cossurf"
echo $outputName
rm $outputName"_rates.txt"
rm $outputName".out"
./test5_conv_mpfao -Grid.File ../grids/test_nonconvex_cossurf_4.dgf -Problem.Name $outputName $inputArgs > $outputName".out"
./test5_conv_mpfao -Grid.File ../grids/test_nonconvex_cossurf_8.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test5_conv_mpfao -Grid.File ../grids/test_nonconvex_cossurf_16.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test5_conv_mpfao -Grid.File ../grids/test_nonconvex_cossurf_32.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"

outputName="test5_conv_mpfao_rand"
echo $outputName
rm $outputName"_rates.txt"
rm $outputName".out"
./test5_conv_mpfao -Grid.File ../grids/RandMesh4.msh.dgf -Problem.Name $outputName $inputArgs > $outputName".out"
./test5_conv_mpfao -Grid.File ../grids/RandMesh8.msh.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test5_conv_mpfao -Grid.File ../grids/RandMesh16.msh.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test5_conv_mpfao -Grid.File ../grids/RandMesh32.msh.dgf -Problem.Name $outputName  $inputArgs >> $outputName".out"