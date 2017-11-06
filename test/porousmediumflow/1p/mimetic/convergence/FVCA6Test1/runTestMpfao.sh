#!/bin/sh

make test1_conv_mpfao

inputArgs="-Newton.UseLineSearch 0"

outputName="test1_fvca6_mpfao_nonconvex"
echo $outputName
rm $outputName"_rates.txt"
rm $outputName".out"
./test1_conv_mpfao -Grid.File ../grids/nonconvex_4.dgf -Problem.Name $outputName $inputArgs > $outputName".out"
./test1_conv_mpfao -Grid.File ../grids/nonconvex_8.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test1_conv_mpfao -Grid.File ../grids/nonconvex_16.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test1_conv_mpfao -Grid.File ../grids/nonconvex_32.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test1_conv_mpfao -Grid.File ../grids/nonconvex_64.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"

outputName="test1_fvca6_mpfao_cossurf"
echo $outputName
rm $outputName"_rates.txt"
rm $outputName".out"
./test1_conv_mpfao -Grid.File ../grids/nonconvex_cossurf_4.dgf -Problem.Name $outputName $inputArgs > $outputName".out"
./test1_conv_mpfao -Grid.File ../grids/nonconvex_cossurf_8.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test1_conv_mpfao -Grid.File ../grids/nonconvex_cossurf_16.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test1_conv_mpfao -Grid.File ../grids/nonconvex_cossurf_32.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test1_conv_mpfao -Grid.File ../grids/nonconvex_cossurf_64.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"

outputName="test1_fvca6_mpfao_rand"
echo $outputName
rm $outputName"_rates.txt"
rm $outputName".out"
./test1_conv_mpfao -Grid.File ../grids/random_4.dgf -Problem.Name $outputName $inputArgs > $outputName".out"
./test1_conv_mpfao -Grid.File ../grids/random_8.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test1_conv_mpfao -Grid.File ../grids/random_16.dgf -Problem.Name $outputName $inputArgs >> $outputName".out"
./test1_conv_mpfao -Grid.File ../grids/random_32.dgf -Problem.Name $outputName  $inputArgs >> $outputName".out"
./test1_conv_mpfao -Grid.File ../grids/random_64.dgf -Problem.Name $outputName  $inputArgs >> $outputName".out"
