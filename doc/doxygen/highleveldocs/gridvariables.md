## GridVariables {#gridvariables}
<!-- @page gridVariables GridVariables -->

gridVariables is a concept which includes all variables needed to solve a PDE. Such a variable could be for exampe the density of your fluid of interest.
Those variables are needed to compute a system of equations.
Dumux has discretizations that have degrees of freedoms inside a volume and some have degrees of freedom on the face.
Hence, when creating gridVariables, @ref gridVolumeVariables and @ref gridFaceVariables are created.
Since Dumux uses a element-wise assembly for the system of equations, depending on your discretization @ref elementVolumeVariables and @ref elementFaceVariables are created.
