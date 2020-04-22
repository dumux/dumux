# Part 1: Implementation of the 2p set-up

The two-phase flow set-up is implemented in the files `properties.hh`,
`problem.hh` and `spatialparams.hh`. In the file `problem.hh` a description of the boundary and initial conditions can be found. Additionally, you can see how to implement an injection well as a point source and how to read the initial solution from a text file. In the file `spatialparams.hh` we create a lens within the porous medium, that has different spatial parameters than the surrounding material and set the parameters for the
$`p_c - S_w`$ relationship.

[[_TOC_]]
