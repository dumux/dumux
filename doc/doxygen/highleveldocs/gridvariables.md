@page gridvariables GridVariables

gridVariables is a concept which includes all variables needed to solve a PDE. Such a variable could be for exampe the density of your fluid of interest.
In your main file you need to create a gridVariables object. Afterwards you have to call the init() function.
After those steps, you will, dependent on your problem, you will create different objects of classes.
For example if your discretization has variables that are defined inside your volume, a object of type @ref gridVolumeVariables will be created. Is you have a instationary problem (i.e., you need to discretize a time derivate) there will be current and previous object respectivly. Another object that might be created is a gridFluxVariablesCache or gridFaceVariables object. However, this depends on you chosen discretization.

@subpage gridvolumevariables