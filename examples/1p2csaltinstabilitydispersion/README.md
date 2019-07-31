This tutorial addresses a problem of robustness with dispersion effects in instability problems

We provide here two exemplary approaches to include dispersion within the Box context. We have observed that the updated calculation of the dispersion tensor within a Newton iteration can lead to severe problems with numerical robustness when dispersivity values are sufficiently large.
In such a case, it helps to update dispersion tensors only after each time step. 
This is documented in the two examples here.
