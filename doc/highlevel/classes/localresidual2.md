% Introduction to DuMux - Overview and Available Models
% The DuMux Development Team, IWS-LH2, University of Stuttgart
% January 26, 2023

# History and Structure
<!--![DuMux logo](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/raw/master/doc/logo/dumux_logo_hires_whitebg.png)-->

## DuMu$^\mathsf{x}$ is a DUNE module

![](fig/dumux_dune_module.png)

# Available Models

## Models: 1p -- one-phase

* Standard Darcy approach for the conservation of momentum

$$
v = - \frac{\textbf{K}}{\mu}
 \left(\textbf{grad}\, p - \varrho {\textbf{g}} \right)
$$

* Balance equation

$$
\phi \frac{\partial \varrho}{\partial t} + \text{div} \left\lbrace
 - \varrho \frac{\textbf{K}}{\mu} \left( \textbf{grad}\, p -\varrho {\textbf{g}} \right) \right\rbrace = q
$$

* Primary variable: $p$


# Model Components

## Model components: VTK output fields

* Each model provides `...VtkOutputFields`, a class where the default output fields are set by employing a `VtkOutputModule vtk` via, e.g.,
```cpp
   vtk.addVolumeVariable( [](const auto& v)
                          { return v.saturation(FS::phase1Idx); },
                          "S_n" );
```

* In addition to the default fields, custom fields can be added to the output by defining a function `void addFieldsToWriter(VtkWriter& vtk)` in the problem class.

