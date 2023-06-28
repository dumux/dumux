# SpatialParams

Depending on the utilized model, a `SpatialParams` class might be necessary. Here, spatial distributions of material parameters or parameters in general are specified, for example, the temperature distribution in an isothermal simulation setup. For a porous medium model, one would define the distribution of the porosity or permeability here. For a water retention model, one could also set the distribution of the van Genuchten parameters within the `SpatialParams`. All instances of `SpatialParams` inherit from the [base SpatialParams](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/dumux/common/fvspatialparams.hh).
### Key functionalities

* extrusionFactor(element, scv, ...):
    - Return how much the domain is extruded at a given `element` or sub-control volume `scv`. With extrusion e.g. a 2D domain can be given depth via the `extrusionFactor` thus making the 2D domain a pseudo 3D domain. The `extrusionFactor` has a default value of 1.
* extrusionFactorAtPos(pos):
    - Return how much the domain is extruded at a given position `pos`.
* temperature(element, scv, ...):
    - Return the temperature in the given `element` or sub-control volume `scv`.
* temperatureAtPos(pos):
    - Return the temperature at a given position `pos`.
* gravity(pos):
    - Return the acceleration vector due to gravity. Can be customized at given positions `pos` if necessary.
* gridGeometry():
    - Return the finite volume grid geometry belonging to the spatial parameters object.

### Overview

@mermaid{spatialparams}
