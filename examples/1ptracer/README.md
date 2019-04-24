This tutorial was copied from /home/felix/Dokumente/Aktuelles_U/19_03_27_DumuxTag/dumux/dumux/test/porousmediumflow/tracer/1ptracer.

# One-phase flow with with random permeability distribution (and a tracer model)

## Problem set-up
This example combines a stationary One-phase flow in a porous medium with randomly generated permeabilities with a tracer model. In the first step, the groundwater-velocity is evaluated under stationary conditions. Based on the velocity field the tracer model is running instationary.

## Random permeability generation
### spatialparams_1p.hh
The follwing code can be found in lines 64-72.

```C++
// generate random fields
        std::mt19937 rand(0);
        std::lognormal_distribution<Scalar> K(std::log(permeability_), std::log(permeability_)*0.1);
        std::lognormal_distribution<Scalar> KLens(std::log(permeabilityLens_), std::log(permeabilityLens_)*0.1);
        for (const auto& element : elements(fvGridGeometry->gridView()))
        {
            const auto eIdx = fvGridGeometry->elementMapper().index(element);
            const auto globalPos = element.geometry().center();
            K_[eIdx] = isInLens_(globalPos) ? KLens(rand) : K(rand);
        }
```
Two lognormal distributions are defined: K and KLens with different means and standard deviations. In this case the mean is defined as the log of the permeability given from the input file and the standard deviation is 10% of the mean.
Within a loop through all elements the randomly generated permeabilities are assigned according to their position in the domain (inside or outside the lense).

## Tracer model
...