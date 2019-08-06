// ## The spatialparams for the tracer problem

// the header guard
#ifndef DUMUX_TRACER_TEST_SPATIAL_PARAMS_HH
#define DUMUX_TRACER_TEST_SPATIAL_PARAMS_HH

// In the file properties.hh all properties are declared.
#include <dumux/porousmediumflow/properties.hh>
// As in the 1p spatialparams we inherit from the spatial parameters for single-phase, finite volumes, which we include here.
#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux {

// In the TracerTestSpatialParams class, we define all functions needed to describe spatially dependent parameters for the tracer_problem.

template<class FVGridGeometry, class Scalar>
class TracerTestSpatialParams
: public FVSpatialParamsOneP<FVGridGeometry, Scalar,
                             TracerTestSpatialParams<FVGridGeometry, Scalar>>
{
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParamsOneP<FVGridGeometry, Scalar,
                                           TracerTestSpatialParams<FVGridGeometry, Scalar>>;

    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Dune::FieldVector<Scalar, dimWorld>;

public:

    TracerTestSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry) {}

    // ### Properties of the porous matrix
    // We define the same porosity for the whole domain as in the 1p spatialparams.
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.2; }

    // We don't consider dispersivity for the tracer transport. So we set the dispersivity coefficient to zero.
    template<class ElementSolution>
    Scalar dispersivity(const Element &element,
                        const SubControlVolume& scv,
                        const ElementSolution& elemSol) const
    { return 0; }

    // ### Properties of the fluid system
    // In the following we define fluid properties that are spatial parameters in the tracer model. They can possible vary with space but are usually constants. Further spatially constant values of the fluid system are defined in the TracerFluidSystem class in problem.hh.

    // We define the fluid density to a constant value of 1000 $`\frac{kg}{m^3}`$.
    Scalar fluidDensity(const Element &element,
                        const SubControlVolume& scv) const
    { return 1000; }

    // We define the fluid molar mass.
    Scalar fluidMolarMass(const Element &element,
                          const SubControlVolume& scv) const
    { return 18.0; }

    Scalar fluidMolarMass(const GlobalPosition &globalPos) const
    { return 18.0; }

    // ### The volume fluxes
    // We define a function which returns the field of volume fluxes. This is e.g. used to calculate the transport of the tracer.
    template<class ElementVolumeVariables>
    Scalar volumeFlux(const Element &element,
                      const FVElementGeometry& fvGeometry,
                      const ElementVolumeVariables& elemVolVars,
                      const SubControlVolumeFace& scvf) const
    {
        return volumeFlux_[scvf.index()];
    }

    // We define a function to set the volume flux. This is used in the main function to set the volume flux to the calculated value based on the solution of the 1p problem.
    void setVolumeFlux(const std::vector<Scalar>& f)
    { volumeFlux_ = f; }

private:
    std::vector<Scalar> volumeFlux_;
};

} // end namespace Dumux

#endif
