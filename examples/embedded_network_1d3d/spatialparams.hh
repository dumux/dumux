// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_NETWORK_TISSUE_SPATIAL_PARAMS_HH
#define DUMUX_NETWORK_TISSUE_SPATIAL_PARAMS_HH

// # Spatial parameters (`spatialparams.hh`)
//
// This file contains the __Spatial parameter classes__ for each of the
// sub-problems. This header implements two classes: `TissueSpatialParams`
// and `NetworkSpatialParams`.
//
// [[content]]
//
// ### Include headers
//
#include <dune/common/fvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>
//
// # The `TissueSpatialParams` class
//
// The tissue spatial parameter are relatively simple.
// All parameters are constant and are configured via the configuration file.
// This means the parameters can be changed at runtime (after compilation).
//
// [[codeblock]]
namespace Dumux {
template<class FVGridGeometry, class Scalar>
class TissueSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<FVGridGeometry, Scalar, TissueSpatialParams<FVGridGeometry, Scalar>>
{
// [[/codeblock]]
// [[details]] alias definitions
// [[codeblock]]
    using ThisType = TissueSpatialParams<FVGridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<FVGridGeometry, Scalar, ThisType>;
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
// [[/codeblock]]
// [[/details]] alias definitions
    //
    // We read parameter from the configuration file. In this case
    // all parameters are constants and are configurable at runtime
    // through the configuration file.
    //
    // [[codeblock]]
public:
    TissueSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        porosity_ = getParam<Scalar>("Tissue.SpatialParams.Porosity");
        fluidDensity_ = getParam<Scalar>("Fluid.LiquidDensity", 1000);
        fluidMolarMass_ = getParam<Scalar>("Fluid.AverageMolarMass", 0.018);
    }

    // porosity of the tissue
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }
    // [[/codeblock]]
    //
    // The volume flux (interstitial fluid flow) is zero
    // because we assume no relevant advection takes place.
    //
    // [[codeblock]]
    template<class ElementVolumeVariables>
    Scalar volumeFlux(const Element &element,
                      const FVElementGeometry& fvGeometry,
                      const ElementVolumeVariables& elemVolVars,
                      const SubControlVolumeFace& scvf) const
    {
        return 0.0;
    }
    // [[/codeblock]]
    //
    // The following interfaces are needed by the (relatively general)
    // tracer model to compute the fluid state and the variables needed
    // to evaluate the discrete model equation in the local residual.
    // The constant temperature is needed for isothermal problem since
    // some interfaces (generic material interfaces) sometimes require temperature
    // as input. However, this does not mean that this temperature is in
    // fact used.
    //
    // [[details]] fluid density, molar mass, temperature and private variables
    // [[codeblock]]
    // fluid density
    Scalar fluidDensity(const Element& element, const SubControlVolume& scv) const
    { return fluidDensity_; }

    // fluid molar mass
    Scalar fluidMolarMass(const Element& element, const SubControlVolume& scv) const
    { return fluidMolarMass_; }

    // This unused here but a required interface for isothermal problems
    // passed to fluid property interfaces that may in general depend on temperature
    Scalar temperature() const
    { return 273.15 + 37.0; }

private:
    Scalar porosity_;
    Scalar fluidDensity_, fluidMolarMass_;
};
//[[/codeblock]]
//[[/details]]
//
// # The `NetworkSpatialParams` class
//
// The network parameters contain two additional interfaces.
// First, the parameters from the grid are read in `readGridParams`
// called by the problem class that stores an owning point to the
// spatial parameters object. Second, the `extrusionFactor` interface
// is used to implement extrusion of a line to a cylinder by multiplying
// all areas with the cross-sectional area. Since the cross-sectional
// are depends on the radius this is a spatially-varying parameter.
// `extrusionFactor` always exists as a member function of the base class
// `FVPorousMediumFlowSpatialParamsOneP` and defaults to $`1.0`$ (no extrusion).
//
// [[codeblock]]
template<class GridGeometry, class Scalar>
class NetworkSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, NetworkSpatialParams<GridGeometry, Scalar>>
{
    // [[/codeblock]]
    //
    // [[details]] alias definitions
    // [[codeblock]]
    using ThisType = NetworkSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    static constexpr int dimWorld = GridView::dimensionworld;
    // [[/codeblock]]
    // [[/details]]
    //
    // In the constructor, read several parameters from the configuration file
    // [[codeblock]]
public:
    NetworkSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        fluidDensity_ = getParam<Scalar>("Fluid.LiquidDensity", 1000);
        fluidMolarMass_ = getParam<Scalar>("Fluid.AverageMolarMass", 0.018);
        membraneDiffusionFactor_ = getParam<Scalar>("Tracer.MembraneDiffusionFactor", 1e4); // 1/m
    }
    // [[/codeblock]]
    //
    // Some interface functions that can be used, for example, in the problem class, to
    // access spatial parameter. The interface `radius` is required by the 1D-3D coupling
    // manager and is used to create the integration points for the coupling term
    //
    // [[codeblock]]
    // this is used in the coupling so we use the outer radius of the network
    Scalar radius(std::size_t eIdxGlobal) const
    { return outerRadius_[eIdxGlobal];}

    // the vessel radius (as given in the grid data)
    Scalar vesselRadius(std::size_t eIdxGlobal) const
    { return vesselRadius_[eIdxGlobal];}

    // The membrane diffusion factor (1/m) to be multiplied with the free diffusion coefficient
    Scalar membraneFactor(std::size_t) const
    { return membraneDiffusionFactor_; }

    // The blood volume flux (m^3/s)
    Scalar volumeFlux(std::size_t scvfIndex) const
    { return volumeFluxes_[scvfIndex]; }

    // All space is available for flow
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

    // This unused here but a required interface for isothermal problems
    // passed to fluid property interfaces that may in general depend on temperature
    Scalar temperature() const
    { return 273.15 + 37.0; }

    // fluid density
    Scalar fluidDensity(const Element& element, const SubControlVolume& scv) const
    { return fluidDensity_; }

    // fluid molar mass
    Scalar fluidMolarMass(const Element& element, const SubControlVolume& scv) const
    { return fluidMolarMass_; }
    // [[/codeblock]]

    // Cross-sectional area of network (computed via vessel lumen cross-sectional area)
    // [[codeblock]]
    template<class ElementSolution>
    Scalar extrusionFactor(const Element& element,
                           const SubControlVolume& scv,
                           const ElementSolution& elemSol) const
    {
        const auto eIdx = this->gridGeometry().elementMapper().index(element);
        const auto lumenRadius = vesselRadius(eIdx);
        return M_PI*lumenRadius*lumenRadius;
    }
    // [[/codeblock]]
    //
    // The following interface functions are used in the `main.cc` file
    // for VTK output. For visualization in ParaView, it is nice to
    // write out the vessel radius for ParaView's `Tube` filter.
    //
    // [[codeblock]]
    // Get the vessel lumen radii for output
    const std::vector<Scalar>& getVesselRadii() const
    { return vesselRadius_; }

    // Get vessel segment length (for post-processing)
    const std::vector<Scalar>& getVesselSegmentLengths() const
    { return vesselSegmentLength_; }

    // Get the outer radii for output
    const std::vector<Scalar>& getOuterRadii() const
    { return outerRadius_; }
    // [[/codeblock]]

    // Read params from dgf
    // the grid file contains the vessel radius for each element
    // [[codeblock]]
    template<class GridData>
    void readGridParams(const GridData& gridData)
    {
        const auto& gg = this->gridGeometry();
        const auto numElements = gg.gridView().size(0);
        vesselRadius_.resize(numElements);
        // the vessel segment length
        vesselSegmentLength_.resize(numElements);
        // outer radius including layers outside lumen that are considered part of the membrane
        outerRadius_.resize(numElements);

        // The grid view is a leafGridView. Parameters are only set on level 0.
        // Elements have to inherit spatial parameters from their father.
        for (const auto& element : elements(gg.gridView()))
        {
            auto level0element = element;
            while (level0element.hasFather())
                level0element = level0element.father();

            const auto eIdx = gg.elementMapper().index(element);
            vesselRadius_[eIdx] = gridData.parameters(level0element)[0];
            vesselSegmentLength_[eIdx] = element.geometry().volume();
        }

        // outer radius adding endothelium and basement membrane
        for (std::size_t eIdx = 0; eIdx < gg.gridView().size(0); ++eIdx)
            outerRadius_[eIdx] = vesselRadius_[eIdx] + 0.4e-6;
    }
    // [[/codeblock]]
    //
    // The blood volume flux for each (sub)-control-volume face
    // is provided from the outside (from a blood-flow solver).
    // The advection-diffusion model requires a volume flux
    // at each (sub)-control-volume face. In case, for example,
    // a velocity field is given instead, volume fluxes can be
    // computed using the normal vector `scvf.unitOuterNormal()`
    // of the (sub)-control-volume face.
    //
    // [[codeblock]]
    // set the volume fluxes computed by a blood flow model
    void setVolumeFluxes(const std::vector<Scalar>& fluxes)
    { volumeFluxes_ = fluxes; }

    template<class ElementVolumeVariables>
    Scalar volumeFlux(const Element &element,
                      const FVElementGeometry& fvGeometry,
                      const ElementVolumeVariables& elemVolVars,
                      const SubControlVolumeFace& scvf) const
    {
        assert(!volumeFluxes_.empty());
        return volumeFluxes_[scvf.index()];
    }
    // [[/codeblock]]
    //
// [[details]] private member variables
// [[codeblock]]
private:
    Scalar porosity_;
    Scalar fluidDensity_, fluidMolarMass_;
    Scalar membraneDiffusionFactor_;

    std::vector<Scalar> vesselRadius_;
    std::vector<Scalar> vesselSegmentLength_;
    std::vector<Scalar> outerRadius_;
    std::vector<Scalar> volumeFluxes_;
};

} // end namespace Dumux
// [[/codeblock]]
// [[/details]]
// [[/content]]
#endif
