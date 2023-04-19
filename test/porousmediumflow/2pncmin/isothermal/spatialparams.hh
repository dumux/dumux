// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPNCMinTests
 * \brief Spatial parameters for the dissolution problem where water is injected in a for flushing precipitated salt clogging a gas reservoir.
 */

#ifndef DUMUX_INJECTION_SPATIAL_PARAMETERS_HH
#define DUMUX_INJECTION_SPATIAL_PARAMETERS_HH

#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/porosityprecipitation.hh>
#include <dumux/material/fluidmatrixinteractions/permeabilitykozenycarman.hh>

namespace Dumux {

/*!
 * \ingroup TwoPNCMinTests
 * \brief Spatial parameters for the dissolution problem
 * where water is injected in a for flushing precipitated salt clogging a gas reservoir.
 */
template<class GridGeometry, class Scalar>
class DissolutionSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar,
                                           DissolutionSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;

    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar,
                                                         DissolutionSpatialParams<GridGeometry, Scalar>>;

    using PcKrSwCurve = FluidMatrix::BrooksCoreyDefault<Scalar>;

    using GlobalPosition = typename SubControlVolume::GlobalPosition;

public:
    // type used for the permeability (i.e. tensor or scalar)
    using PermeabilityType = Scalar;

    DissolutionSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , pcKrSwCurve_("SpatialParams")
    {
        solubilityLimit_       = getParam<Scalar>("SpatialParams.SolubilityLimit", 0.26);
        referencePorosity_     = getParam<Scalar>("SpatialParams.referencePorosity", 0.11);
        referencePermeability_ = getParam<Scalar>("SpatialParams.referencePermeability", 2.23e-14);
        temperature_           = getParam<Scalar>("Problem.Temperature");
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     *  \param element The finite volume element
     *  \param scv The sub-control volume
     *  \param elemSol The element solution
     */
    template<class ElementSolution>
    Scalar temperature(const Element& element,
                       const SubControlVolume& scv,
                       const ElementSolution& elemSol) const
    { return temperature_; }

    /*!
     *  \brief Defines the minimum porosity \f$[-]\f$ distribution
     *
     *  \param element The finite element
     *  \param scv The sub-control volume
     */
    Scalar minimalPorosity(const Element& element, const SubControlVolume &scv) const
    { return 1e-5; }

    /*!
     *  \brief Defines the volume fraction of the inert component
     *
     *  \param globalPos The global position in the domain
     *  \param compIdx The index of the inert solid component
     */
    template<class SolidSystem>
    Scalar inertVolumeFractionAtPos(const GlobalPosition& globalPos, int compIdx) const
    { return 1.0-referencePorosity_; }

    /*!
     *  \brief Defines the reference porosity \f$[-]\f$ distribution.
     *
     *  This is the porosity of the porous medium without any of the
     *  considered solid phases.
     *
     *  \param element The finite element
     *  \param scv The sub-control volume
     */
    Scalar referencePorosity(const Element& element, const SubControlVolume &scv) const
    { return referencePorosity_; }

    /*! Intrinsic permeability tensor K \f$[m^2]\f$ depending
     *  on the position in the domain
     *
     *  \param element The finite volume element
     *  \param scv The sub-control volume
     *  \param elemSol The element solution
     *
     *  Solution dependent permeability function
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        auto priVars = evalSolution(element, element.geometry(), elemSol, scv.center(), /*ignoreState=*/true);

        Scalar sumPrecipitates = priVars[/*numComp*/3];

        using std::max;
        const auto poro = max(/*minPoro*/1e-5, referencePorosity_ - sumPrecipitates);
        return permLaw_.evaluatePermeability(referencePermeability_, referencePorosity_, poro);
    }

    Scalar solubilityLimit() const
    { return solubilityLimit_; }

    Scalar theta(const SubControlVolume &scv) const
    { return 10.0; }

    /*!
     * \brief Returns the fluid-matrix interaction law at a given location
     * \param globalPos A global coordinate vector
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition &globalPos) const
    {
        return makeFluidMatrixInteraction(pcKrSwCurve_);
    }

    // define which phase is to be considered as the wetting phase
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::H2OIdx; }

private:

    PermeabilityKozenyCarman<PermeabilityType> permLaw_;

    Scalar solubilityLimit_;
    Scalar referencePorosity_;
    PermeabilityType referencePermeability_ = 0.0;
    const PcKrSwCurve pcKrSwCurve_;
    Scalar temperature_;
};

} // end namespace Dumux

#endif
