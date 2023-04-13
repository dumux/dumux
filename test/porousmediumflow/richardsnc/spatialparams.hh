// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup RichardsNCTests
 * \brief The spatial parameters for the RichardsWellTracerProblem.
 */

#ifndef DUMUX_RICHARDS_LENS_SPATIAL_PARAMETERS_HH
#define DUMUX_RICHARDS_LENS_SPATIAL_PARAMETERS_HH

#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>
#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/porousmediumflow/richards/model.hh>

namespace Dumux {

/*!
 * \ingroup RichardsModel
 * \ingroup ImplicitTestProblems
 * \brief The spatial parameters for the RichardsWellTracerProblem.
 */
template<class GridGeometry, class Scalar>
class RichardsWellTracerSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar,
                                           RichardsWellTracerSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar,
                                                         RichardsWellTracerSpatialParams<GridGeometry, Scalar>>;

    enum { dimWorld=GridView::dimensionworld };
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PcKrSwCurve = FluidMatrix::VanGenuchtenNoReg<Scalar>;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    RichardsWellTracerSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , lensPcKrSwCurve_("SpatialParams.Lens")
    , outerPcKrSwCurve_("SpatialParams.OuterDomain")
    {
        lensLowerLeft_ = getParam<GlobalPosition>("Problem.LensLowerLeft");
        lensUpperRight_ = getParam<GlobalPosition>("Problem.LensUpperRight");

        lensK_ = 1e-14;
        outerK_ = 5e-12;
    }

    /*!
     * \brief Returns the intrinsic permeability tensor [m^2] at a given location
     *
     * \param globalPos The global position where we evaluate
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    {
        if (isInLens_(globalPos))
            return lensK_;
        return outerK_;
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
    { return 273.15 + 10; }; // -> 10°C

    /*!
     * \brief Returns the porosity [] at a given location
     *
     * \param globalPos The global position where we evaluate
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        if (isInLens_(globalPos))
            return 0.2;
        return 0.4;
    }

    /*!
     * \brief Returns the fluid-matrix interaction law at a given location
     *
     * \param globalPos The global coordinates for the given location
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    {
        if (isInLens_(globalPos))
            return makeFluidMatrixInteraction(lensPcKrSwCurve_);
        else
            return makeFluidMatrixInteraction(outerPcKrSwCurve_);
    }

private:
    bool isInLens_(const GlobalPosition &globalPos) const
    {
        for (int i = 0; i < dimWorld; ++i)
            if (globalPos[i] < lensLowerLeft_[i] - eps_ || globalPos[i] > lensUpperRight_[i] + eps_)
                return false;

        return true;
    }

    static constexpr Scalar eps_ = 1.5e-7;

    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;

    Scalar lensK_;
    Scalar outerK_;

    const PcKrSwCurve lensPcKrSwCurve_;
    const PcKrSwCurve outerPcKrSwCurve_;
};

} // end namespace Dumux

#endif
