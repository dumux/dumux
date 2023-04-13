// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup RichardsTests
 * \brief Spatial parameters for the RichardsLensProblem.
 */

#ifndef DUMUX_RICHARDS_LENS_SPATIAL_PARAMETERS_HH
#define DUMUX_RICHARDS_LENS_SPATIAL_PARAMETERS_HH

#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>

#include <dumux/porousmediumflow/richards/model.hh>

namespace Dumux {

/*!
 * \ingroup RichardsTests
 * \brief The spatial parameters for the RichardsLensProblem.
 */
template<class GridGeometry, class Scalar>
class RichardsLensSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, RichardsLensSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = RichardsLensSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    enum { dimWorld = GridView::dimensionworld };

    using PcKrSwCurve = FluidMatrix::VanGenuchtenDefault<Scalar>;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    RichardsLensSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , pcKrSwCurveLens_("SpatialParams.Lens")
    , pcKrSwCurveOuterDomain_("SpatialParams.OuterDomain")
    {
        lensLowerLeft_ = {1.0, 2.0};
        lensUpperRight_ = {4.0, 3.0};
        lensK_ = 1e-12;
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
     * \brief Returns the porosity [] at a given location
     *
     * \param globalPos The global position where we evaluate
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.4; }

    /*!
     * \brief Returns the fluid-matrix interaction law for the sub-control volume
     *
     * This method is not actually required by the Richards model, but provided
     * for the convenience of the RichardsLensProblem
     *
     * \param element The current finite element
     * \param scv The sub-control volume
     * \param elemSol The current element solution
     */
    template<class ElementSolution>
    auto fluidMatrixInteraction(const Element& element,
                                const SubControlVolume& scv,
                                const ElementSolution& elemSol) const
    {
        const auto& globalPos = scv.dofPosition();
        return fluidMatrixInteractionAtPos(globalPos);
    }

    /*!
     * \brief Returns the fluid-matrix interaction law at a given location
     * \param globalPos The global coordinates for the given location
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    {
        if (isInLens_(globalPos))
            return makeFluidMatrixInteraction(pcKrSwCurveLens_);
        return makeFluidMatrixInteraction(pcKrSwCurveOuterDomain_);
    }

    /*!
     * \brief Returns the temperature [K] at a given location
     * \param globalPos The global coordinates for the given location
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    { return 273.15 + 10.0; }; // -> 10°C

private:
    bool isInLens_(const GlobalPosition &globalPos) const
    {
        for (int i = 0; i < dimWorld; ++i)
            if (globalPos[i] < lensLowerLeft_[i] - eps_ || globalPos[i] > lensUpperRight_[i] + eps_)
                return false;

        return true;
    }

    static constexpr Scalar eps_ = 1e-6;

    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;

    Scalar lensK_;
    Scalar outerK_;
    const PcKrSwCurve pcKrSwCurveLens_;
    const PcKrSwCurve pcKrSwCurveOuterDomain_;
};

} // end namespace Dumux

#endif
