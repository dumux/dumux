// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoundaryTests
 * \brief The spatial parameters for the incompressible test.
 */

#ifndef DUMUX_INCOMPRESSIBLE_ONEP_TEST_SPATIAL_PARAMS_HH
#define DUMUX_INCOMPRESSIBLE_ONEP_TEST_SPATIAL_PARAMS_HH

#include <limits>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {
namespace LensSpatialParams {
/*!
 * \brief Returns if a point is in a lens with a given bounding box
 *
 * \param globalPos The position of the point
 * \param lensLowerLeft The lower left corner of the lens
 * \param lensUpperRight The upper right corner of the lens
 */
template<class GlobalPosition>
bool pointInLens(const GlobalPosition& globalPos,
                 const GlobalPosition& lensLowerLeft,
                 const GlobalPosition& lensUpperRight)
{
    const auto eps = 1e-8*(lensUpperRight - lensLowerLeft).two_norm();
    for (int i = 0; i < GlobalPosition::size(); ++i)
        if (globalPos[i] < lensLowerLeft[i] + eps || globalPos[i] > lensUpperRight[i] - eps)
            return false;

    return true;
}
} // end namespace LensSpatialParams

/*!
 * \ingroup BoundaryTests
 * \brief The spatial parameters class for the test problem using the
 *        incompressible 1p model.
 */
template<class GridGeometry, class Scalar>
class OnePTestSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, OnePTestSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = OnePTestSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using PermeabilityType = Scalar;
    OnePTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , lensLowerLeft_(std::numeric_limits<Scalar>::max())
    , lensUpperRight_(std::numeric_limits<Scalar>::lowest())
    {
        permeability_ = getParam<Scalar>("SpatialParams.Permeability");
        permeabilityLens_ = getParam<Scalar>("SpatialParams.PermeabilityLens");
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
     *
     * \return The intrinsic permeability
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return isInLens_(globalPos) ? permeabilityLens_ : permeability_; }

    /*!
     * \brief Defines the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.4; }

    /*!
     * \brief Optionally set a lens
     */
    void setLens(const GlobalPosition& lowerLeft, const GlobalPosition& upperRight)
    { lensLowerLeft_ = lowerLeft; lensUpperRight_ = upperRight; }

private:
    bool isInLens_(const GlobalPosition &globalPos) const
    { return LensSpatialParams::pointInLens(globalPos, lensLowerLeft_, lensUpperRight_); }

    GlobalPosition lensLowerLeft_, lensUpperRight_;
    Scalar permeability_, permeabilityLens_;
};

} // end namespace Dumux

#endif
