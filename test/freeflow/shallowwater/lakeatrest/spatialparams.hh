// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief The spatial parameters for the lake-at-rest problem.
 */
#ifndef DUMUX_LAKE_AT_REST_SPATIAL_PARAMETERS_HH
#define DUMUX_LAKE_AT_REST_SPATIAL_PARAMETERS_HH

#include <cmath>

#include <dumux/common/parameters.hh>
#include <dumux/freeflow/spatialparams.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/frictionlaw.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/nofriction.hh>

namespace Dumux {

/*!
 * \ingroup ShallowWaterTests
 * \brief The spatial parameters class for the lake-at-rest test.
 *
 * The bed consists of a smooth Gaussian hump superposed with a rectangular sill
 * with discontinuous flanks
 * \f[
 *   z(x,y) = a \exp\left(-5(x-0.9)^2 - 50(y-0.5)^2\right)
 *          + \begin{cases} b & 1.4 < x < 1.6 \\ 0 & \text{else} \end{cases}
 * \f]
 * on the domain \f$ [0,2] \times [0,1] \f$. The hump is the bed of the two-dimensional
 * still water test of Randall LeVeque, "Balancing source terms and flux gradients in
 * high-resolution Godunov methods: the quasi-steady wave-propagation algorithm", Journal
 * of Computational Physics, 146(1):346-365, 1998,
 * doi: https://doi.org/10.1006/jcph.1998.6058.
 * The sill adds a bed discontinuity, across which the hydrostatic reconstruction has to
 * be well-balanced as well.
 */
template<class GridGeometry, class Scalar, class VolumeVariables>
class LakeAtRestSpatialParams
: public FreeFlowSpatialParams<GridGeometry, Scalar, LakeAtRestSpatialParams<GridGeometry, Scalar, VolumeVariables>>
{
    using ThisType = LakeAtRestSpatialParams<GridGeometry, Scalar, VolumeVariables>;
    using ParentType = FreeFlowSpatialParams<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    LakeAtRestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , gravity_(getParam<Scalar>("Problem.Gravity"))
    , humpHeight_(getParam<Scalar>("Problem.HumpHeight"))
    , sillHeight_(getParam<Scalar>("Problem.SillHeight"))
    , frictionLaw_(std::make_unique<FrictionLawNoFriction<VolumeVariables>>())
    {}

    //! Define the gravitation
    Scalar gravity(const GlobalPosition& globalPos) const
    { return gravity_; }

    //! Get the friction law, which already includes the friction value
    const FrictionLaw<VolumeVariables>& frictionLaw(const Element& element,
                                                    const SubControlVolume& scv) const
    { return *frictionLaw_; }

    //! Define the bed surface
    Scalar bedSurface(const Element& element,
                      const SubControlVolume& scv) const
    { return bedSurfaceAtPos(element.geometry().center()); }

    //! The bed surface as a function of the position
    Scalar bedSurfaceAtPos(const GlobalPosition& globalPos) const
    {
        using std::exp;
        const auto x = globalPos[0];
        const auto y = globalPos[1];

        auto bedSurface = humpHeight_*exp(-5.0*(x - 0.9)*(x - 0.9) - 50.0*(y - 0.5)*(y - 0.5));
        if (x > 1.4 && x < 1.6)
            bedSurface += sillHeight_;

        return bedSurface;
    }

private:
    Scalar gravity_;
    Scalar humpHeight_;
    Scalar sillHeight_;
    std::unique_ptr<FrictionLaw<VolumeVariables>> frictionLaw_;
};

} // end namespace Dumux

#endif
