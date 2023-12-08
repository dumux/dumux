// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Darcy Brinkman model for a single-domain evaluation of coupled freeflow and porous medium flows
 */

#ifndef DUMUX_BRINKMAN_SPATIAL_PARAMS_HH
#define DUMUX_BRINKMAN_SPATIAL_PARAMS_HH

#include <cmath>
#include <numeric>
#include <algorithm>
#include <dune/common/fmatrix.hh>
#include <dumux/freeflow/spatialparams.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief The spatial parameters class for the Darcy-Brinkman model test
 */
template<class GridGeometry, class Scalar>
class BrinkmanTestSpatialParams
: public BrinkmanSpatialParams<GridGeometry, Scalar, BrinkmanTestSpatialParams<GridGeometry, Scalar>>
{
    using ParentType = BrinkmanSpatialParams<GridGeometry, Scalar, BrinkmanTestSpatialParams<GridGeometry, Scalar>>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    using PermeabilityType = DimWorldMatrix;

    BrinkmanTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , permeability_(0.0)
    , inversePermeability_(0.0)
    {
        initPermeability_();

        pmLowerLeft_ = getParam<GlobalPosition>("SpatialParams.PorousMediumLowerLeft");
        pmUpperRight_ = getParam<GlobalPosition>("SpatialParams.PorousMediumUpperRight");
        pmCenter_ = 0.5*(pmLowerLeft_ + pmUpperRight_);
        pmDimensions_ = pmUpperRight_ - pmLowerLeft_;
        pmRoundedness_ = getParam<Scalar>("SpatialParams.PorousMediumRoundedness", 0.0);
        pmTransitionHalfLength_ = 0.5*getParam<Scalar>("SpatialParams.PorousMediumTransitionLength");
    }

    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return permeability_ * brinkmanEpsilonAtPos(globalPos); }

    PermeabilityType inversePermeability(const Element& element, const FVElementGeometry& fvGeometry, const SubControlVolume& scv) const
    { return inversePermeability_; }

    Scalar brinkmanEpsilon(const Element& element, const FVElementGeometry& fvGeometry, const SubControlVolume& scv) const
    { return brinkmanEpsilonAtPos(scv.center()); }

    Scalar brinkmanEpsilonAtPos(const GlobalPosition& globalPos) const
    {
        const auto p = globalPos - pmCenter_;
        GlobalPosition q(0.0);
        for (int i = 0; i < dimWorld; ++i)
            q[i] += std::abs(p[i]) - 0.5*pmDimensions_[i] + pmRoundedness_;

        const auto length = std::hypot(std::max(q[0], 0.0), std::max(q[1], 0.0));
        const auto sdf = std::min(std::max(q[0], q[1]), 0.0) + length - pmRoundedness_;
        const auto clampedSDF = -std::clamp(sdf, -pmTransitionHalfLength_, pmTransitionHalfLength_) + pmTransitionHalfLength_;
        return 0.5*clampedSDF/pmTransitionHalfLength_;
    }

private:
    void initPermeability_()
    {
        Scalar k = getParam<Scalar>("SpatialParams.Permeability");
        Scalar anisotropyRatio = getParam<Scalar>("SpatialParams.AnisotropyRatio", 0.0);
        permeability_[0][0] = k;
        permeability_[1][1] = k * anisotropyRatio;

        // rotate the tensor by theta
        Scalar theta_ = getParam<Scalar>("SpatialParams.PermeabilityRotation", 0.0);

        // degrees to Radians for the rotation angle, store rotation entries
        Scalar radTheta = theta_ * M_PI / 180.0;
        Scalar cosTheta = std::cos(radTheta);
        Scalar sinTheta = std::sin(radTheta);

        // create a rotation matrix according to the rotation angle
        PermeabilityType rotationMatrix;
        rotationMatrix[0][0] = cosTheta;
        rotationMatrix[0][1] = sinTheta * -1.0;
        rotationMatrix[1][0] = sinTheta;
        rotationMatrix[1][1] = cosTheta;

        // rotate the permeability tensor
        PermeabilityType originalPermeability = permeability_;
        permeability_ = rotationMatrix * originalPermeability * getTransposed(rotationMatrix);

        inversePermeability_ = permeability_;
        inversePermeability_.invert();
    }

    PermeabilityType permeability_;
    PermeabilityType inversePermeability_;

    PermeabilityType ffPermeability_;

    GlobalPosition pmLowerLeft_;
    GlobalPosition pmUpperRight_;
    GlobalPosition pmCenter_;
    GlobalPosition pmDimensions_;
    Scalar pmRoundedness_;
    Scalar pmTransitionHalfLength_;

    static constexpr Scalar eps_ = 1e-7;
};

} // end namespace Dumux

#endif
