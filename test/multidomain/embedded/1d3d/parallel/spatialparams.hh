// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief Spatial parameters of the tissue and the network of the embedded coupling test on distributed grids
 */
#ifndef DUMUX_TEST_EMBEDDED_PARALLEL_SPATIALPARAMS_HH
#define DUMUX_TEST_EMBEDDED_PARALLEL_SPATIALPARAMS_HH

#include <cmath>
#include <vector>

#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

template<class GridGeometry, class Scalar>
class EmbeddedParallelTissueSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, EmbeddedParallelTissueSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = EmbeddedParallelTissueSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using GlobalPosition = typename GridGeometry::GlobalCoordinate;

public:
    using PermeabilityType = Scalar;

    using ParentType::ParentType;

    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }
};

/*!
 * \brief The network with a radius that changes along the network, so that point sources of
 *        different elements have different values
 */
template<class GridGeometry, class Scalar>
class EmbeddedParallelNetworkSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, EmbeddedParallelNetworkSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = EmbeddedParallelNetworkSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using GlobalPosition = typename GridGeometry::GlobalCoordinate;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;

public:
    using PermeabilityType = Scalar;

    EmbeddedParallelNetworkSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        const auto& gridView = gridGeometry->gridView();
        radius_.resize(gridView.size(0));
        for (const auto& element : elements(gridView))
            radius_[gridGeometry->elementMapper().index(element)] = radiusAtPos(element.geometry().center());
    }

    static Scalar radiusAtPos(const GlobalPosition& globalPos)
    { return 0.02*(1.0 + 0.5*std::sin(7.0*globalPos[2] + 3.0*globalPos[0])); }

    Scalar radius(std::size_t eIdx) const
    { return radius_[eIdx]; }

    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

    template<class ElementSolution>
    Scalar extrusionFactor(const Element& element, const SubControlVolume& scv, const ElementSolution& elemSol) const
    {
        const auto r = radius(this->gridGeometry().elementMapper().index(element));
        return M_PI*r*r;
    }

private:
    std::vector<Scalar> radius_;
};

} // end namespace Dumux

#endif
