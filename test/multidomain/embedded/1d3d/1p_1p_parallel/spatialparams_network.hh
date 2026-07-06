// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief Spatial parameters for a 1d network with a per-element radius supplied as a vector.
 *
 * Unlike the analytical blood-flow spatial params (constant radius from a parameter), the radius
 * here is given per local element. This is what a DGF-loaded network needs in parallel: after the
 * network is distributed by the bulk partition, the per-element radius read from the DGF on the
 * full grid is migrated onto the local elements (by global id) and injected here.
 */
#ifndef DUMUX_NETWORK_SPATIALPARAMS_HH
#define DUMUX_NETWORK_SPATIALPARAMS_HH

#include <vector>
#include <memory>
#include <cmath>

#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup EmbeddedTests
 * \brief Spatial parameters for a 1d network with an injected per-element radius vector.
 */
template<class GridGeometry, class Scalar>
class NetworkSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, NetworkSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = NetworkSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using GlobalPosition = typename GridGeometry::GlobalCoordinate;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;

public:
    using PermeabilityType = Scalar;

    //! \param radii the per-(local-)element radius, indexed by the leaf element index
    NetworkSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                         std::vector<Scalar> radii)
    : ParentType(gridGeometry)
    , radii_(std::move(radii))
    {}

    //! constant intrinsic permeability of the xylem (well-posed; value irrelevant to the test)
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

    //! the radius of the circular pipe for element eIdxGlobal in [m]
    Scalar radius(std::size_t eIdxGlobal) const
    { return radii_[eIdxGlobal]; }

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

    //! extrude the 1d line to a circular tube of cross-section area pi*r^2
    template<class ElementSolution>
    Scalar extrusionFactor(const Element& element,
                           const SubControlVolume& scv,
                           const ElementSolution& elemSol) const
    {
        const auto eIdx = this->gridGeometry().elementMapper().index(element);
        const auto r = radius(eIdx);
        return M_PI*r*r;
    }

    //! access to the radii (e.g. for vtk output)
    const std::vector<Scalar>& getRadii() const
    { return radii_; }

private:
    std::vector<Scalar> radii_;
};

} // end namespace Dumux

#endif
