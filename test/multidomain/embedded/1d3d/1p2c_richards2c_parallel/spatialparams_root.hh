// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief Root xylem spatial parameters with the per-element radius supplied as a vector.
 *
 * Identical physics to the sequential RootSpatialParams (constant Kx/Kr, radius drives
 * permeability and extrusion), but the radius is injected per local element. This is what the
 * distributed network needs: the radius is read from the DGF on the full grid and migrated onto
 * the local elements (by global id) before being injected here.
 */
#ifndef DUMUX_ROOT_SPATIALPARAMS_PARALLEL_HH
#define DUMUX_ROOT_SPATIALPARAMS_PARALLEL_HH

#include <vector>
#include <memory>
#include <cmath>

#include <dumux/common/parameters.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup EmbeddedTests
 * \brief Root xylem spatial parameters with an injected per-element radius vector.
 */
template<class GridGeometry, class Scalar>
class RootSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, RootSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = RootSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using PermeabilityType = Scalar;

    //! \param radii the per-(local-)element radius, indexed by the leaf element index
    RootSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                      std::vector<Scalar> radii)
    : ParentType(gridGeometry)
    , radii_(std::move(radii))
    {
        porosity_ = getParam<Scalar>("Root.SpatialParams.Porosity", 0.4);
        constantKx_ = getParam<Scalar>("Root.SpatialParams.Kx", 5.0968e-17);
        constantKr_ = getParam<Scalar>("Root.SpatialParams.Kr", 2.04e-13);
    }

    //! axial conductivity translated into an intrinsic permeability
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        const Scalar r = radius(this->gridGeometry().gridView().indexSet().index(element));
        return constantKx_ / (M_PI*r*r) * Components::SimpleH2O<Scalar>::liquidViscosity(285.15, 1e5);
    }

    //! the radius of the root segment eIdxGlobal in [m]
    Scalar radius(std::size_t eIdxGlobal) const
    { return radii_[eIdxGlobal]; }

    //! the radial conductivity
    Scalar Kr(std::size_t eIdxGlobal) const
    { return constantKr_; }

    const std::vector<Scalar>& getRadii() const
    { return radii_; }

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }

    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    { return 273.15 + 10.0; }

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

private:
    std::vector<Scalar> radii_;
    Scalar porosity_, constantKx_, constantKr_;
};

} // end namespace Dumux

#endif
