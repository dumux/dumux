// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FacetTests
 * \brief The spatial parameters for the single-phase facet coupling test.
 */

#ifndef DUMUX_TEST_TPFAFACETCOUPLING_THREEDOMAIN_ONEP_SPATIALPARAMS_HH
#define DUMUX_TEST_TPFAFACETCOUPLING_THREEDOMAIN_ONEP_SPATIALPARAMS_HH

#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup FacetTests
 * \brief The spatial parameters for the single-phase facet coupling test.
 */
template<class GridGeometry, class Scalar>
class OnePSpatialParams : public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                                                                 OnePSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using ThisType = OnePSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;

public:
    //! Export the type used for permeability
    using PermeabilityType = Scalar;

    OnePSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry, const std::string& paramGroup = "")
    : ParentType(gridGeometry)
    {
        permeability_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Permeability");
        extrusion_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture", 1.0);
    }

    //! Function for defining the (intrinsic) permeability \f$[m^2]\f$.
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return permeability_; }

    //! Returns the porosity.
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

    //! Returns the extrusion factor.
    Scalar extrusionFactorAtPos(const GlobalPosition& globalPos) const
    { return extrusion_; }

private:
    Scalar permeability_;
    Scalar extrusion_;
};

} // end namespace Dumux

#endif
