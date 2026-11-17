// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief The spatial parameters of the tissue and vessel problems stated with the variables of the local dofs
 */
#ifndef DUMUX_EMBEDDED_1P1P_VARIABLES_SPATIALPARAMS_HH
#define DUMUX_EMBEDDED_1P1P_VARIABLES_SPATIALPARAMS_HH

#include <cmath>
#include <memory>
#include <type_traits>

#include <dumux/common/parameters.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>
#include <dumux/porousmediumflow/spatialparams.hh>

namespace Dumux {

template<class GridGeometry, class Scalar>
class TissueSpatialParams;

namespace Detail {
//! The base class with sub-control volumes for the schemes that have them, the general one otherwise
template<class GridGeometry, class Scalar>
using TissueSpatialParamsBase = std::conditional_t<
    requires { typename GridGeometry::SubControlVolume; },
    FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, TissueSpatialParams<GridGeometry, Scalar>>,
    Experimental::PorousMediumFlowSpatialParams<GridGeometry, Scalar, TissueSpatialParams<GridGeometry, Scalar>>
>;
} // end namespace Detail

//! Constant permeability and porosity of the tissue
template<class GridGeometry, class Scalar>
class TissueSpatialParams
: public Detail::TissueSpatialParamsBase<GridGeometry, Scalar>
{
    using ParentType = Detail::TissueSpatialParamsBase<GridGeometry, Scalar>;
    using GlobalPosition = typename GridGeometry::GlobalCoordinate;
public:
    using PermeabilityType = Scalar;

    TissueSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , permeability_(getParam<Scalar>("SpatialParams.PermeabilityTissue"))
    {}

    PermeabilityType permeabilityAtPos(const GlobalPosition&) const
    { return permeability_; }

    Scalar porosityAtPos(const GlobalPosition&) const
    { return 1.0; }

    using ParentType::permeability;
    using ParentType::porosity;

    template<class ElementDiscretization, class IpData, class ElementSolution>
    PermeabilityType permeability(const ElementDiscretization&, const IpData& ipData, const ElementSolution&) const
    { return permeabilityAtPos(ipData.global()); }

    template<class ElementDiscretization, class IpData, class ElementSolution>
    Scalar porosity(const ElementDiscretization&, const IpData& ipData, const ElementSolution&) const
    { return porosityAtPos(ipData.global()); }

private:
    Scalar permeability_;
};

//! Permeability varying along the vessel and extrusion to a tube of constant radius
template<class GridGeometry, class Scalar>
class BloodFlowSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, BloodFlowSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = BloodFlowSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using GlobalPosition = typename GridGeometry::GlobalCoordinate;
public:
    using PermeabilityType = Scalar;

    BloodFlowSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , radius_(getParam<Scalar>("SpatialParams.Radius"))
    , permeabilityVariation_(getParam<Scalar>("SpatialParams.PermeabilityVariation", 1.0))
    {}

    PermeabilityType permeabilityAtPos(const GlobalPosition& pos) const
    { return (1 + permeabilityVariation_*(pos[2] + 0.5*pos[2]*pos[2]))/(M_PI*radius_*radius_); }

    Scalar porosityAtPos(const GlobalPosition&) const
    { return 1.0; }

    Scalar extrusionFactorAtPos(const GlobalPosition&) const
    { return M_PI*radius_*radius_; }

    //! The radius of the vessel in the element with the index eIdx
    Scalar radius(unsigned int eIdx) const
    { return radius_; }

    using ParentType::permeability;
    using ParentType::porosity;
    using ParentType::extrusionFactor;

    template<class ElementDiscretization, class IpData, class ElementSolution>
    PermeabilityType permeability(const ElementDiscretization&, const IpData& ipData, const ElementSolution&) const
    { return permeabilityAtPos(ipData.global()); }

    template<class ElementDiscretization, class IpData, class ElementSolution>
    Scalar porosity(const ElementDiscretization&, const IpData& ipData, const ElementSolution&) const
    { return porosityAtPos(ipData.global()); }

private:
    Scalar radius_;
    Scalar permeabilityVariation_;
};

} // end namespace Dumux

#endif
