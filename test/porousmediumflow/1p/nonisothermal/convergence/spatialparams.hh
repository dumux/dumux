// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Spatial parameters of the convergence test of the non-isothermal single-phase model
 */
#ifndef DUMUX_ONEPNI_CONVERGENCE_SPATIALPARAMS_HH
#define DUMUX_ONEPNI_CONVERGENCE_SPATIALPARAMS_HH

#include <dumux/common/parameters.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief Spatial parameters of the convergence test of the non-isothermal single-phase model
 *
 * The porosity \f$ \phi = \phi_0 (1 + x y) \f$ varies in space, so that the effective thermal
 * conductivity \f$ \phi \lambda_f + (1 - \phi) \lambda_s \f$ does as well.
 */
template<class GridGeometry, class Scalar>
class OnePNIConvergenceSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                                             OnePNIConvergenceSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = OnePNIConvergenceSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using GlobalPosition = typename GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;

public:
    using PermeabilityType = Scalar;

    OnePNIConvergenceSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , permeability_(getParam<Scalar>("SpatialParams.Permeability"))
    , referencePorosity_(getParam<Scalar>("SpatialParams.ReferencePorosity"))
    {}

    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return permeability_; }

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return referencePorosity_*(1.0 + globalPos[0]*globalPos[1]); }

    //! The analytical fields are evaluated at the position of the local dof
    template<class ElementDiscretization, class IpData, class ElementSolution>
    PermeabilityType permeability(const ElementDiscretization&, const IpData& ipData, const ElementSolution&) const
    { return permeabilityAtPos(ipData.global()); }

    //! The analytical fields are evaluated at the position of the local dof
    template<class ElementDiscretization, class IpData, class ElementSolution>
    Scalar porosity(const ElementDiscretization&, const IpData& ipData, const ElementSolution&) const
    { return porosityAtPos(ipData.global()); }

    using ParentType::permeability;
    using ParentType::porosity;

private:
    Scalar permeability_;
    Scalar referencePorosity_;
};

} // end namespace Dumux

#endif
