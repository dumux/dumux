// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_COOKS_MEMBRANE_PRESSURE_MIXED_PROBLEM_HH
#define DUMUX_COOKS_MEMBRANE_PRESSURE_MIXED_PROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/problemwithspatialparams.hh>
#include <dumux/common/properties.hh>

namespace Dumux {

/*!
 * \brief Pressure subproblem of Cook's membrane with the mixed u–p formulation.
 *
 * The local constraint \f$ p_s = \lambda (J^2 - 1)/(2J) \f$, the derivative of the volumetric energy
 * \f$ \frac{\lambda}{4} (J^2 - 1 - 2 \ln J) \f$, determines \f$ p_s \f$ without boundary conditions.
 */
template<class TypeTag>
class CooksMembranePressureMixedProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::Experimental::BoundaryTypes<PrimaryVariables::size()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    CooksMembranePressureMixedProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                                      std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    {}

    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    Scalar volumetricPressure(Scalar J) const
    { return this->spatialParams().firstLameParameter() * (J*J - 1.0) / (2.0*J); }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition&) const
    { return BoundaryTypes{}; }

    PrimaryVariables initialAtPos(const GlobalPosition&) const
    { return PrimaryVariables(0.0); }

private:
    std::shared_ptr<CouplingManager> couplingManager_;
};

} // end namespace Dumux

#endif
