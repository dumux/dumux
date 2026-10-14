// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_COOKS_MEMBRANE_PRESSURE_MIXED_PROBLEM_HH
#define DUMUX_COOKS_MEMBRANE_PRESSURE_MIXED_PROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

namespace Dumux {

// Pressure subdomain problem for the MINI mixed formulation.
// The only equation is the local bulk pressure constraint — no BCs needed
// (the constraint (J²-1)/(2J) = p_s/λ uniquely determines p_s everywhere).
template<class TypeTag>
class CooksMembranePressureMixedProblem : public FVProblemWithSpatialParams<TypeTag>
{
    using ParentType = FVProblemWithSpatialParams<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    CooksMembranePressureMixedProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                                       std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    {}

    const CouplingManager& couplingManager() const { return *couplingManager_; }

    // Volumetric pressure for ψ₁ = µ/2·(I_C−3) + λ/4·(J²−1) − (λ/2+µ)·lnJ: p_s^eq(J) = λ(J²−1)/(2J)
    Scalar volumetricPressure(Scalar J) const
    { return this->spatialParams().firstLameParameter() * (J*J - 1.0) / (2.0*J); }

    // pressure constraint is purely algebraic: no spatial flux
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition&) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    // zero boundary flux (pressure constraint has no spatial boundary flux)
    template<class ElementVariables>
    NumEqVector boundaryFluxIntegral(const auto&, const ElementVariables&, const auto&) const
    { return NumEqVector(0.0); }

    PrimaryVariables initialAtPos(const GlobalPosition&) const
    { return PrimaryVariables(0.0); }

private:
    std::shared_ptr<CouplingManager> couplingManager_;
};

} // end namespace Dumux
#endif
