// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief One-sided Rayleigh-Bénard dissolution problem (pressure-based Boussinesq).
 *
 * Primary variables (OnePNC, 2 components, replaceCompEqIdx = 0):
 *   - p   (pressure, index 0): driven by total mass balance ∇·u = 0
 *   - C   (solute mass fraction, index 1): solute transport
 *
 * Boundary conditions:
 *   - C = 1 Dirichlet at top
 *   - ∂C/∂n = 0 elsewhere
 *   - Pressure: Neumann (no-flow) everywhere
 *   - Pressure reference fixed weakly at top via Robin BC:
 *       neumann[0] = (p_inside - p_ref) * penalty
 */
#ifndef DUMUX_BOUSSINESQ_PRESSURE_PROBLEM_HH
#define DUMUX_BOUSSINESQ_PRESSURE_PROBLEM_HH

#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

template<class TypeTag>
class BoussinesqPressureProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType   = PorousMediumFlowProblem<TypeTag>;
    using Scalar       = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView     = typename GridGeometry::GridView;
    using FluidSystem  = GetPropType<TypeTag, Properties::FluidSystem>;
    using ModelTraits  = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices      = typename ModelTraits::Indices;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;

    using BoundaryTypes    = Dumux::BoundaryTypes<ModelTraits::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector      = Dumux::NumEqVector<PrimaryVariables>;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ElementVolumeVariables =
        typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache =
        typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int dimWorld    = GridView::dimensionworld;
    static constexpr int pressureIdx = Indices::pressureIdx;
    static constexpr int conti0EqIdx = Indices::conti0EqIdx;
    // solute mass fraction is primary variable 1 (component 1 in OnePNC with 2 components)
    static constexpr int soluteIdx   = FluidSystem::soluteIdx;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    BoussinesqPressureProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        pRef_    = getParam<Scalar>("Problem.PRef", 0.0);
        penalty_ = getParam<Scalar>("Problem.PenaltyFactor", 1e6);
    }

    std::string name() const
    { return getParam<std::string>("Problem.Name", "boussinesq_pressure"); }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        // pressure eq: Neumann everywhere (Robin penalty applied in neumann())
        values.setNeumann(conti0EqIdx);
        const Scalar zMax = this->gridGeometry().bBoxMax()[dimWorld-1];
        if (globalPos[dimWorld-1] > zMax - eps_)
            values.setDirichlet(conti0EqIdx + soluteIdx); // C = 1 at top
        else
            values.setNeumann(conti0EqIdx + soluteIdx);   // ∂C/∂n = 0 elsewhere
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition&) const
    {
        PrimaryVariables values(0.0);
        values[pressureIdx]           = pRef_;
        values[conti0EqIdx + soluteIdx] = 1.0;
        return values;
    }

    // Robin penalty BC for the pressure equation at the bottom face, used to weakly fix the
    // pressure level (the pure-Neumann elliptic pressure equation is otherwise singular up
    // to a constant). Pinned at a single vertex (x=0 corner) -- a real vertex sits exactly
    // there, so the `dofPosition()[0] < eps_` check selects exactly one DOF.
    NumEqVector neumann(const Element&,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache&,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);
        const Scalar zMin = this->gridGeometry().bBoxMin()[dimWorld-1];
        if (scvf.center()[dimWorld-1] < zMin + eps_)
        {
            const bool pinHere = fvGeometry.scv(scvf.insideScvIdx()).dofPosition()[0] < eps_;
            if (pinHere)
            {
                const Scalar p = elemVolVars[scvf.insideScvIdx()].pressure(0);
                values[conti0EqIdx] = (p - pRef_) * penalty_;
            }
        }
        return values;
    }

    // Called by BoussinesqCVFEDarcyLaw: returns C so that ρ_buoy = 1 + C
    Scalar boussinesq_term(const VolumeVariables& volVars) const
    {
        return volVars.massFraction(0, soluteIdx);
    }

    PrimaryVariables initialAtPos(const GlobalPosition&) const
    {
        PrimaryVariables values(0.0);
        return values; // p = 0, C = 0; erfc profile applied in main.cc
    }

private:
    static constexpr Scalar eps_ = 1e-7;
    Scalar pRef_;
    Scalar penalty_;
};

} // namespace Dumux

#endif
