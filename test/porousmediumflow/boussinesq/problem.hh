// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief One-sided Rayleigh-Bénard dissolution problem (streamfunction-Boussinesq).
 *
 * Primary variables: ψ (streamfunction, slot 0) and C (solute concentration, slot 1).
 *
 * Boundary conditions:
 *   - ψ = 0 on ALL walls  (no-flow, exact — no penalty needed)
 *   - C = 1 at top        (saturated dissolving interface)
 *   - ∂C/∂n = 0 elsewhere (no solute flux through solid walls)
 *
 * The velocity field is u = (∂ψ/∂y, −∂ψ/∂x), so ψ = 0 at a wall means
 * zero normal AND tangential velocity there by construction.
 */
#ifndef DUMUX_BOUSSINESQ_ONESIDED_RB_PROBLEM_HH
#define DUMUX_BOUSSINESQ_ONESIDED_RB_PROBLEM_HH

#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

template<class TypeTag>
class BoussinesqOneSidedRBProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType  = PorousMediumFlowProblem<TypeTag>;
    using Scalar      = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView    = typename GridGeometry::GridView;
    using Indices     = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using BoundaryTypes   = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector     = Dumux::NumEqVector<PrimaryVariables>;

    static constexpr int dimWorld   = GridView::dimensionworld;
    static constexpr int psiIdx     = Indices::pressureIdx; // streamfunction slot
    static constexpr int soluteIdx  = FluidSystem::soluteIdx;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    BoussinesqOneSidedRBProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        FluidSystem::init();
        std::cout << "Boussinesq streamfunction problem: Ra = "
                  << FluidSystem::rayleighNumber() << std::endl;
    }

    std::string name() const
    { return getParam<std::string>("Problem.Name", "boussinesq_onesided_rb"); }

    /*!
     * \brief All walls: ψ = 0 Dirichlet (no-flow).
     *        Top:       also C = 1 Dirichlet (dissolving interface).
     *        Other walls: ∂C/∂n = 0 (natural Neumann for transport eq).
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        // ψ = 0 on every wall
        values.setDirichlet(psiIdx);
        // C = 1 only at the top; natural (zero-flux) elsewhere
        const Scalar yMax = this->gridGeometry().bBoxMax()[dimWorld-1];
        if (globalPos[dimWorld-1] > yMax - eps_)
            values.setDirichlet(soluteIdx);
        else
            values.setNeumann(soluteIdx);
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        values[psiIdx] = 0.0; // ψ = 0
        const Scalar yMax = this->gridGeometry().bBoxMax()[dimWorld-1];
        values[soluteIdx] = (globalPos[dimWorld-1] > yMax - eps_) ? 1.0 : 0.0;
        return values;
    }

    NumEqVector neumannAtPos(const GlobalPosition&) const
    { return NumEqVector(0.0); }

    /*!
     * \brief Initial condition: no flow (ψ = 0), pure solvent (C = 0).
     *
     * The diffusive erfc-profile starting condition is applied in main.cc
     * by overwriting the C DOFs after applyInitialSolution().
     */
    PrimaryVariables initialAtPos(const GlobalPosition&) const
    {
        PrimaryVariables values(0.0);
        values[psiIdx]   = 0.0;
        values[soluteIdx] = 0.0;
        return values;
    }

private:
    static constexpr Scalar eps_ = 1e-7;
};

} // namespace Dumux

#endif
