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
 *   - ψ = 0 Dirichlet on ALL walls  (exact no-flow via streamfunction)
 *   - C = 1 Dirichlet at top        (saturated dissolving interface)
 *   - ∂C/∂n = 0 elsewhere           (no solute flux through solid walls)
 */
#ifndef DUMUX_BOUSSINESQ_ONESIDED_RB_PROBLEM_HH
#define DUMUX_BOUSSINESQ_ONESIDED_RB_PROBLEM_HH

#include <iostream>
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
    using ParentType   = PorousMediumFlowProblem<TypeTag>;
    using Scalar       = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView     = typename GridGeometry::GridView;
    using Indices      = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using BoundaryTypes    = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector      = Dumux::NumEqVector<PrimaryVariables>;

    static constexpr int dimWorld          = GridView::dimensionworld;
    static constexpr int nPot              = Indices::numPotentialEqs;
    static constexpr int concentrationIdx  = Indices::concentrationIdx;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    BoussinesqOneSidedRBProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        std::cout << "Boussinesq streamfunction problem: Ra = "
                  << getParam<Scalar>("DimensionlessNumbers.Ra") << "\n";
    }

    std::string name() const
    { return getParam<std::string>("Problem.Name", "boussinesq_onesided_rb"); }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        // A_k = 0 on every wall (exact no-flow for any gravity direction)
        for (int k = 0; k < nPot; ++k)
            values.setDirichlet(Indices::vectorPotentialIdx(k));
        const Scalar zMax = this->gridGeometry().bBoxMax()[dimWorld-1];
        if (globalPos[dimWorld-1] > zMax - eps_)
            values.setDirichlet(concentrationIdx); // C = 1 at top
        else
            values.setNeumann(concentrationIdx);   // ∂C/∂n = 0 elsewhere
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        for (int k = 0; k < nPot; ++k)
            values[Indices::vectorPotentialIdx(k)] = 0.0;
        const Scalar zMax = this->gridGeometry().bBoxMax()[dimWorld-1];
        values[concentrationIdx] = (globalPos[dimWorld-1] > zMax - eps_) ? 1.0 : 0.0;
        return values;
    }

    NumEqVector neumannAtPos(const GlobalPosition&) const
    { return NumEqVector(0.0); }

    PrimaryVariables initialAtPos(const GlobalPosition&) const
    {
        PrimaryVariables values(0.0);
        return values; // ψ=0, C=0; erfc profile applied in main.cc
    }

private:
    static constexpr Scalar eps_ = 1e-7;
};

} // namespace Dumux

#endif
