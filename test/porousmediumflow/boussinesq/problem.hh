// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief One-sided Rayleigh-Bénard dissolution problem (dimensionless Boussinesq).
 *
 * A 2D porous layer of dimensionless height H = 1 and aspect ratio Γ = L/H.
 * The top boundary is held at saturation concentration (c = 1) and a reference
 * pressure.  All other boundaries are no-flow / no-flux.
 *
 * The dissolved solute makes the fluid denser:
 *   ρ = 1 + Ra · c
 * with dimensionless gravity g = 1 and Ra = DimensionlessNumbers.Ra from the input file.
 * Convection sets in above the critical Rayleigh number Ra_c ≈ 4π².
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

    using BoundaryTypes  = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector    = Dumux::NumEqVector<PrimaryVariables>;

    static constexpr int dimWorld   = GridView::dimensionworld;
    static constexpr int pressureIdx = Indices::pressureIdx;
    static constexpr int soluteIdx   = FluidSystem::soluteIdx;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    BoussinesqOneSidedRBProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        FluidSystem::init();
        std::cout << "Boussinesq one-sided RB problem: Ra = "
                  << FluidSystem::rayleighNumber() << std::endl;
    }

    std::string name() const
    { return getParam<std::string>("Problem.Name", "boussinesq_onesided_rb"); }

    /*!
     * \brief Top boundary: Dirichlet (saturated solute, reference pressure).
     *        All other boundaries: Neumann (zero flux / no-flow).
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        const Scalar yMax = this->gridGeometry().bBoxMax()[dimWorld-1];
        if (globalPos[dimWorld-1] > yMax - eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();
        return values;
    }

    //! Saturated solute concentration and reference pressure at the top
    PrimaryVariables dirichletAtPos(const GlobalPosition&) const
    {
        PrimaryVariables values(0.0);
        values[pressureIdx] = pRef_;
        values[soluteIdx]   = 1.0; // dimensionless saturation concentration
        return values;
    }

    //! Zero flux everywhere except the top (which is Dirichlet)
    NumEqVector neumannAtPos(const GlobalPosition&) const
    { return NumEqVector(0.0); }

    /*!
     * \brief Initial condition: pure solvent with hydrostatic pressure.
     *
     * Hydrostatic pressure for ρ = 1, g = 1:
     *   p(y) = p_ref + (y_max - y)
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        const Scalar yMax = this->gridGeometry().bBoxMax()[dimWorld-1];
        PrimaryVariables values(0.0);
        values[pressureIdx] = pRef_ + (yMax - globalPos[dimWorld-1]);
        values[soluteIdx]   = 0.0; // pure solvent
        return values;
    }

private:
    static constexpr Scalar eps_  = 1e-7;
    static constexpr Scalar pRef_ = 1.0; // dimensionless reference pressure at top
};

} // namespace Dumux

#endif