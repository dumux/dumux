// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_EXAMPLES_FREEFLOW_CHANNEL_PROBLEM_HH
#define DUMUX_EXAMPLES_FREEFLOW_CHANNEL_PROBLEM_HH

// ## Initial and boundary conditions (`problem.hh`)
//
// This file contains the __problem class__ which defines the initial and boundary
// conditions for the Navier-Stokes single-phase flow simulation.
//
// [[content]]
//
// ### Include files
//
// Include the `NavierStokesBoundaryTypes` class which specifies the boundary types set in this problem.
#include <dumux/freeflow/navierstokes/boundarytypes.hh>
//
// Include helper functions to compute values for boundary conditions
#include <dumux/freeflow/navierstokes/momentum/fluxhelper.hh>
#include <dumux/freeflow/navierstokes/scalarfluxhelper.hh>
#include <dumux/freeflow/navierstokes/mass/1p/advectiveflux.hh>

// ### The problem class
// We enter the problem class `ChannelExampleProblem` where all necessary boundary conditions and initial conditions are set for our simulation.
// As we are solving a problem related to free flow using a coupled model, we inherit from the base
// class of the respective subproblem for momentum or mass balance (see `properties.hh`).
// [[codeblock]]
namespace Dumux {

template <class TypeTag, class BaseProblem>
class ChannelExampleProblem : public BaseProblem
{
    // A few convenience aliases used throughout this class.
    using ParentType = BaseProblem;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using InitialValues = typename ParentType::InitialValues;
    using Sources = typename ParentType::Sources;
    using DirichletValues = typename ParentType::DirichletValues;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    // This is the constructor of our problem class:
    // Within the constructor, we set the inlet velocity to a run-time specified value.
    // If no run-time value is specified, we set the outlet pressure to 1.1e5 Pa.
    ChannelExampleProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                          std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        inletVelocity_ = getParam<Scalar>("Problem.InletVelocity");
        outletPressure_ = getParam<Scalar>("Problem.OutletPressure", 1.1e5);
    }
    // [[/codeblock]]

    // #### Boundary conditions
    // With the following function we define the __type of boundary conditions__ depending on the location.
    // Two types of boundary conditions can be specified: Dirichlet or Neumann. On
    // Dirichlet boundaries, the values of the primary variables need to be fixed. On Neumann boundaries,
    // the flux needs to be fixed.
    // To set different conditions for the two subproblems, use `constexpr ParentType::isMomentumProblem()`
    // to distinguish between momentum and mass problem.
    // To set Dirichlet conditions for the pressure, instead specify a solution-dependent Neumann
    // condition for the momentum balance, which depends on the pressure, using the helper function.
    // [[codeblock]]
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;

        if constexpr(ParentType::isMomentumProblem())
        {
            // We specify Dirichlet boundary conditions for the velocity on most boundaries of our domain
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);

            if (isOutlet_(globalPos))
            {
                // We fix the pressure on the right side of the domain, for the momentum balance we compute the resulting flux
                values.setAllNeumann();
            }
        }
        else
        {
            if (isInlet_(globalPos))
            {
                // We specify Dirichlet boundary conditions for the velocity on the left of our
                // domain, the corresponding pressure can be obtained from the coupling manager
                values.setDirichlet(Indices::pressureIdx);
            }
            else if (isOutlet_(globalPos))
            {
                // We fix the pressure on the right side of the domain through the momentum outflow,
                // for the mass balance we may prescribe a pressure or a mass outflow computed from velocity fields
                values.setNeumann(Indices::conti0EqIdx);
            }
            else
            {
                // We specify no-flow Neumann boundary conditions for the mass balance on the remaining boundaries (lower and upper wall)
                // in addition to the Dirichlet boundary conditions for the velocity (momentum balance)
                values.setAllNeumann();
            }
        }

        return values;
    }
    // [[/codeblock]]

    // The following function specifies the __fluxes on Neumann boundaries__.
    // We need to define fluxes for the balance equations (momentum or mass).
    // [[codeblock]]
    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    BoundaryFluxes neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVariablesCache& elemFluxVarsCache,
                           const SubControlVolumeFace& scvf) const
    {
        // No flow as default
        BoundaryFluxes values(0.0);

        if constexpr (ParentType::isMomentumProblem())
        {
            // Compute the solution-dependent momentum flux for the specified pressure and zero normal velocity gradient
            using FluxHelper = NavierStokesMomentumBoundaryFlux<typename GridGeometry::DiscretizationMethod>;
            values = FluxHelper::fixedPressureMomentumFlux(*this, fvGeometry, scvf, elemVolVars, elemFluxVarsCache, outletPressure_, true /*zeroNormalVelocityGradient*/);
        }
        else
        {
            // Compute the solution-dependent mass flux based on velocity fields
            using FluxHelper = NavierStokesScalarBoundaryFluxHelper<AdvectiveFlux<ModelTraits>>;
            if (isOutlet_(scvf.ipGlobal()))
                values = FluxHelper::scalarOutflowFlux(*this, element, fvGeometry, scvf, elemVolVars);
        }

        return values;
    }
    // [[/codeblock]]

    // The following function specifies the __values on Dirichlet boundaries__.
    // We need to define values for the primary variables (velocity or pressure).
    // [[codeblock]]
    DirichletValues dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    {
        const auto& globalPos = scvf.ipGlobal();
        // Use the initial values as default Dirichlet values
        DirichletValues values = initialAtPos(globalPos);

        if constexpr (ParentType::isMomentumProblem())
        {
            // Set a no-slip condition at the top and bottom wall of the channel
            if (!isInlet_(globalPos))
                values[Indices::velocityXIdx] = 0.0;
        }
        else
        {
            if (isInlet_(globalPos))
                values = this->couplingManager().cellPressure(element, scvf);
        }

        return values;
    }
    // [[/codeblock]]

    // The following function defines the initial conditions.
    // [[codeblock]]
    InitialValues initialAtPos(const GlobalPosition& globalPos) const
    {
        InitialValues values;

        // Set the pressure and velocity values
        if constexpr (ParentType::isMomentumProblem())
        {
            values[Indices::velocityXIdx] = inletVelocity_;
            values[Indices::velocityYIdx] = 0.0;
        }
        else
        {
            //std::cout << "setting outlet pressure at " << globalPos << std::endl;
            values[Indices::pressureIdx] = outletPressure_;
        }

        return values;
    }
    // [[/codeblock]]

// The inlet is on the left side of the physical domain.
// [[codeblock]]
private:
    bool isInlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] < eps_; }
    // [[/codeblock]]

    // The outlet is on the right side of the physical domain.
    bool isOutlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    // Finally, private variables are declared:
    // [[codeblock]]
    static constexpr Scalar eps_ = 1e-6;
    Scalar inletVelocity_;
    Scalar outletPressure_;

}; // end class definition of ChannelExampleProblem
} // end namespace Dumux
// [[/codeblock]]
// [[/content]]
#endif
