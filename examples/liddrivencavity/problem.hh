// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_LIDDRIVENCAVITY_EXAMPLE_PROBLEM_HH
#define DUMUX_LIDDRIVENCAVITY_EXAMPLE_PROBLEM_HH

// ## Initial and boundary conditions (`problem.hh`)
//
// This file contains the __problem class__ which defines the initial and boundary
// conditions for the Navier-Stokes single-phase flow simulation.
//
// [[content]]
//
// ### Include files
//
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

// Include the `NavierStokesBoundaryTypes` class which specifies the boundary types set in this problem.
#include <dumux/freeflow/navierstokes/boundarytypes.hh>

// ### The problem class
// As we are solving a problem related to free flow, we create a new class called `LidDrivenCavityExampleProblem`
// and let it inherit from a base class for the momentum and mass subproblems (selected in properties.hh).
// [[codeblock]]
namespace Dumux {
template <class TypeTag, class BaseProblem>
class LidDrivenCavityExampleProblem : public BaseProblem
{
    using ParentType = BaseProblem;

    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using InitialValues = typename ParentType::InitialValues;
    using Sources = typename ParentType::Sources;
    using DirichletValues = typename ParentType::DirichletValues;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename FVElementGeometry::Element;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    // Within the constructor, we set the lid velocity to a run-time specified value.
    LidDrivenCavityExampleProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        lidVelocity_ = getParam<Scalar>("Problem.LidVelocity");
    }
    // [[/codeblock]]

    // #### Boundary conditions
    // With the following function we define the __type of boundary conditions__ depending on the location.
    // Three types of boundary conditions can be specified: Dirichlet or Neumann boundary conditions. On
    // Dirichlet boundaries, the values of the primary variables need to be fixed. On a Neumann boundaries,
    // values for derivatives need to be fixed.
    // [[codeblock]]
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        // We set Dirichlet values for the velocity at each boundary. At the same time,
        // Neumann (no-flow) conditions hold at the boundaries for the mass model.
        if constexpr (ParentType::isMomentumProblem())
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        return values;
    }
    // [[/codeblock]]

    // The following function specifies the __values on Dirichlet boundaries__.
    // We need to define values for the primary variables (velocity).
    // [[codeblock]]
    DirichletValues dirichletAtPos(const GlobalPosition &globalPos) const
    {
        DirichletValues values(0.0);

        if constexpr (ParentType::isMomentumProblem())
        {
            if (globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_)
                values[Indices::velocityXIdx] = lidVelocity_;
        }

        return values;
    }
    // [[/codeblock]]

    // The following function specifies the __values on Neumann boundaries__.
    // We define a (zero) mass flux here.
    // [[codeblock]]
    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    BoundaryFluxes neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVariablesCache& elemFluxVarsCache,
                           const SubControlVolumeFace& scvf) const
    {
        BoundaryFluxes values(0.0);

        if constexpr (!ParentType::isMomentumProblem())
        {
            // Density is constant, so inside or outside does not matter.
            const auto insideDensity = elemVolVars[scvf.insideScvIdx()].density();

            // The resulting flux over the boundary is zero anyway (velocity is zero), but this will add some non-zero derivatives to the
            // Jacobian and makes the BC more general.
            values[Indices::conti0EqIdx] = this->faceVelocity(element, fvGeometry, scvf) * insideDensity * scvf.unitOuterNormal();
        }

        return values;
    }
    // [[/codeblock]]

    // The problem setup considers closed boundaries everywhere. In order to have a defined pressure level, we impose an __internal Dirichlet
    // constraint for pressure__ in a single cell.
    // [[codeblock]]

    // Use internal Dirichlet constraints for the mass problem.
    static constexpr bool enableInternalDirichletConstraints()
    { return !ParentType::isMomentumProblem(); }

    // Set a fixed pressure a the lower-left cell.
    std::bitset<DirichletValues::dimension> hasInternalDirichletConstraint(const Element& element, const SubControlVolume& scv) const
    {
        std::bitset<DirichletValues::dimension> values;

        if constexpr (!ParentType::isMomentumProblem())
        {
            const bool isLowerLeftCell = (scv.dofIndex() == 0);
            if (isLowerLeftCell)
                values.set(0);
        }

        return values;
    }

    // Specify the pressure value in the internal Dirichlet cell.
    DirichletValues internalDirichlet(const Element& element, const SubControlVolume& scv) const
    { return DirichletValues(1.1e5); }
    // [[/codeblock]]

    // Setting a __reference pressure__ can help to improve the Newton convergence rate by making the numerical derivatives more exact.
    // This is related to floating point arithmetic as pressure values are usually much higher than velocities.
    // [[codeblock]]
    Scalar referencePressure(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const SubControlVolumeFace& scvf) const
    { return 1.0e5; }
    // [[/codeblock]]

    // The following function defines the initial conditions.
    // [[codeblock]]
    InitialValues initialAtPos(const GlobalPosition &globalPos) const
    {
        InitialValues values(0.0);

        if constexpr (!ParentType::isMomentumProblem())
            values[Indices::pressureIdx] = 1.0e+5;

        return values;
    }
    // [[/codeblock]]
    // the data members of the problem class
    // [[codeblock]]
private:
    static constexpr Scalar eps_ = 1e-6;
    Scalar lidVelocity_;
};

} // end namespace Dumux
// [[/codeblock]]
// [[/content]]
#endif
