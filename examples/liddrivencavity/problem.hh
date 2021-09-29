// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/

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
#include <dumux/common/numeqvector.hh>

// Include the `NavierStokesProblem` class, the base
// class from which we will derive.
#include <dumux/freeflow/navierstokes/problem.hh>

// Include the `NavierStokesBoundaryTypes` class which specifies the boundary types set in this problem.
#include <dumux/freeflow/navierstokes/boundarytypes.hh>

// ### The problem class
// As we are solving a problem related to free flow, we create a new class called `LidDrivenCavityExampleProblem`
// and let it inherit from the class `NavierStokesProblem`.
// [[codeblock]]
namespace Dumux {
template <class TypeTag>
class LidDrivenCavityExampleProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PrimaryVariables = typename ParentType::PrimaryVariables;
    using NumEqVector = typename ParentType::NumEqVector;
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

    // #### Temperature distribution
    // We need to specify a constant temperature for our isothermal problem.
    // Fluid properties that depend on temperature will be calculated with this value.
    // This would be important if another fluidsystem was used.
    Scalar temperature() const
    { return 273.15 + 10; } // 10°C

    // #### Boundary conditions
    // With the following function we define the __type of boundary conditions__ depending on the location.
    // Three types of boundary conditions can be specified: Dirichlet, Neumann or outflow boundary conditions. On
    // Dirichlet boundaries, the values of the primary variables need to be fixed. On a Neumann boundaries,
    // values for derivatives need to be fixed. Outflow conditions set a gradient of zero in normal direction towards the boundary
    // for the respective primary variables (excluding pressure).
    // When Dirichlet conditions are set for the pressure, the velocity gradient
    // with respect to the direction normal to the boundary is automatically set to zero.
    // [[codeblock]]
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        // We set Dirichlet values for the velocity at each boundary
        if constexpr (ParentType::isMomentumProblem())
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        return values;
    }
    // [[/codeblock]]

    // The following function specifies the __values on Dirichlet boundaries__.
    // We need to define values for the primary variables (velocity and pressure).
    // [[codeblock]]
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);

        if constexpr (ParentType::isMomentumProblem())
        {
            if (globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_)
                values[Indices::velocityXIdx] = lidVelocity_;
        }
        else
            values[Indices::pressureIdx] = 1.1e+5;

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     * \param scvf The boundary sub control volume face
     */
    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        if constexpr (!ParentType::isMomentumProblem())
        {
            // Density is constant, so inside or outside does not matter.
            const auto insideDensity = elemVolVars[scvf.insideScvIdx()].density();

            // The resulting flux over the boundary is zero anyway, but this will add some non-zero derivatives to the
            // Jacobian which makes the BC more general.
            values[Indices::conti0EqIdx] = this->faceVelocity(element, fvGeometry, scvf) * insideDensity * scvf.unitOuterNormal();
        }

        return values;
    }
    // [[/codeblock]]

    /*!
     * \brief Returns a reference pressure at a given sub control volume face.
     *        This pressure is substracted from the actual pressure for the momentum balance
     *        which potentially helps to improve numerical accuracy by avoiding issues related do floating point arithmetic.
     */
    Scalar referencePressure(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const SubControlVolumeFace& scvf) const
    { return 1.0e5; }

    //! Enable internal Dirichlet constraints
    static constexpr bool enableInternalDirichletConstraints()
    { return !ParentType::isMomentumProblem(); }

    // We define a function for setting a fixed Dirichlet pressure value at a given internal cell.
    // This is required for having a defined pressure level in our closed system domain.
    /*!
     * \brief Tag a degree of freedom to carry internal Dirichlet constraints.
     *        If true is returned for a dof, the equation for this dof is replaced
     *        by the constraint that its primary variable values must match the
     *        user-defined values obtained from the function internalDirichlet(),
     *        which must be defined in the problem.
     *
     * \param element The finite element
     * \param scv The sub-control volume
     */
    std::bitset<PrimaryVariables::dimension> hasInternalDirichletConstraint(const Element& element, const SubControlVolume& scv) const
    {
        std::bitset<PrimaryVariables::dimension> values;

        const bool isLowerLeftCell = (scv.dofIndex() == 0);

        if constexpr (!ParentType::isMomentumProblem())
        {
            if (isLowerLeftCell)
                values.set(0);
        }

        return values;
    }

    /*!
     * \brief Define the values of internal Dirichlet constraints for a degree of freedom.
     * \param element The finite element
     * \param scv The sub-control volume
     */
    PrimaryVariables internalDirichlet(const Element& element, const SubControlVolume& scv) const
    { return PrimaryVariables(1.1e5); }

    // The following function defines the initial conditions.
    // [[codeblock]]
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);

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
