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
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Test for the staggered grid Navier-Stokes model with analytical solution (Kovasznay 1948, \cite Kovasznay1948)
 */

#ifndef DUMUX_KOVASZNAY_TEST_PROBLEM_HH
#define DUMUX_KOVASZNAY_TEST_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/problem.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the staggered grid (Kovasznay 1948, \cite Kovasznay1948)
 *
 * A two-dimensional Navier-Stokes flow with a periodicity in one direction
 * is considered. The set-up represents a wake behind a two-dimensional grid
 * and is chosen in a way such that an exact solution is available.
 */
template <class TypeTag>
class KovasznayTestProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using PrimaryVariables = typename ParentType::PrimaryVariables;

    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    static constexpr auto upwindSchemeOrder = getPropValue<TypeTag, Properties::UpwindSchemeOrder>();

public:
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    KovasznayTestProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                         std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        std::cout<< "upwindSchemeOrder is: " << upwindSchemeOrder << "\n";
        rho_ = getParam<Scalar>("Component.LiquidDensity", 1.0);
        kinematicViscosity_ = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0);
        Scalar reynoldsNumber = 1.0 / kinematicViscosity_;
        lambda_ = 0.5 * reynoldsNumber
                        - std::sqrt(reynoldsNumber * reynoldsNumber * 0.25 + 4.0 * M_PI * M_PI);
    }

   /*!
     * \brief Returns the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 298.0; }

    // \}
   /*!
     * \name Boundary conditions
     */
    // \{

   /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        if constexpr (ParentType::isMomentumProblem())
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
        }
        else
            values.setAllNeumann();

        return values;
    }

   /*!
     * \brief Returns Dirichlet boundary values at a given position.
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition & globalPos) const
    {
        // use the values of the analytical solution
        return analyticalSolution(globalPos);
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
            const auto insideDensity = elemVolVars[scvf.insideScvIdx()].density();
            values[Indices::conti0EqIdx] = insideDensity * (this->faceVelocity(element, fvGeometry, scvf) * scvf.unitOuterNormal());
        }
        return values;
    }
   /*!
     * \brief Returns the analytical solution of the problem at a given position.
     *
     * \param globalPos The global position
     * \param time A parameter for consistent signatures. It is ignored here as this is a stationary test
     */
    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos, Scalar time = 0.0) const
    {
        const Scalar x = globalPos[0];
        PrimaryVariables values;

        if constexpr (ParentType::isMomentumProblem())
        {
            const Scalar y = globalPos[1];
            values[Indices::velocityXIdx] = 1.0 - std::exp(lambda_ * x) * std::cos(2.0 * M_PI * y);
            values[Indices::velocityYIdx] = 0.5 * lambda_ / M_PI * std::exp(lambda_ * x) * std::sin(2.0 * M_PI * y);
        }
        else
        {   values[Indices::pressureIdx] = rho_ *0.5 * (1.0 - std::exp(2.0 * lambda_ * x));
        }


        return values;
    }

    // \}

   /*!
     * \name Volume terms
     */
    // \{

   /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        return PrimaryVariables(0.0);
    }

    //! Enable internal Dirichlet constraints
    static constexpr bool enableInternalDirichletConstraints()
    { return !ParentType::isMomentumProblem(); }

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

        const auto fvGeometry = localView(this->gridGeometry()).bindElement(element);
        const bool isAtLeftBoundary = [&]{
            if (fvGeometry.hasBoundaryScvf())
                 for (const auto& scvf : scvfs(fvGeometry))
                    if (scvf.boundary() && scvf.center()[0] < this->gridGeometry().bBoxMin()[0] + eps_)
                        return true;
            return false;
        }();

        if (isAtLeftBoundary)
            values.set(0);

        // TODO: only use one cell or pass fvGeometry to hasInternalDirichletConstraint
        // the pure Neumann problem is only defined up to a constant
        // we create a well-posed problem by fixing the pressure at one dof
        return values;
    }

    /*!
     * \brief Define the values of internal Dirichlet constraints for a degree of freedom.
     * \param element The finite element
     * \param scv The sub-control volume
     */
    PrimaryVariables internalDirichlet(const Element& element, const SubControlVolume& scv) const
    { return PrimaryVariables(analyticalSolution(scv.center())[Indices::pressureIdx]); }

private:
    static constexpr Scalar eps_=1e-6;

    Scalar rho_;
    Scalar kinematicViscosity_;
    Scalar lambda_;
};
} // end namespace Dumux

#endif
