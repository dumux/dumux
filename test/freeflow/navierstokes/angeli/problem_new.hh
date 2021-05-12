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
 * \brief Test for the instationary staggered grid Navier-Stokes model
 *        with analytical solution (Angeli et al. 2017, \cite Angeli2017).
 */

#ifndef DUMUX_ANGELI_TEST_PROBLEM_HH
#define DUMUX_ANGELI_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "../l2error.hh"

namespace Dumux {
template <class TypeTag>
class AngeliTestProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct AngeliTest {};
struct AngeliTestMomentum { using InheritsFrom = std::tuple<AngeliTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct AngeliTestMass { using InheritsFrom = std::tuple<AngeliTest, NavierStokesMassOneP, CCTpfaModel>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::AngeliTest>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::AngeliTest> { using type = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::AngeliTest> { using type = Dumux::AngeliTestProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::AngeliTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::AngeliTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::AngeliTest> { static constexpr bool value = true; };
} // end namespace Properties

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the staggered grid (Angeli et al. 2017, \cite Angeli2017).
 *
 * The unsteady, 2D, incompressible Navier-Stokes equations for a zero source and a Newtonian
 * flow is solved and compared to an analytical solution (sums/products of trigonometric functions).
 * The velocities and pressures decay exponentially. The Dirichlet boundary conditions are
 * time-dependent and consistent with the analytical solution.
 */
template <class TypeTag>
class AngeliTestProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using NumEqVector = typename ParentType::NumEqVector;
    using PrimaryVariables = typename ParentType::PrimaryVariables;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using TimeLoopPtr = std::shared_ptr<TimeLoop<Scalar>>;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using FVElementGeometry = typename GridGeometry::LocalView;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    AngeliTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        kinematicViscosity_ = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0);
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
            // set Dirichlet values for the velocity everywhere
            if (globalPos[0] < this->gridGeometry().bBoxMax()[0] - 1e-8)
            {
                values.setDirichlet(Indices::velocityXIdx);
                values.setDirichlet(Indices::velocityYIdx);
            }
            else
                values.setAllNeumann();
        }
        else
            values.setAllNeumann();

        return values;
    }

    //! Enable internal Dirichlet constraints
    static constexpr bool enableInternalDirichletConstraints()
    { return false; }
    // { return !ParentType::isMomentumProblem(); }

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

        // auto fvGeometry = localView(this->gridGeometry());
        // fvGeometry.bindElement(element);

        // auto isAtBoundary = [&](const FVElementGeometry& fvGeometry)
        // {
        //     if (fvGeometry.hasBoundaryScvf())
        //         return true;
        //     return false;
        // };

        // if (isAtBoundary(fvGeometry))
        //     values.set(Indices::pressureIdx);
        return values;
    }

    /*!
     * \brief Define the values of internal Dirichlet constraints for a degree of freedom.
     * \param element The finite element
     * \param scv The sub-control volume
     */
    PrimaryVariables internalDirichlet(const Element& element, const SubControlVolume& scv) const
    { return dirichletAtPos(scv.center()); }


   /*!
     * \brief Returns Dirichlet boundary values at a given position.
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        // use the values of the analytical solution
        return analyticalSolution(globalPos, time_);
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
            values[Indices::conti0EqIdx] = this->faceVelocity(element, fvGeometry, scvf) * insideDensity * scvf.unitOuterNormal();
        }
        else
        {
            using std::sin;
            using std::cos;
            Dune::FieldMatrix<Scalar, 2, 2> momentumFlux(0.0);
            const auto x = scvf.ipGlobal()[0];
            const auto y = scvf.ipGlobal()[1];
            const Scalar mu = kinematicViscosity_;
            const Scalar t = time_;
            momentumFlux[0][0] = M_PI*M_PI *(-4.0*mu*exp(10.0*M_PI*M_PI*mu*t)*sin(M_PI*x)*sin(2*M_PI*y) + (4.0*sin(2*M_PI*y)*sin(2*M_PI*y)*cos(M_PI*x)*cos(M_PI*x) - 1.0*cos(2*M_PI*x) - 0.25*cos(4*M_PI*y))*exp(5.0*M_PI*M_PI*mu*t))*exp(-15.0*M_PI*M_PI*mu*t);
            momentumFlux[0][1] = M_PI*M_PI *( 3.0*mu*exp(10.0*M_PI*M_PI*mu*t) - 2.0*exp(5.0*M_PI*M_PI*mu*t)*sin(M_PI*x)*sin(2*M_PI*y))*exp(-15.0*M_PI*M_PI*mu*t)*cos(M_PI*x)*cos(2*M_PI*y);
            momentumFlux[1][0] = M_PI*M_PI *( 3.0*mu*exp(10.0*M_PI*M_PI*mu*t) - 2.0*exp(5.0*M_PI*M_PI*mu*t)*sin(M_PI*x)*sin(2*M_PI*y))*exp(-15.0*M_PI*M_PI*mu*t)*cos(M_PI*x)*cos(2*M_PI*y);
            momentumFlux[1][1] = M_PI*M_PI *( 4.0*mu*exp(10.0*M_PI*M_PI*mu*t)*sin(M_PI*x)*sin(2*M_PI*y) + (sin(M_PI*x)*sin(M_PI*x)*cos(2*M_PI*y)*cos(2*M_PI*y) - 1.0*cos(2*M_PI*x) - 0.25*cos(4*M_PI*y))*exp(5.0*M_PI*M_PI*mu*t))*exp(-15.0*M_PI*M_PI*mu*t);

            const auto normal = scvf.unitOuterNormal();
            momentumFlux.mv(normal, values);
        }

        return values;
    }

    /*!
     * \brief Returns the analytical solution of the problem at a given time and position.
     *
     * \param globalPos The global position
     * \param time The current simulation time
     */
    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos, const Scalar time) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        const Scalar t = time;

        PrimaryVariables values;

        if constexpr (ParentType::isMomentumProblem())
        {
            values[Indices::velocityXIdx] = - 2.0 * M_PI * std::exp(- 5.0 * kinematicViscosity_ * M_PI * M_PI * t) * std::cos(M_PI * x) * std::sin(2.0 * M_PI * y);
            values[Indices::velocityYIdx] = M_PI * std::exp(- 5.0 * kinematicViscosity_ * M_PI * M_PI * t) * std::sin(M_PI * x) * std::cos(2.0 * M_PI * y);
        }
        else
        {
            values[Indices::pressureIdx] = - 0.25 * std::exp(-10.0 * kinematicViscosity_ * M_PI * M_PI * t) * M_PI * M_PI * (4.0 * std::cos(2.0 * M_PI * x) + std::cos(4.0 * M_PI * y));
        }

        return values;
    }

    /*!
     * \brief Returns the analytical solution of the problem at a given position.
     *
     * \param globalPos The global position
     */
    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos) const
    { return analyticalSolution(globalPos, time_); }

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
    { return analyticalSolution(globalPos, 0.0); }

    /*!
     * \brief Updates the time
     */
    void updateTime(const Scalar time)
    { time_ = time; }

private:
    Scalar kinematicViscosity_;
    Scalar time_ = 0.0;
};
} // end namespace Dumux

#endif
