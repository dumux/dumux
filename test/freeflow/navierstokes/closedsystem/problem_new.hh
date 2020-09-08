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
 * \brief A test problem for the staggered (Navier-) Stokes model.
 */

#ifndef DUMUX_CLOSEDSYSTEM_TEST_PROBLEM_HH
#define DUMUX_CLOSEDSYSTEM_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "../l2error.hh"

namespace Dumux {
template <class TypeTag>
class ClosedSystemTestProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct ClosedSystemTest {};
struct ClosedSystemTestMomentum { using InheritsFrom = std::tuple<ClosedSystemTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct ClosedSystemTestMass { using InheritsFrom = std::tuple<ClosedSystemTest, NavierStokesMassOneP, CCTpfaModel>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::ClosedSystemTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::ClosedSystemTest> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::ClosedSystemTest> { using type = Dumux::ClosedSystemTestProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::ClosedSystemTest> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::ClosedSystemTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::ClosedSystemTest> { static constexpr bool value = true; };
} // end namespace Properties

/*!
 * \ingroup NavierStokesTests
 * \brief Test problem for the one-phase (Navier-) Stokes model.
 *
 * Here, a quadratic two-dimensional domain with closed walls at all sides is considered.
 * If all walls are immobile and gravity is switched on, a hydrostatic pressure
 * gradient will develop. When assigning a fixed velocity to the top wall (without gravity),
 * this test corresponds to a lid-driven cavity problem.
 */
template <class TypeTag>
class ClosedSystemTestProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using NumEqVector = typename ParentType::NumEqVector;
    using PrimaryVariables = typename ParentType::PrimaryVariables;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    ClosedSystemTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        lidVelocity_ = getParam<Scalar>("Problem.LidVelocity");
    }

   /*!
     * \name Problem parameters
     */
    // \{

   /*!
     * \brief Returns the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10C

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

        // set Dirichlet values for the velocity everywhere
        values.setAllDirichlet();

        return values;
    }


   /*!
     * \brief Returns Dirichlet boundary values at a given position.
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);

        if constexpr (ParentType::isMomentumProblem())
            if (globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_)
                values[Indices::velocityXIdx] = lidVelocity_;

        return values;
    }

    /*!
     * \brief Returns a reference pressure at a given sub control volume face.
     *        This pressure is substracted from the actual pressure for the momentum balance
     *        which potentially helps to improve numerical accuracy by avoiding issues related do floating point arithmetic.
     * \note  Overload this for reference pressures other than zero.
     */
    Scalar referencePressure(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const SubControlVolumeFace& scvf) const
    { return 1e5; }

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

        auto isLowerLeftCell = [&](const SubControlVolume& scv)
        { return scv.dofIndex() == 0; };

        if (isLowerLeftCell(scv))
            values.set(0);

        return values;
    }

    /*!
     * \brief Define the values of internal Dirichlet constraints for a degree of freedom.
     * \param element The finite element
     * \param scv The sub-control volume
     */
    PrimaryVariables internalDirichlet(const Element& element, const SubControlVolume& scv) const
    { return PrimaryVariables(1.1e5); }

   /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;

        if constexpr (!ParentType::isMomentumProblem())
            values[Indices::pressureIdx] = 1.0e+5;

        return values;
    }

    // \}

private:

    static constexpr Scalar eps_=1e-6;
    Scalar lidVelocity_;
};
} // end namespace Dumux

#endif
