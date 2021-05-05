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

#ifndef DUMUX_KOVASZNAY_TEST_PROBLEM_NEW_HH
#define DUMUX_KOVASZNAY_TEST_PROBLEM_NEW_HH


#include <dune/grid/spgrid.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>

namespace Dumux {

template <class TypeTag>
class PeriodicTestProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct PeriodicTest {};
struct PeriodicTestMomentum { using InheritsFrom = std::tuple<PeriodicTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct PeriodicTestMass { using InheritsFrom = std::tuple<PeriodicTest, NavierStokesMassOneP, CCTpfaModel>; };
} // end namespace TTag


// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PeriodicTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::PeriodicTest> { using type = Dune::SPGrid<double, 2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PeriodicTest> { using type = Dumux::PeriodicTestProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::PeriodicTest> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::PeriodicTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::PeriodicTest> { static constexpr bool value = true; };

} // end namespace Properties


/*!
 * \ingroup NavierStokesTests
 * \brief  Periodic test problem for the staggered grid
 *
 * A two-dimensional Navier-Stokes flow with a periodicity in one direction
 * is considered.
 */
template <class TypeTag>
class PeriodicTestProblem :  public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;
    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using NumEqVector = typename ParentType::NumEqVector;
    using PrimaryVariables = typename ParentType::PrimaryVariables;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    PeriodicTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        usePressureDifference_ = getParam<bool>("Problem.UsePressureDifference", false);
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

        values.setAllDirichlet(); // does not really make sense for the mass balance
        return values;
    }

   /*!
     * \brief Returns Dirichlet boundary values at a given position.
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition & globalPos) const
    {
        return PrimaryVariables(0.0);
    }

    template<class ElementVolumeVariables>
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        NumEqVector source;

        if constexpr (ParentType::isMomentumProblem())
        {

            if (usePressureDifference_ && scv.dofPosition()[1] < this->gridGeometry().bBoxMin()[1] + eps_)
            {
                const auto& frontalScvf = (*scvfs(fvGeometry, scv).begin());
                source[Indices::momentumYBalanceIdx] = 100 * frontalScvf.area() / scv.volume();
            }
        }

        return source;
    }

    // \}

   /*!
     * \name Volume terms
     */
    // \{

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

        for (const auto& intersection : intersections(this->gridGeometry().gridView(), element))
        {
            if (intersection.boundary() && intersection.neighbor() && intersection.geometry().center()[1] > this->gridGeometry().bBoxMax()[1] - eps_)
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
    {
        return PrimaryVariables(1.0);
    }

private:
    static constexpr Scalar eps_ = 1e-6;
    bool usePressureDifference_;
};
} // end namespace Dumux

#endif
