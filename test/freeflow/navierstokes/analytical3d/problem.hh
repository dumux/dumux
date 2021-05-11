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
 * \brief 3D transient test for the Navier-Stokes model (Ethier and Steinmann, 1994).
 */

#ifndef DUMUX_NAVIERSTOKES_ANALYTICAL3D_PROBLEM_HH
#define DUMUX_NAVIERSTOKES_ANALYTICAL3D_PROBLEM_HH


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
class AnalyticalThreeDProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct Analytical3DTest {};
struct Analytical3DTestMomentum { using InheritsFrom = std::tuple<Analytical3DTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct Analytical3DTestMass { using InheritsFrom = std::tuple<Analytical3DTest, NavierStokesMassOneP, CCTpfaModel>; };
} // end namespace TTag


// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Analytical3DTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Analytical3DTest> { using type = Dune::SPGrid<double, 3>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Analytical3DTest> { using type = Dumux::AnalyticalThreeDProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::Analytical3DTest> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::Analytical3DTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::Analytical3DTest> { static constexpr bool value = true; };

} // end namespace Properties


/*!
 * \ingroup NavierStokesTests
 * \brief  Periodic test problem for the staggered grid
 *
 * A two-dimensional Navier-Stokes flow with a periodicity in one direction
 * is considered.
 */
template <class TypeTag>
class AnalyticalThreeDProblem :  public NavierStokesProblem<TypeTag>
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
    AnalyticalThreeDProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {}

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
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        return analyticalSolution(globalPos);
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
        const Scalar z = globalPos[2];

        using std::exp;
        const Scalar t = time;
        static const Scalar nu = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0);
        static constexpr Scalar a = M_PI/4.0;
        static constexpr Scalar d = M_PI/2.0;
        const Scalar f = exp(-(d*d)*nu*t);

        const Scalar ax = a*x;
        const Scalar ay = a*y;
        const Scalar az = a*z;
        const Scalar dx = d*x;
        const Scalar dy = d*y;
        const Scalar dz = d*z;

        PrimaryVariables values;

        if constexpr (ParentType::isMomentumProblem())
        {
            const auto eax = exp(ax);
            const auto eay = exp(ay);
            const auto eaz = exp(az);
            values[Indices::velocityXIdx] = -a * (eax*sin(ay + dz) + eaz*cos(ax + dy)) * f;
            values[Indices::velocityYIdx] = -a * (eay*sin(az + dx) + eax*cos(ay + dz)) * f;
            values[Indices::velocityZIdx] = -a * (eaz*sin(ax + dy) + eay*cos(az + dx)) * f;
        }
        else
        {
            values[Indices::pressureIdx] = -0.5*a*a * (exp(2*ax) + exp(2*ay) + exp(2*az)
                                                       + 2*sin(ax + dy)*cos(az + dx)*exp(a*(y+z))
                                                       + 2*sin(ay + dz)*cos(ax + dy)*exp(a*(z+x))
                                                       + 2*sin(az + dx)*cos(ay + dz)*exp(a*(x+y))) * exp(-2.0*(d*d)*nu*t);
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

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return analyticalSolution(globalPos, 0.0); }

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

        if (scv.dofIndex() == 0)
            values.set(0);

        return values;
    }

    /*!
     * \brief Define the values of internal Dirichlet constraints for a degree of freedom.
     * \param element The finite element
     * \param scv The sub-control volume
     */
    PrimaryVariables internalDirichlet(const Element& element, const SubControlVolume& scv) const
    {
        return analyticalSolution(scv.center());
    }


    /*!
     * \brief Updates the time
     */
    void updateTime(const Scalar time)
    { time_ = time; }

private:
    Scalar time_ = 0.0;
};
} // end namespace Dumux

#endif
