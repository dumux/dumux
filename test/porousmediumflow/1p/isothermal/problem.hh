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
 * \ingroup OnePTests
 * \brief A test problem for the one-phase model:
 * water is flowing from bottom to top through and around a low permeable lens.
 */

#ifndef DUMUX_1PTEST_PROBLEM_HH
#define DUMUX_1PTEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/boundarytypes.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "spatialparams.hh"

#if FORCHHEIMER
#include <dumux/flux/forchheimerslaw.hh>
#endif

namespace Dumux {

template <class TypeTag>
class OnePTestProblem;

namespace Properties {

// Create new type tags
namespace TTag {
struct OnePTest { using InheritsFrom = std::tuple<OneP>; };
struct OnePTestBox { using InheritsFrom = std::tuple<BoxModel, OnePTest>; };
struct OnePTestCCTpfa { using InheritsFrom = std::tuple<CCTpfaModel, OnePTest>; };
struct OnePTestCCMpfa { using InheritsFrom = std::tuple<CCMpfaModel, OnePTest>; };
} // end namespace TTag

// Specialize the fluid system type for this type tag
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePTest>
{
    using Scalar = GetPropType<TypeTag, Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
};

// Specialize the grid type for this type tag
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePTest>
{ using type = Dune::YaspGrid<2>; };

// Specialize the problem type for this type tag
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePTest>
{ using type = OnePTestProblem<TypeTag>; };

// Specialize the spatial params type for this type tag
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePTest>
{
    using GridGeometry = GetPropType<TypeTag, GridGeometry>;
    using Scalar = GetPropType<TypeTag, Scalar>;
    using type = OnePTestSpatialParams<GridGeometry, Scalar>;
};

#ifdef FORCHHEIMER
// Specialize the advection type for this type tag
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::OnePTest>
{ using type = ForchheimersLaw<TypeTag>; };
#endif

} // end namespace Properties

/*!
 * \ingroup OnePTests
 * \brief  Test problem for the one-phase model:
 * water is flowing from bottom to top through and around a low permeable lens.
 *
 * The domain is box shaped. All sides are closed (Neumann 0 boundary)
 * except the top and bottom boundaries (Dirichlet), where water is
 * flowing from bottom to top.
 *
 * In the middle of the domain, a lens with low permeability (\f$K=10e-12\f$)
 * compared to the surrounding material (\f$ K=10e-10\f$) is defined.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box1p -parameterFile test_box1p.input</tt> or
 * <tt>./test_cc1p -parameterFile test_cc1p.input</tt>
 *
 * The same parameter file can be also used for 3d simulation but you need to change line
 * <tt>using type = Dune::YaspGrid<2>;</tt> to
 * <tt>using type = Dune::YaspGrid<3>;</tt> in the problem file
 * and use <tt>test_1p_3d.dgf</tt> in the parameter file.
 */
template <class TypeTag>
class OnePTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    // copy some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    enum {
        // index of the primary variable
        pressureIdx = Indices::pressureIdx
    };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;

    static constexpr int dimWorld = GridView::dimensionworld;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    OnePTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        name_ = getParam<std::string>("Problem.Name");
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    std::string name() const
    {
        return name_;
    }

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

        if (globalPos[dimWorld-1] < eps_ || globalPos[dimWorld-1] > this->gridGeometry().bBoxMax()[dimWorld-1] - eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0);
        values[pressureIdx] = 1.0e+5*(2.0 - globalPos[dimWorld-1]);
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
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars[pressureIdx] = 1.0e+5;
        return priVars;
    }

    // \}

private:
    std::string name_;
    static constexpr Scalar eps_ = 1.0e-6;
};

} // end namespace Dumux

#endif
