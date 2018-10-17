// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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

#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/cellcentered/mpfa/properties.hh>
#include <dumux/discretization/box/properties.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "1ptestspatialparams.hh"

#if FORCHHEIMER
#include <dumux/discretization/forchheimerslaw.hh>
#endif

namespace Dumux {

template <class TypeTag>
class OnePTestProblem;

namespace Properties {

// Create new type tags
namespace TTag {
struct OnePTestTypeTag { using InheritsFrom = std::tuple<OneP>; };
struct OnePTestBoxTypeTag { using InheritsFrom = std::tuple<OnePTestTypeTag, BoxModel>; };
struct OnePTestCCTpfaTypeTag { using InheritsFrom = std::tuple<OnePTestTypeTag, CCTpfaModel>; };
struct OnePTestCCMpfaTypeTag { using InheritsFrom = std::tuple<OnePTestTypeTag, CCMpfaModel>; };
} // end namespace TTag

// Specialize the fluid system type for this type tag
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePTestTypeTag>
{
    using Scalar = GetPropType<TypeTag, Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
};

// Specialize the grid type for this type tag
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePTestTypeTag>
{ using type = Dune::YaspGrid<2>; };

// Specialize the problem type for this type tag
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePTestTypeTag>
{ using type = OnePTestProblem<TypeTag>; };

// Specialize the spatial params type for this type tag
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePTestTypeTag>
{
    using FVGridGeometry = GetPropType<TypeTag, FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Scalar>;
    using type = OnePTestSpatialParams<FVGridGeometry, Scalar>;
};

#ifdef FORCHHEIMER
// Specialize the advection type for this type tag
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::OnePTestTypeTag>
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
    using GridView = Properties::GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = Properties::GetPropType<TypeTag, Properties::Scalar>;

    // copy some indices for convenience
    using Indices = typename Properties::GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    enum {
        // index of the primary variable
        pressureIdx = Indices::pressureIdx
    };

    using PrimaryVariables = Properties::GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Properties::GetPropType<TypeTag, Properties::BoundaryTypes>;
    using FVGridGeometry = Properties::GetPropType<TypeTag, Properties::FVGridGeometry>;

    static constexpr int dimWorld = GridView::dimensionworld;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    OnePTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
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
     * \brief Return the temperature within the domain in [K].
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
     * \param values The boundary types for the conservation equations
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        if (globalPos[dimWorld-1] < eps_ || globalPos[dimWorld-1] > this->fvGridGeometry().bBoxMax()[dimWorld-1] - eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
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
     * \brief Evaluate the initial value for a control volume.
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

} //end namespace Dumux

#endif
