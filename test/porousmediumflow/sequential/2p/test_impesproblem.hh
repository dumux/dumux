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
 * \ingroup SequentialTwoPTests
 * \brief test problem for the sequential 2p model
 */
#ifndef DUMUX_TEST_IMPES_PROBLEM_HH
#define DUMUX_TEST_IMPES_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/properties.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/simpleh2o.hh>

#include <dumux/porousmediumflow/2p/sequential/diffusion/cellcentered/pressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/properties.hh>
#include <dumux/porousmediumflow/2p/sequential/impes/problem.hh>

//following includes are only needed if a global pressure formulation is chosen!
//Then only a total velocity can be reconstructed for the transport step
#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/capillarydiffusion.hh>
#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/gravitypart.hh>

#include "test_impesspatialparams.hh"

#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/evalcflfluxcoats.hh>

#include <dumux/linear/amgbackend.hh>
#include <dumux/linear/linearsolvertraits.hh>

namespace Dumux
{

template<class TypeTag>
class IMPESTestProblem;

//////////
// Specify the properties
//////////
namespace Properties
{
// Create new type tags
namespace TTag {
struct IMPESTest { using InheritsFrom = std::tuple<TestIMPESSpatialParams, FVTransportTwoP, IMPESTwoP, FVPressureTwoP>; };

// set up an additional problem where the AMG backend is used
struct IMPESTestWithAMG { using InheritsFrom = std::tuple<IMPESTest>; };

} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::IMPESTest> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::IMPESTest> { using type = IMPESTestProblem<TypeTag>; };

////////////////////////////////////////////////////////////////////////
//Switch to a p_n-S_w formulation
//
// template<class TypeTag>
// struct Formulation<TypeTag, TTag::IMPESTest> { static constexpr int value = SequentialTwoPCommonIndices::pnsn; };
//
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//Switch to a p_global-S_w formulation
//
// template<class TypeTag>
// struct Formulation<TypeTag, TTag::IMPESTest> { static constexpr int value = SequentialTwoPCommonIndices::pGlobalSw; };
//
//Define the capillary pressure term in the transport equation -> only needed in case of a p_global-S_w formulation!
template<class TypeTag>
struct CapillaryFlux<TypeTag, TTag::IMPESTest> { using type = CapillaryDiffusion<TypeTag>; };
//
//Define the gravity term in the transport equation -> only needed in case of a p_global-S_w formulation!
template<class TypeTag>
struct GravityFlux<TypeTag, TTag::IMPESTest> { using type = GravityPart<TypeTag>; };
//
////////////////////////////////////////////////////////////////////////

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::IMPESTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

template<class TypeTag>
struct EvalCflFluxFunction<TypeTag, TTag::IMPESTest> { using type = EvalCflFluxCoats<TypeTag>; };

// use the AMG backend for the corresponding test
template<class TypeTag>
struct LinearSolver<TypeTag, TTag::IMPESTestWithAMG>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using type = AMGBiCGSTABBackend<LinearSolverTraits<GridGeometry>>;

};
// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::IMPESTestWithAMG> { using type = Dune::YaspGrid<2>; };

} // end namespace Properties

/*!
 * \ingroup SequentialTwoPTests
 * \brief test problem for the sequential 2p model
 *
 * Water is injected from the left side into a rectangular 2D domain also
 * filled with water. Upper and lower boundary is closed (Neumann = 0),
 * and there is free outflow on the right side.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_impes -parameterFile ./test_impes.input</tt>,
 * where the arguments define the parameter file..
 */
template<class TypeTag>
class IMPESTestProblem: public IMPESProblem2P<TypeTag>
{
using ParentType = IMPESProblem2P<TypeTag>;
using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
using Grid = typename GridView::Grid;

using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

using WettingPhase = typename GetProp<TypeTag, Properties::FluidSystem>::WettingPhase;

using TimeManager = GetPropType<TypeTag, Properties::TimeManager>;

enum
{
    dim = GridView::dimension, dimWorld = GridView::dimensionworld
};

enum
{
    nPhaseIdx = Indices::nPhaseIdx,
    pwIdx = Indices::pwIdx,
    swIdx = Indices::swIdx,
    eqIdxPress = Indices::pressureEqIdx,
    eqIdxSat = Indices::satEqIdx
};

using Scalar = GetPropType<TypeTag, Properties::Scalar>;

using Element = typename GridView::Traits::template Codim<0>::Entity;
using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;
using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
using PrimaryVariables = typename SolutionTypes::PrimaryVariables;

public:
IMPESTestProblem(TimeManager &timeManager, Grid &grid) :
ParentType(timeManager, grid)
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
const std::string name() const
{
    return name_;
}

bool shouldWriteRestartFile() const
{
    return false;
}

/*!
 * \brief Returns the temperature within the domain.
 *
 * This problem assumes a temperature of 10 degrees Celsius.
 */
Scalar temperatureAtPos(const GlobalPosition& globalPos) const
{
    return 273.15 + 10; // -> 10°C
}

// \}

//! Returns the reference pressure for evaluation of constitutive relations
Scalar referencePressureAtPos(const GlobalPosition& globalPos) const
{
    return 1e5; // -> 10°C
}

void source(PrimaryVariables &values,const Element& element) const
{
    values = 0;
}

/*!
* \brief Returns the type of boundary condition.
*
* BC for pressure equation can be dirichlet (pressure) or neumann (flux).
*
* BC for saturation equation can be dirichlet (saturation), neumann (flux), or outflow.
*/
void boundaryTypesAtPos(BoundaryTypes &bcTypes, const GlobalPosition& globalPos) const
{
        if (globalPos[0] < eps_)
        {
            bcTypes.setAllDirichlet();
        }
        else if (globalPos[0] > this->bBoxMax()[0] - eps_)
        {
            bcTypes.setNeumann(eqIdxPress);
            bcTypes.setOutflow(eqIdxSat);
        }
        // all other boundaries
        else
        {
            bcTypes.setAllNeumann();
        }
}

//! set dirichlet condition  (pressure [Pa], saturation [-])
void dirichletAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
{
    values = 0;
    if (globalPos[0] < eps_)
    {
        if (getParam<bool>("Problem.EnableGravity"))
        {
            Scalar pRef = referencePressureAtPos(globalPos);
            Scalar temp = temperatureAtPos(globalPos);

            values[pwIdx] = (2e5 + (this->bBoxMax()[dim-1] - globalPos[dim-1])
                                   * WettingPhase::density(temp, pRef)
                                   * this->gravity().two_norm());
        }
        else
        {
            values[pwIdx] = 2e5;
        }
        values[swIdx] = 0.8;
    }
    else
    {
        values[pwIdx] = 2e5;
        values[swIdx] = 0.2;
    }
}

//! set neumann condition for phases (flux, [kg/(m^2 s)])
void neumannAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
{
    values = 0;
    if (globalPos[0] > this->bBoxMax()[0] - eps_)
    {
        values[nPhaseIdx] = 3e-4;
    }
}
//! return initial solution -> only saturation values have to be given!
void initial(PrimaryVariables &values,
        const Element& element) const
{
    values[pwIdx] = 0;
    values[swIdx] = 0.2;
}

private:

static constexpr Scalar eps_ = 1e-6;
std::string name_;
};
} //end namespace

#endif
