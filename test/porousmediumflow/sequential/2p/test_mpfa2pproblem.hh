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
 * \brief test problem for sequential 2p models
 */

#ifndef DUMUX_TEST_MPFA2P_PROBLEM_HH
#define DUMUX_TEST_MPFA2P_PROBLEM_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/properties.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/trichloroethene.hh>

#include <dumux/porousmediumflow/2p/sequential/diffusion/cellcentered/pressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/cellcentered/pressurepropertiesadaptive.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/mpfa/omethod/2dpressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/mpfa/lmethod/2dpressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/mpfa/lmethod/2dpressurepropertiesadaptive.hh>
#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/properties.hh>
#include <dumux/porousmediumflow/2p/sequential/impes/problem.hh>


#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/evalcflfluxcoats.hh>
#include <dumux/porousmediumflow/2p/sequential/impes/gridadaptionindicatorlocal.hh>

#include "test_mpfa2pspatialparams.hh"
#include "buckleyleverettanalyticsolution.hh"
#include "mcwhorteranalyticsolution.hh"

namespace Dumux
{

template<class TypeTag>
class MPFATwoPTestProblem;

//////////
// Specify the properties
//////////
namespace Properties
{

// Create new type tags
namespace TTag {
struct MPFATwoPTest { using InheritsFrom = std::tuple<Test2PSpatialParams>; };
struct FVTwoPTest { using InheritsFrom = std::tuple<MPFATwoPTest, IMPESTwoP, FVTransportTwoP, FVPressureTwoP>; };
struct FVAdaptiveTwoPTest { using InheritsFrom = std::tuple<MPFATwoPTest, IMPESTwoPAdaptive, FVTransportTwoP, FVPressureTwoPAdaptive>; };
struct MPFAOTwoPTest { using InheritsFrom = std::tuple<MPFATwoPTest, IMPESTwoP, FVTransportTwoP, FvMpfaO2dPressureTwoP>; };
struct MPFALTwoPTest { using InheritsFrom = std::tuple<MPFATwoPTest, IMPESTwoP, FVTransportTwoP, FvMpfaL2dPressureTwoP>; };
struct MPFALAdaptiveTwoPTest { using InheritsFrom = std::tuple<MPFATwoPTest, IMPESTwoPAdaptive, FVTransportTwoP, FvMpfaL2dPressureTwoPAdaptive>; };

} // end namespace TTag

// Set the grid type
#if HAVE_UG
template<class TypeTag>
struct Grid<TypeTag, TTag::MPFATwoPTest> { using type = Dune::UGGrid<2>; };
#else
template<class TypeTag>
struct Grid<TypeTag, TTag::MPFATwoPTest> { using type = Dune::YaspGrid<2>; };
#endif

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::MPFATwoPTest> { using type = MPFATwoPTestProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::MPFATwoPTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
#if PROBLEM == 2
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Trichloroethene<Scalar> >;
#else
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
#endif
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

#if PROBLEM == 1
template<class TypeTag>
struct Formulation<TypeTag, TTag::MPFATwoPTest> { static constexpr int value = SequentialTwoPCommonIndices::pnsw; };
#endif

template<class TypeTag>
struct EvalCflFluxFunction<TypeTag, TTag::MPFATwoPTest> { using type = EvalCflFluxCoats<TypeTag>; };
template<class TypeTag>
struct AdaptionIndicator<TypeTag, TTag::MPFATwoPTest> { using type = GridAdaptionIndicator2PLocal<TypeTag>; };
}

/*!
 * \ingroup SequentialTwoPTests
 * \brief test problem for sequential 2p models
 *
 * Trichloroethene (a DNAPL) is injected from the top into a rectangular 2D domain saturated by water.
 * The remaining upper and the lower boundary is closed (Neumann = 0). At the sides a hydrostatic pressure condition
 * and free outflow for saturation are set. The domain is heterogeneous with a backround material and three lenses.
 *
 * To run the simulation execute the following line in shell:
 *
 * <tt>./test_mpfa2p</tt>,
 *
 * Additionally, the numerical model can be switched by executing with the parameter "ModelType":
 *
 * <tt>./test_mpfa2p --ModelType=...</tt>,
 *
 * where ModelType can be:
 * - FV (standard finite volume)
 * - FVAdaptive (adaptive finite volume)
 * - MPFAO (MPFA o-method)
 * - MPFAL (MPFA l-method)
 * - MPFALAdaptive (adaptive MPFA l-method)
 */
template<class TypeTag>
class MPFATwoPTestProblem: public IMPESProblem2P<TypeTag>
{
using ParentType = IMPESProblem2P<TypeTag>;
using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
using Grid = typename GridView::Grid;

using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

using WettingPhase = typename GetProp<TypeTag, Properties::FluidSystem>::WettingPhase;

using TimeManager = GetPropType<TypeTag, Properties::TimeManager>;

using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;
using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
using PrimaryVariables = typename SolutionTypes::PrimaryVariables;

enum
{
    dim = GridView::dimension, dimWorld = GridView::dimensionworld
};

enum
{
    nPhaseIdx = Indices::nPhaseIdx,
#if PROBLEM == 1
    pnIdx = Indices::pnIdx,
#else
    pwIdx = Indices::pwIdx,
#endif
    swIdx = Indices::swIdx,
    eqIdxPress = Indices::pressureEqIdx,
    eqIdxSat = Indices::satEqIdx
};

using Scalar = GetPropType<TypeTag, Properties::Scalar>;

using Element = typename GridView::Traits::template Codim<0>::Entity;
using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
MPFATwoPTestProblem(TimeManager &timeManager, Grid &grid) :
ParentType(timeManager, grid)
#if PROBLEM != 2
, analyticSolution_(*this)
#endif
{
    int refinementFactor = getParam<Scalar>("Grid.RefinementFactor", 0);
    this->grid().globalRefine(refinementFactor);

    Scalar inletWidth = getParam<Scalar>("Problem.InletWidth", 1.0);
    GlobalPosition inletCenter = this->bBoxMax();
    inletCenter[0] *= 0.5;

    inletLeftCoord_ = inletCenter;
    inletLeftCoord_[0] -=0.5*inletWidth;
    inletRightCoord_ = inletCenter;
    inletRightCoord_[0] +=0.5*inletWidth;

    inFlux_ = getParam<Scalar>("Problem.InjectionFlux", 1e-4);

    int outputInterval = getParam<int>("Problem.OutputInterval", 0);
    this->setOutputInterval(outputInterval);

    Scalar outputTimeInterval = getParam<Scalar>("Problem.OutputTimeInterval", 1e6);
    this->setOutputTimeInterval(outputTimeInterval);
}

#if PROBLEM != 2
void init()
{
    ParentType::init();

#if PROBLEM == 0
    Scalar vTot = 3e-6;
    analyticSolution_.initialize(vTot);
#endif
#if PROBLEM == 1
    analyticSolution_.initialize();
#endif
}

void addOutputVtkFields()
{
    analyticSolution_.calculateAnalyticSolution();

    ParentType::addOutputVtkFields();
    analyticSolution_.addOutputVtkFields(this->resultWriter());
}
#endif

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
    return getParam<std::string>("Problem.Name");
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
#if PROBLEM == 2
    if (isInlet(globalPos))
    {
        bcTypes.setNeumann(eqIdxPress);
        bcTypes.setDirichlet(swIdx);
    }
    else if (isBottom(globalPos) || isTop(globalPos))
    {
        bcTypes.setAllNeumann();
    }
    else
    {
        bcTypes.setDirichlet(pwIdx);
        bcTypes.setOutflow(eqIdxSat);
    }
#elif  PROBLEM == 0
    if (globalPos[0] < eps_)
    {
        bcTypes.setAllDirichlet();
    }
    else if (globalPos[0] > this->bBoxMax()[0] - eps_)
    {
        bcTypes.setNeumann(eqIdxPress);
        bcTypes.setOutflow(eqIdxSat);
    }
    else
    {
        bcTypes.setAllNeumann();
    }
#else
    if (globalPos[0] < eps_)
    {
        bcTypes.setAllDirichlet();
    }
    else
    {
        bcTypes.setAllNeumann();
    }
#endif
}

//! set dirichlet condition  (pressure [Pa], saturation [-])
void dirichletAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
{
    values = 0;
#if PROBLEM == 2
    Scalar pRef = referencePressureAtPos(globalPos);
    Scalar temp = temperatureAtPos(globalPos);

    values[pwIdx] = (1e5 - (this->bBoxMax()- globalPos) * this->gravity() * WettingPhase::density(temp, pRef));
    values[swIdx] = 1.0;

    if (isInlet(globalPos))
    {
        values[swIdx] = 0.0;
    }
#elif PROBLEM == 0
    if (globalPos[0] < eps_)
    {
        values[swIdx] = 0.8;
        values[pwIdx] = 1;
    }
#else
    if (globalPos[0] < eps_)
    {
        values[swIdx] = 1.0;
        values[pnIdx] = 1e5;
    }
#endif
}

//! set neumann condition for phases (flux, [kg/(m^2 s)])
void neumannAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
{
    values = 0;
#if PROBLEM == 2
    if (isInlet(globalPos))
    {
        values[nPhaseIdx] = -inFlux_;
    }
#elif PROBLEM == 0
    if (globalPos[0] > this->bBoxMax()[0] - eps_)
    {
        values[nPhaseIdx] = 3e-3;
    }
#endif
}

//! return initial solution -> only saturation values have to be given!
void initialAtPos(PrimaryVariables &values,
        const GlobalPosition& globalPos) const
{
#if PROBLEM == 2
    values[pwIdx] = 0;
    values[swIdx] = 1.0;
#elif PROBLEM == 0
    values[pwIdx] = 0;
    values[swIdx] = 0.2;
#else
    values[pnIdx] = 0;
    values[swIdx] = 0.0;
#endif
}

private:

bool isInlet(const GlobalPosition& globalPos) const
{
        if (!isTop(globalPos))
            return false;

    for (int i = 0; i < dim; i++)
    {
        if (globalPos[i] < inletLeftCoord_[i] - eps_)
            return false;
        if (globalPos[i] > inletRightCoord_[i] + eps_)
            return false;
    }
    return true;
}

bool isTop(const GlobalPosition& globalPos) const
{
    if (dim == 2)
    {
        if (globalPos[1] > this->bBoxMax()[1] - eps_)
            return true;
    }
    if (dim == 3)
    {
        if (globalPos[2] > this->bBoxMax()[2] - eps_)
            return true;
    }
    return false;
}

bool isBottom(const GlobalPosition& globalPos) const
{
    if (dim == 2)
    {
        if (globalPos[1] < eps_)
            return true;
    }
    if (dim == 3)
    {
        if (globalPos[2] < eps_)
            return true;
    }
    return false;
}

Scalar inFlux_;
GlobalPosition inletLeftCoord_;
GlobalPosition inletRightCoord_;
static constexpr Scalar eps_ = 1e-6;
#if PROBLEM == 0
BuckleyLeverettAnalytic<TypeTag> analyticSolution_;
#endif
#if PROBLEM == 1
McWhorterAnalytic<TypeTag> analyticSolution_;
#endif
};
} //end namespace

#endif
