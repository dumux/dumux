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
 * \brief test problem for sequential 2p models in 3d
 */

#ifndef DUMUX_TEST_3D2P_PROBLEM_HH
#define DUMUX_TEST_3D2P_PROBLEM_HH

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dumux/common/properties.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/components/simpleh2o.hh>

#include <dumux/porousmediumflow/2p/sequential/diffusion/mpfa/lmethod/3dpressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/mpfa/lmethod/3dpressurepropertiesadaptive.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/cellcentered/pressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/cellcentered/pressurepropertiesadaptive.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/mimetic/pressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/mimetic/pressurepropertiesadaptive.hh>
#include <dumux/porousmediumflow/2p/sequential/impes/problem.hh>
#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/properties.hh>

#include "test_3d2pspatialparams.hh"

#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/evalcflfluxcoats.hh>
#include <dumux/porousmediumflow/2p/sequential/impes/gridadaptionindicatorlocal.hh>


namespace Dumux
{
/*!
 * \ingroup SequentialTwoPTests
 * \brief test problem for sequential 2p models in 3d
 */
template<class TypeTag>
class Test3D2PProblem;

//////////
// Specify the properties
//////////
namespace Properties
{
// Create new type tags
namespace TTag {
struct ThreeDTwoPTest { using InheritsFrom = std::tuple<Test3d2pSpatialParams>; };
struct FVTwoPTest { using InheritsFrom = std::tuple<ThreeDTwoPTest, IMPESTwoP, FVTransportTwoP, FVPressureTwoP>; };
struct FVAdaptiveTwoPTest { using InheritsFrom = std::tuple<ThreeDTwoPTest, IMPESTwoPAdaptive, FVTransportTwoP, FVPressureTwoPAdaptive>; };
struct MPFALTwoPTest { using InheritsFrom = std::tuple<ThreeDTwoPTest, IMPESTwoP, FVTransportTwoP, FvMpfaL3dPressureTwoP>; };
struct MPFALAdaptiveTwoPTest { using InheritsFrom = std::tuple<ThreeDTwoPTest, IMPESTwoPAdaptive, FVTransportTwoP, FvMpfaL3dPressureTwoPAdaptive>; };
struct MimeticTwoPTest { using InheritsFrom = std::tuple<ThreeDTwoPTest, IMPESTwoP, FVTransportTwoP, MimeticPressureTwoP>; };
struct MimeticAdaptiveTwoPTest { using InheritsFrom = std::tuple<ThreeDTwoPTest, IMPESTwoPAdaptive, FVTransportTwoP, MimeticPressureTwoPAdaptive>; };
} // end namespace TTag

// Set the grid type
#if HAVE_DUNE_ALUGRID
template<class TypeTag>
struct Grid<TypeTag, TTag::ThreeDTwoPTest> { using type = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>; };
#endif

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::ThreeDTwoPTest> { using type = Test3D2PProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::ThreeDTwoPTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

#if PROBLEM == 1
template<class TypeTag>
struct Formulation<TypeTag, TTag::ThreeDTwoPTest> { static constexpr int value = SequentialTwoPCommonIndices::pnSw; };
template<class TypeTag>
struct EvalCflFluxFunction<TypeTag, TTag::ThreeDTwoPTest> { using type = EvalCflFluxCoats<TypeTag>; };
#endif

template<class TypeTag>
struct AdaptionIndicator<TypeTag, TTag::ThreeDTwoPTest> { using type = GridAdaptionIndicator2PLocal<TypeTag>; };
}

/*!
 * \ingroup IMPETtests
 *
 * \brief test problem for the sequential 2p models in 3d
 *
 * Standard setting: one eighth of a nine-spot waterflood problem. The domain of size 1 x 1 x 1 m is initially saturated by a nonwetting fluid.
 * Water is injected in the corner of the origin (0,0,0) and nonwetting fluid produced in the upper corner at (1,1,1).
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_3d2p</tt>.
 */
template<class TypeTag>
class Test3D2PProblem: public IMPESProblem2P<TypeTag>
{
using ThisType = Test3D2PProblem<TypeTag>;
using ParentType = IMPESProblem2P<TypeTag>;
using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
using Grid = typename GridView::Grid;

using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
using FluidState = GetPropType<TypeTag, Properties::FluidState>;

using TimeManager = GetPropType<TypeTag, Properties::TimeManager>;

using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;
using PrimaryVariables = typename GetProp<TypeTag, Properties::SolutionTypes>::PrimaryVariables;

enum
{
    dim = GridView::dimension, dimWorld = GridView::dimensionworld
};

enum
{
    wPhaseIdx = Indices::wPhaseIdx,
    nPhaseIdx = Indices::nPhaseIdx,
#if PROBLEM == 1
    pNIdx = Indices::pnIdx,
#else
    pWIdx = Indices::pwIdx,
#endif
    swIdx = Indices::swIdx,
    pressureEqIdx = Indices::pressureEqIdx,
    satEqIdx = Indices::satEqIdx
};

using Scalar = GetPropType<TypeTag, Properties::Scalar>;

using Element = typename GridView::Traits::template Codim<0>::Entity;
using Intersection = typename GridView::Intersection;
using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
using LocalPosition = Dune::FieldVector<Scalar, dim>;

public:

Test3D2PProblem(TimeManager& timeManager, Grid& grid) :
ParentType(timeManager, grid), inflowEdge_(0), outflowEdge_(0)
{
    int refinementFactor = 0;
    if (hasParam("Grid.RefinementFactor") && !getPropValue<TypeTag, Properties::AdaptiveGrid>())
    {
        refinementFactor = getParam<Scalar>("Grid.RefinementFactor");
        grid.globalRefine(refinementFactor);
    }

    Scalar minDist = this->bBoxMax().two_norm();
    Scalar maxDist = this->bBoxMin().two_norm();

    // calculate the bounding box of the grid view
    for (const auto& vertex : vertices(grid.leafGridView())) {
        GlobalPosition vertexCoord(vertex.geometry().center());

        Scalar dist = vertexCoord.two_norm();
        if (dist > maxDist)
        {
            maxDist = dist;
            outflowEdge_ = vertex.geometry().center();
        }
        if (dist < minDist)
        {
            minDist = dist;
            inflowEdge_ = vertex.geometry().center();
        }
    }

    int outputInterval = 0;
    if (hasParam("Problem.OutputInterval"))
    {
        outputInterval = getParam<int>("Problem.OutputInterval");
    }
    this->setOutputInterval(outputInterval);

    Scalar outputTimeInterval = 1e6;
    if (hasParam("Problem.OutputTimeInterval"))
    {
        outputTimeInterval = getParam<Scalar>("Problem.OutputTimeInterval");
    }
    this->setOutputTimeInterval(outputTimeInterval);
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
    return getParam<std::string>("Problem.OutputName", "test_3d2p");
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

Scalar referencePressureAtPos(const GlobalPosition& globalPos) const
{
    return 1e5; // -> 10°C
}

void sourceAtPos(PrimaryVariables &values,const GlobalPosition& globalPos) const
{
    values = 0;

#if PROBLEM == 2 //Nine-Spot
        if (globalPos[0] < 0.1 * outflowEdge_[0] + eps_ && globalPos[1] < 0.1 * outflowEdge_[1] && globalPos[2] < 0.1 * outflowEdge_[2])
        {
            values[wPhaseIdx] = 1;
        }
#endif
}

void source(PrimaryVariables &values,const Element& element) const
{
    values = 0;

#if PROBLEM == 2 //Nine-Spot
    int numVertices = element.geometry().corners();
    for(int i = 0; i < numVertices; i++)
    {
        GlobalPosition globalPos(element.template subEntity<dim>(i).geometry().center());

        if (globalPos[0] < inflowEdge_[0] + eps_ && globalPos[1] < inflowEdge_[1] + eps_ && globalPos[2] < inflowEdge_[2] + eps_)
        {
            values[wPhaseIdx] = 1;
            break;
        }
    }
#endif
}

void boundaryTypesAtPos(BoundaryTypes &bcTypes, const GlobalPosition& globalPos) const
{
#if PROBLEM == 2 //Nine-Spot
        if (globalPos[0] > 0.9 * outflowEdge_[0]  && globalPos[1] > 0.9 * outflowEdge_[1] && globalPos[2] > 0.9 * outflowEdge_[2])
        {
            bcTypes.setAllDirichlet();
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
        else if (globalPos[0] > this->bBoxMax()[0] - eps_)
        {
            bcTypes.setNeumann(pressureEqIdx);
            bcTypes.setOutflow(satEqIdx);
        }
        // all other boundaries
        else
        {
            bcTypes.setAllNeumann();
        }
#endif
}

void boundaryTypes(BoundaryTypes &bcTypes, const Intersection& intersection) const
{
#if PROBLEM == 2 //Nine-Spot
    auto element = intersection.inside();
    int numVertices = element.geometry().corners();
    for(int i = 0; i < numVertices; i++)
    {
        GlobalPosition globalPos(element.template subEntity<dim>(i).geometry().center());

        if (globalPos[0] > outflowEdge_[0] - eps_ && globalPos[1] > outflowEdge_[1] - eps_ && globalPos[2] > outflowEdge_[2] - eps_)
        {
            bcTypes.setAllDirichlet();
            break;
        }
        else
        {
            bcTypes.setAllNeumann();
        }
    }
#else
    GlobalPosition globalPos(intersection.geometry().center());

    boundaryTypesAtPos(bcTypes,globalPos);
#endif
}

void dirichletAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
{
    values = 0;

#if PROBLEM == 2 //Nine-Spot
            values[swIdx] = 0.8;
            values[pWIdx] = 1;
#elif PROBLEM == 0
    if (globalPos[0] < eps_)
    {
        values[swIdx] = 0.8;
        values[pWIdx] = 1;
    }
#elif PROBLEM == 1
    if (globalPos[0] < eps_)
    {
        values[swIdx] = 1.0;
        values[pNIdx] = 1e5;
    }
#endif
}

void neumannAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
{
    values = 0;

#if PROBLEM == 0 //Buckley-Leverett
    if (globalPos[0] > this->bBoxMax()[0] - eps_)
    {
        values[nPhaseIdx] = 3e-3;
    }
#endif
}

void initialAtPos(PrimaryVariables &values,
        const GlobalPosition& globalPos) const
{
#if PROBLEM == 1
    values[pNIdx] = 0;
    values[swIdx] = 0.0;
#else
    values[pWIdx] = 0;
    values[swIdx] = 0.2;
#endif
}

private:
static constexpr Scalar eps_ = 1e-6;
GlobalPosition inflowEdge_;
GlobalPosition outflowEdge_;
};
} //end namespace

#endif
