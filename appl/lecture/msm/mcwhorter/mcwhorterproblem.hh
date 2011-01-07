/*****************************************************************************
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#ifndef DUMUX_MCWHORTERPROBLEM_HH
#define DUMUX_MCWHORTERPROBLEM_HH

#include <dune/grid/sgrid.hh>

#include <dumux/material/fluidsystems/liquidphase.hh>
//#include <dumux/material/components/simpleh2o.hh>
//#include <dumux/material/components/oil.hh>

#include "../buckleyleverett/pseudooil.hh"
#include "../buckleyleverett/pseudoh2o.hh"

#include <dumux/decoupled/2p/impes/impesproblem2p.hh>
#include <dumux/decoupled/2p/diffusion/fv/fvvelocity2p.hh>
#include <dumux/decoupled/2p/transport/fv/fvsaturation2p.hh>
#include <dumux/decoupled/2p/transport/fv/capillarydiffusion.hh>
#include <dumux/decoupled/2p/transport/fv/gravitypart.hh>
#include<dumux/decoupled/2p/transport/fv/evalcflflux_coats.hh>

#include "mcwhorter_spatialparams.hh"
#include "mcwhorter_analytic.hh"

namespace Dumux
{
template <class TypeTag>
class McWhorterProblem;

//////////
// Specify the properties
//////////
namespace Properties
{
NEW_TYPE_TAG(McWhorterProblem, INHERITS_FROM(DecoupledTwoP, Transport));

// Set the grid type
SET_PROP(McWhorterProblem, Grid)
{
    //    typedef Dune::YaspGrid<2> type;
    typedef Dune::SGrid<1, 1> type;
};

// Set the problem property
SET_PROP(McWhorterProblem, Problem)
{
public:
    typedef Dumux::McWhorterProblem<TypeTag> type;
};

// Set the model properties
SET_PROP(McWhorterProblem, TransportModel)
{
    typedef Dumux::FVSaturation2P<TTAG(McWhorterProblem)> type;
};
SET_TYPE_PROP(McWhorterProblem, DiffusivePart, Dumux::CapillaryDiffusion<TypeTag>);
SET_TYPE_PROP(McWhorterProblem, ConvectivePart, Dumux::GravityPart<TypeTag>);

SET_PROP(McWhorterProblem, PressureModel)
{
    typedef Dumux::FVVelocity2P<TTAG(McWhorterProblem)> type;
};

//SET_INT_PROP(McWhorterProblem, VelocityFormulation,
//        GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices))::velocityW);

SET_INT_PROP(McWhorterProblem, PressureFormulation,
        GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices))::pressureNW);

// Set the wetting phase
SET_PROP(McWhorterProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::PseudoOil<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(McWhorterProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::PseudoH2O<Scalar> > type;
};

// Set the spatial parameters
SET_PROP(McWhorterProblem, SpatialParameters)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

public:
    typedef Dumux::McWhorterSpatialParams<TypeTag> type;
};

// Disable gravity
SET_BOOL_PROP(McWhorterProblem, EnableGravity, false);

SET_TYPE_PROP(McWhorterProblem, EvalCflFluxFunction, Dumux::EvalCflFluxCoats<TypeTag>);
SET_SCALAR_PROP(McWhorterProblem, CFLFactor, 0.8);
}

//! \ingroup transportProblems
//! @brief McWhorter transport problem

template<class TypeTag = TTAG(McWhorterProblem)>
class McWhorterProblem: public IMPESProblem2P<TypeTag, McWhorterProblem<TypeTag> >
{
    typedef McWhorterProblem<TypeTag> ThisType;
    typedef IMPESProblem2P<TypeTag, ThisType> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
    };
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;

public:

    McWhorterProblem(const GridView &gridView,
                            const GlobalPosition lowerLeft = 0,
                            const GlobalPosition upperRight = 0,
                            const Scalar pleftbc = 2.0e5)//? initial oil pressure?
     : ParentType(gridView),
       lowerLeft_(lowerLeft),
       upperRight_(upperRight),
       eps_(1e-6),
       pLeftBc_(pleftbc),
       analyticSolution_(*this)
     {}

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const char *name() const
    {
        return "McWhorter";
    }

    bool shouldWriteRestartFile() const
    {
        return false;
    }

    void postTimeStep()
    {
        analyticSolution_.calculateAnalyticSolution();

        ParentType::postTimeStep();
    };

    void addOutputVtkFields()
    {
        ParentType::addOutputVtkFields();
        analyticSolution_.addOutputVtkFields(this->resultWriter());
    }

    bool shouldWriteOutput() const
    {
        if (this->timeManager().timeStepIndex() % 10 == 0)
        {
        return true;
        }
        return false;
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */

    Scalar temperature(const GlobalPosition& globalPos, const Element& element) const
    {return 273.15 + 10; // -> 10°C
    }

    Scalar referencePressure(const GlobalPosition& globalPos, const Element& element) const
    {
        return 1e5;
    }
     std::vector<Scalar> source (const GlobalPosition& globalPos, const Element& element)
    {
        return std::vector<Scalar>(2,0.0);//no source and sink terms
    }

    BoundaryConditions::Flags bctypePress(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        if (globalPos[0] < eps_)//west
        return BoundaryConditions::dirichlet;

        // all other boundaries
        else
        return BoundaryConditions::neumann;
    }

    BoundaryConditions::Flags bctypeSat(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        if (globalPos[0] < eps_)// west and east
        return Dumux::BoundaryConditions::dirichlet;
        else
        return Dumux::BoundaryConditions::neumann;
    }

    Scalar dirichletPress(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        return pLeftBc_;
    }

    Scalar dirichletSat(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        if (globalPos[0] < eps_)//west
        return 1.0;
        else //east, north and south are Neumann
        return 0.0;
    }

    std::vector<Scalar> neumann(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        return std::vector<Scalar>(2,0.0);
    }

    Scalar initSat (const GlobalPosition& globalPos, const Element& element) const
    {
        return 0.0;
    }

  /*  McWhorterProblem(VC& variables, Fluid& wettingphase, Fluid& nonwettingphase,
            Matrix2p<Grid, Scalar>& soil, TwoPhaseRelations<Grid, Scalar>& materialLaw = *(new TwoPhaseRelations<Grid,Scalar>&),
            const GlobalPosition Right = 0)
    : FractionalFlowProblem<GridView, Scalar, VC>(variables, wettingphase, nonwettingphase, soil, materialLaw),
    right_(Right[0]), eps_(1e-6)
    {}*/

private:
    GlobalPosition lowerLeft_;
    GlobalPosition upperRight_;
    Scalar eps_;
    Scalar pLeftBc_;
    Scalar right_;
    McWhorterAnalytic<TypeTag> analyticSolution_;
 };
}
#endif
