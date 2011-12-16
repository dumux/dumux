// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
NEW_TYPE_TAG(McWhorterProblem, INHERITS_FROM(DecoupledTwoP, Transport, McWhorterSpatialParams));

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

SET_PROP(McWhorterProblem, PressureModel)
{
    typedef Dumux::FVVelocity2P<TTAG(McWhorterProblem)> type;
};

SET_INT_PROP(McWhorterProblem, Formulation,
        DecoupledTwoPCommonIndices::pnSw);

// Set the wetting phase
SET_PROP(McWhorterProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::PseudoOil<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(McWhorterProblem, NonWettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::PseudoH2O<Scalar> > type;
};

// Disable gravity
SET_BOOL_PROP(McWhorterProblem, EnableGravity, false);

SET_TYPE_PROP(McWhorterProblem, EvalCflFluxFunction, Dumux::EvalCflFluxCoats<TypeTag>);
SET_SCALAR_PROP(McWhorterProblem, CFLFactor, 0.8);
}

//! \ingroup transportProblems
//! @brief McWhorter transport problem

template<class TypeTag = TTAG(McWhorterProblem)>
class McWhorterProblem: public IMPESProblem2P<TypeTag>
{
    typedef IMPESProblem2P<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;
    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes))::PrimaryVariables PrimaryVariables;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        pNIdx = Indices::pnIdx,
        SwIdx = Indices::SwIdx,
        pressEqIdx = Indices::pressEqIdx,
        satEqIdx = Indices::satEqIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef typename GET_PROP(TypeTag, PTAG(ParameterTree)) Params;

public:

    McWhorterProblem(TimeManager& timeManager, const GridView &gridView,
                            const GlobalPosition lowerLeft = 0,
                            const GlobalPosition upperRight = 0,
                            const Scalar pleftbc = 2.0e5)//? initial oil pressure?
     : ParentType(timeManager, gridView),
       lowerLeft_(lowerLeft),
       upperRight_(upperRight),
       eps_(1e-6),
       pLeftBc_(pleftbc),
       analyticSolution_(*this)
     {
        this->setOutputInterval(10);

        //Write header for ViPLab-Outputfile
        std::ofstream dataFile;
        dataFile.open("vipplot.vgf");
        dataFile.close();
     }

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

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    {
        return 273.15 + 10; // -> 10°C
    }

    Scalar referencePressureAtPos(const GlobalPosition& globalPos) const
    {
        return 1e5; // -> 10°C
    }

    void sourceAtPos(PrimaryVariables &values,const GlobalPosition& globalPos) const
    {
        values = 0;
    }

    void boundaryTypesAtPos(BoundaryTypes &bcTypes, const GlobalPosition& globalPos) const
    {
            if (globalPos[0] < eps_)//west
            {
                bcTypes.setAllDirichlet();
            }
            // all other boundaries
            else
            {
                bcTypes.setAllNeumann();
            }
    }

    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
    {
        if (globalPos[0] < eps_)
        {
            values[pNIdx] = pLeftBc_;
            values[SwIdx] = 1.0;
        }
        else
        {
            values[pNIdx] = pLeftBc_;
            values[SwIdx] = 0.0;
        }
    }

    void neumannAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
    {
        values = 0;
    }

    void initialAtPos(PrimaryVariables &values,
            const GlobalPosition &globalPos) const
    {
        values = 0;
    }

    //Override outputfunction for ViPLab-Output
    void writeOutput()
    {
        double time = this->timeManager().time();
        if (time < 0)
            return;

        double discretizationLength = Params::tree().template get<double>("problem.DiscretizationLength");
        int cellNumberX = static_cast<int>(2/discretizationLength);
        discretizationLength = 2.0/cellNumberX; //Might be slightly different from the input parameter

        std::ofstream dataFile;
        dataFile.open("vipplot.vgf", std::fstream::app);

        std::cout<<"Writing output for time step.\n";

        if (time  > 0)
            dataFile << "# newframe\n";

        dataFile << "# title Mc Whorther Problem\n";
        dataFile << "# x-label x\n";
        dataFile << "# y-label Saturation\n";
        dataFile << "# x-range -0.1 2.1\n";
        dataFile << "# y-range 0 1\n";

        dataFile << "# color 0 0 255\n";

        Scalar curSat = this->variables().saturation()[0];

        // Draw first piece separately, as no vertical line is needed
        dataFile << 0 << " " << curSat << " "
                 << discretizationLength << " " << curSat <<std::endl;
        for (int i=1; i < cellNumberX; i++)
        {
            // Vertical Line
            dataFile << i*discretizationLength << " "
                     << curSat << " ";
            curSat = this->variables().saturation()[i];
            dataFile << i*discretizationLength << " " << curSat << std::endl;
            //Horizontal Line
            dataFile << i*discretizationLength << " " << curSat << " "
                     << (i+1)*discretizationLength << " " << curSat << std::endl;
        }
        dataFile << "# color 255 0 0\n";

        curSat = analyticSolution_.AnalyticSolution()[0];
        dataFile << 0 << " " << curSat << " "
                 << discretizationLength << " " << curSat <<std::endl;
        for (int i=1; i < cellNumberX; i++)
        {
            // Vertical Line
            dataFile << i*discretizationLength << " "
                     << curSat << " ";
            curSat = analyticSolution_.AnalyticSolution()[i];
            dataFile << i*discretizationLength << " " << curSat << std::endl;
            //Horizontal Line
            dataFile << i*discretizationLength << " " << curSat << " "
                     << (i+1)*discretizationLength << " " << curSat << std::endl;
        }
        dataFile << "# legend Numerical Solution,Analytic Solution\n";
        dataFile.close();
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
