/*****************************************************************************
 *   Copyright (C) 2007-2008 by Markus Wolff, Klaus Mosthaf                  *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
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
#ifndef DUMUX_BUCKLEYLEVERETTPROBLEM_HH
#define DUMUX_BUCKLEYLEVERETTPROBLEM_HH

//#include <fstream>

#include <dune/grid/sgrid.hh>

#include <dumux/material/fluidsystems/liquidphase.hh>
//#include <dumux/material/components/simpleh2o.hh>
//#include <dumux/material/components/oil.hh>
#include <pseudooil.hh>
#include <pseudoh2o.hh>

#include <dumux/decoupled/2p/impes/impesproblem2p.hh>
#include <dumux/decoupled/2p/diffusion/fv/fvvelocity2p.hh>
#include <dumux/decoupled/2p/transport/fv/fvsaturation2p.hh>

#include "buckleyleverett_spatialparams.hh"
#include "buckleyleverett_analytic.hh"

namespace Dumux
{
template <class TypeTag>
class BuckleyLeverettProblem;

//////////
// Specify the properties
//////////
namespace Properties
{
NEW_TYPE_TAG(BuckleyLeverettProblem, INHERITS_FROM(DecoupledTwoP, Transport, BuckleyLeverettSpatialParams));

// Set the grid type
SET_PROP(BuckleyLeverettProblem, Grid)
{
    typedef Dune::SGrid<2, 2> type;
};

// Set the problem property
SET_PROP(BuckleyLeverettProblem, Problem)
{
public:
    typedef Dumux::BuckleyLeverettProblem<TypeTag> type;
};

// Set the model properties
SET_PROP(BuckleyLeverettProblem, TransportModel)
{
    typedef Dumux::FVSaturation2P<TTAG(BuckleyLeverettProblem)> type;
};

SET_PROP(BuckleyLeverettProblem, PressureModel)
{
    typedef Dumux::FVVelocity2P<TTAG(BuckleyLeverettProblem)> type;
};

// Set the wetting phase
SET_PROP(BuckleyLeverettProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::PseudoOil<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(BuckleyLeverettProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::PseudoH2O<Scalar> > type;
};

// Enable gravity
SET_BOOL_PROP(BuckleyLeverettProblem, EnableGravity, false);

SET_SCALAR_PROP(BuckleyLeverettProblem, CFLFactor, 0.95); //0.95
}

/*!
 * \ingroup DecoupledProblems
 */
template<class TypeTag = TTAG(BuckleyLeverettProblem)>
class BuckleyLeverettProblem: public IMPESProblem2P<TypeTag>
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
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
        pWIdx = Indices::pwIdx,
        SwIdx = Indices::SwIdx,
        pressEqIdx = Indices::pressEqIdx,
        satEqIdx = Indices::satEqIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef typename GET_PROP(TypeTag, PTAG(ParameterTree)) Params;

public:
    BuckleyLeverettProblem(TimeManager& timeManager, const GridView &gridView,
                           const GlobalPosition lowerLeft = 0,
                           const GlobalPosition upperRight = 0,
                           const Scalar pleftbc = 2e5)
    : ParentType(timeManager, gridView),
      lowerLeft_(lowerLeft),
      upperRight_(upperRight),
      eps_(1e-8),
      pLeftBc_(pleftbc),
      analyticSolution_(*this, 3e-7)
    {
        densityNonWetting_ = Params::tree().template get<double>("Fluid.densityNW");

        //Write header for ViPLab-Outputfile
    	std::ofstream dataFile;
		dataFile.open("vipplot.vgf");
		dataFile.close();
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
    const char *name() const
    {
        std::string simName = "buckleyleverett_linear_run";

        return (str(boost::format("%s-")
                %simName).c_str());
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

    // \}

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
            if (globalPos[0] < eps_)
            {
                bcTypes.setAllDirichlet();
            }
            else if (globalPos[0] > upperRight_[0] - eps_)
            {
                bcTypes.setNeumann(pressEqIdx);
                bcTypes.setOutflow(satEqIdx);
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
            values[pWIdx] = pLeftBc_;
            values[SwIdx] = 0.8;
        }
        else
        {
            values[pWIdx] = pLeftBc_;
            values[SwIdx] = 0.2;
        }
    }

    void neumannAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
    {
        values = 0;
        if (globalPos[0]> upperRight_[0] - eps_)
        {
            // the volume flux should remain constant, when density is changed
            // here, we multiply by the density of the NonWetting Phase
            const Scalar referenceDensity = 1000.0;

            values[nPhaseIdx] = 3e-4 * densityNonWetting_/referenceDensity;
        }
    }

    void initialAtPos(PrimaryVariables &values,
            const GlobalPosition &globalPos) const
    {
        values[pWIdx] = 0;
        values[SwIdx] = 0.2;
    }


    //Override outputfunction for ViPLab-Output
    void writeOutput()
    {
    	double time = this->timeManager().time();
    	if (time < 0)
    		return;

        double discretizationLength = Params::tree().template get<double>("Geometry.discretizationLength");
		int cellNumberX = static_cast<int>(300/discretizationLength);
		discretizationLength = 300/cellNumberX; //Might be slightly different from the input parameter

    	std::ofstream dataFile;
		dataFile.open("vipplot.vgf", std::fstream::app);

    	std::cout<<"Writing output for time step.\n";

    	if (time  > 0)
    		dataFile << "# newframe\n";

		dataFile << "# title Buckley-Leverett-Problem\n";
		dataFile << "# x-label x\n";
		dataFile << "# y-label Saturation\n";
		dataFile << "# x-range -10 310\n";
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
//		ElementIterator eItEnd = this->gridView().template end<0>();
//	    for (ElementIterator eIt = this->gridView().template begin<0>(); eIt != eItEnd; ++eIt)
//	    {
//	    	int currentIdx = this->variables().index(*eIt);
//	    	Scalar sat = this->variables().saturation()[currentIdx];
//	    	Scalar leftPos = eIt->geometry().corner(0)[0];
//	    	Scalar rightPos = eIt->geometry().corner(1)[0];
//	    	dataFile << leftPos <<" "<<sat<<" "<<rightPos<<" "<<sat<<"\n";
//	    }


//		for (ElementIterator eIt = this->gridView().template begin<0>(); eIt != eItEnd; ++eIt)
//	    {
//	    	int currentIdx = this->variables().index(*eIt);
//	    	Scalar sat = analyticSolution_.AnalyticSolution()[currentIdx];
//	    	Scalar leftPos = eIt->geometry().corner(0)[0];
//	    	Scalar rightPos = eIt->geometry().corner(1)[0];
//	    	dataFile << leftPos <<" "<<sat<<" "<<rightPos<<" "<<sat<<"\n";
//	    }

		dataFile << "# legend Numerical Solution,Analytic Solution\n";
		dataFile.close();
    }

private:
    GlobalPosition lowerLeft_;
    GlobalPosition upperRight_;
    Scalar eps_;
    Scalar pLeftBc_;
    Scalar densityNonWetting_;
    BuckleyLeverettAnalytic<TypeTag> analyticSolution_;
};
}
#endif
