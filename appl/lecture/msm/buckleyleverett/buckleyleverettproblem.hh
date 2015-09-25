// $Id: test_impes_problem.hh 4209 2010-09-01 12:31:45Z mwolff $
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

#include <dune/grid/sgrid.hh>

#include <dumux/material/fluidsystems/liquidphase.hh>
//#include <dumux/material/components/simpleh2o.hh>
//#include <dumux/material/components/oil.hh>
#include <pseudooil.hh>
#include <pseudoh2o.hh>

#include <dumux/decoupled/2p/impes/impesproblem2p.hh>
#include <dumux/decoupled/2p/diffusion/fv/fvvelocity2p.hh>
#include <dumux/decoupled/2p/diffusion/fvmpfa/fvmpfaovelocity2p.hh>
#include <dumux/decoupled/2p/transport/fv/fvsaturation2p.hh>
#include <dumux/decoupled/2p/transport/fv/capillarydiffusion.hh>
#include <dumux/decoupled/2p/transport/fv/gravitypart.hh>

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
NEW_TYPE_TAG(BuckleyLeverettProblem, INHERITS_FROM(DecoupledTwoP, MPFAProperties, Transport));

// Set the grid type
SET_PROP(BuckleyLeverettProblem, Grid)
{
    //    typedef Dune::YaspGrid<2> type;
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
SET_TYPE_PROP(BuckleyLeverettProblem, DiffusivePart, Dumux::CapillaryDiffusion<TypeTag>);
SET_TYPE_PROP(BuckleyLeverettProblem, ConvectivePart, Dumux::GravityPart<TypeTag>);

SET_PROP(BuckleyLeverettProblem, PressureModel)
{
    typedef Dumux::FVVelocity2P<TTAG(BuckleyLeverettProblem)> type;
//    typedef Dumux::FVMPFAOVelocity2P<TTAG(BuckleyLeverettProblem)> type;
};

//SET_INT_PROP(BuckleyLeverettProblem, VelocityFormulation,
//        GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices))::velocityW);

//SET_INT_PROP(BuckleyLeverettProblem, PressureFormulation,
//        GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices))::pressureGlobal);

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

// Set the spatial parameters
SET_PROP(BuckleyLeverettProblem, SpatialParameters)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

public:
    typedef Dumux::BuckleyLeverettSpatialParams<TypeTag> type;
};

// Enable gravity
SET_BOOL_PROP(BuckleyLeverettProblem, EnableGravity, false);

SET_SCALAR_PROP(BuckleyLeverettProblem, CFLFactor, 0.95); //0.95
}

/*!
 * \ingroup DecoupledProblems
 */
template<class TypeTag = TTAG(BuckleyLeverettProblem)>
class BuckleyLeverettProblem: public IMPESProblem2P<TypeTag, BuckleyLeverettProblem<TypeTag> >
{
    typedef BuckleyLeverettProblem<TypeTag> ThisType;
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
        wetting = 0, nonwetting = 1
    };
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    BuckleyLeverettProblem(const GridView &gridView,
                           const GlobalPosition lowerLeft = 0,
                           const GlobalPosition upperRight = 0,
                           const Scalar pleftbc = 2e5)
    : ParentType(gridView),
      lowerLeft_(lowerLeft),
      upperRight_(upperRight),
      eps_(1e-8),
      pLeftBc_(pleftbc),
      analyticSolution_(*this, 3e-7)
    {
        //load interface-file
        Dumux::InterfaceFluidProperties interfaceFluidProps("interface_BL.xml");

        densityNonWetting_ = interfaceFluidProps.IFP_DensityNonWettingFluid;
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
        Dumux::InterfaceProblemProperties interfaceProbProps("interface_BL.xml");
        Scalar simNum =  interfaceProbProps.IPP_SimulationNumber;

        return (str(boost::format("%s-%02d")
                %simName%simNum).c_str());
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

    //! Write the fields current solution into an VTK output file.
    void writeOutput()
    {
//        if (this->timeManager().time() > 43.1999e6)
//        {
            if (this->gridView().comm().rank() == 0)
                std::cout << "Writing result file for current time step\n";

            this->resultWriter().beginTimestep(this->timeManager().time() + this->timeManager().timeStepSize(),
                                        this->gridView());
            this->asImp_().addOutputVtkFields();
            this->resultWriter().endTimestep();
//        }
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature(const GlobalPosition& globalPos, const Element& element) const
    {
        return 273.15 + 10; // -> 10°C
    }

    // \}


    Scalar referencePressure(const GlobalPosition& globalPos, const Element& element) const
    {
        return 1e5;
    }

    std::vector<Scalar> source(const GlobalPosition& globalPos, const Element& element)
    {
        return std::vector<Scalar>(2, 0.0);
    }

    typename BoundaryConditions::Flags bctypePress(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        if (globalPos[0] < eps_)
            return BoundaryConditions::dirichlet;

        // all other boundaries
        else
            return BoundaryConditions::neumann;
    }

    BoundaryConditions::Flags bctypeSat(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        if (globalPos[0] < eps_)
            return Dumux::BoundaryConditions::dirichlet;
        else if (globalPos[0] > upperRight_[0] - eps_)
            return Dumux::BoundaryConditions::outflow;
        else
            return Dumux::BoundaryConditions::neumann;
    }

    Scalar dirichletPress(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        return pLeftBc_;
    }

    Scalar dirichletSat(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        if (globalPos[0] < eps_)
        return 0.8;

        // all other boundaries

        else
        return 0.2;
    }

    std::vector<Scalar> neumann(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        std::vector<Scalar> neumannFlux(2, 0.0);
        if (globalPos[0]> upperRight_[0] - eps_)
        {
            // the volume flux should remain constant, when density is changed
            // here, we multiply by the density of the NonWetting Phase
            const Scalar referenceDensity = 1000.0;
            neumannFlux[nonwetting] = 3e-4 * densityNonWetting_/referenceDensity;
        }
        return neumannFlux;
    }

    Scalar initSat(const GlobalPosition& globalPos, const Element& element) const
    {
        if (globalPos[0] < eps_)
        return 0.8;

        else
        return 0.2;
    }

private:
    GlobalPosition lowerLeft_;
    GlobalPosition upperRight_;
    Scalar eps_;
    Scalar pLeftBc_;
    Scalar simulationNumber_;
    Scalar densityNonWetting_;
    BuckleyLeverettAnalytic<TypeTag> analyticSolution_;
};
}
#endif
