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
/**
 * \file
 * \brief Definition of a problem, for the 1p2c problem:
 * Component transport of nitrogen dissolved in the water phase.
 */
#ifndef DUMUX_1P2C_CONDUCTION_PROBLEM_HH
#define DUMUX_1P2C_CONDUCTION_PROBLEM_HH

//#include <dune/grid/alugrid/2d/alugrid.hh>
#if HAVE_UG
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#endif
#include <math.h>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/implicit/1p2c/1p2cmodel.hh>
#include <dumux/implicit/common/implicitporousmediaproblem.hh>

#include <dumux/material/fluidsystems/h2on2liquidphasefluidsystem.hh>
#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>
#include "1p2coutflowspatialparams.hh"

#define NONISOTHERMAL 1

namespace Dumux
{

template <class TypeTag>
class OnePTwoCConductionProblem;

namespace Properties
{
#if NONISOTHERMAL
NEW_TYPE_TAG(OnePTwoCConductionProblem, INHERITS_FROM(OnePTwoCNI));
NEW_TYPE_TAG(OnePTwoCConductionBoxProblem, INHERITS_FROM(BoxModel, OnePTwoCConductionProblem));
NEW_TYPE_TAG(OnePTwoCConductionCCProblem, INHERITS_FROM(CCModel, OnePTwoCConductionProblem));
#else
NEW_TYPE_TAG(OnePTwoCConductionProblem, INHERITS_FROM(OnePTwoC));
NEW_TYPE_TAG(OnePTwoCConductionBoxProblem, INHERITS_FROM(BoxModel, OnePTwoCConductionProblem));
NEW_TYPE_TAG(OnePTwoCConductionCCProblem, INHERITS_FROM(CCModel, OnePTwoCConductionProblem));
#endif

// Set the grid type
SET_PROP(OnePTwoCConductionProblem, Grid)
{
#if HAVE_UG
    typedef Dune::UGGrid<2> type;
#else
//    typedef Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming> type;
    typedef Dune::YaspGrid<2> type;
#endif
};

// Set the problem property
SET_PROP(OnePTwoCConductionProblem, Problem)
{
    typedef Dumux::OnePTwoCConductionProblem<TypeTag> type;
};

// Set fluid configuration
SET_PROP(OnePTwoCConductionProblem, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::FluidSystems::H2ON2LiquidPhase<Scalar, false> type;
};

// Set the spatial parameters
SET_TYPE_PROP(OnePTwoCConductionProblem,
              SpatialParams,
              Dumux::OnePTwoCOutflowSpatialParams<TypeTag>);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(OnePTwoCConductionProblem, UseMoles, true);

// Enable velocity output
SET_BOOL_PROP(OnePTwoCConductionProblem, VtkAddVelocity, true);

// Disable gravity
SET_BOOL_PROP(OnePTwoCConductionProblem, ProblemEnableGravity, false);
}


/*!
 * \ingroup OnePTwoCModel
 * \ingroup ImplicitTestProblems
 *
 * \brief Definition of a problem, for the 1p2c problem:
 * Nitrogen is dissolved in the water phase and
 * is transported with the water flow from the left side to the right.
 *
 * The model domain is 1 m times 1 m with a discretization length of 0.05 m
 * and homogeneous soil properties (\f$ \mathrm{K=10e-10, \Phi=0.4, \tau=0.28}\f$).
 * Initially the domain is filled with pure water.
 *
 * At the left side, a Dirichlet condition defines a nitrogen mole fraction
 * of 0.3 mol/mol.
 * The water phase flows from the left side to the right due to the applied pressure
 * gradient of 1e5 Pa/m. The nitrogen is transported with the water flow
 * and leaves the domain at the right boundary
 * where an outflow boundary condition is applied.
 * 
 * The model is able to use either mole or mass fractions. The property useMoles can be set to either true or false in the
 * problem file. Make sure that the according units are used in the problem setup. The default setting for useMoles is true.
 *
 * This problem uses the \ref OnePTwoCModel model.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box1p2c -parameterFile ./test_box1p2c.input</tt> or 
 * <tt>./test_cc1p2c -parameterFile ./test_cc1p2c.input</tt>
 */
template <class TypeTag>
class OnePTwoCConductionProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel) ThermalConductivityModel;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;





    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        // world dimension
        dimWorld = GridView::dimensionworld
    };

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dimWorld : 0 };

    enum {
        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,
        massOrMoleFracIdx = Indices::massOrMoleFracIdx,
#if NONISOTHERMAL
        temperatureIdx = Indices::temperatureIdx
#endif
    };
    enum {
        // index of the transport equation
        conti0EqIdx = Indices::conti0EqIdx,
        transportEqIdx = Indices::transportEqIdx,
#if NONISOTHERMAL
        energyEqIdx = Indices::energyEqIdx
#endif
    };


    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::Intersection Intersection;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    //! property that defines whether mole or mass fractions are used
        static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

public:
    OnePTwoCConductionProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
        , eps_(1e-6)
    {
        //initialize fluid system
        FluidSystem::init();

        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, 
                                             std::string, 
                                             Problem, 
                                             Name);
        outputInterval_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                  int, Problem, OutputInterval);

        temperatureHigh_ = 300.;

        //stateing in the console whether mole or mass fractions are used
        if(!useMoles)
        {
        	std::cout<<"problem uses mass-fractions"<<std::endl;
        }
        else
        {
        	std::cout<<"problem uses mole-fractions"<<std::endl;
        }
    }


    bool shouldWriteOutput() const
    {
        return
            this->timeManager().timeStepIndex() == 0 ||
            this->timeManager().timeStepIndex() % outputInterval_ == 0 ||
            this->timeManager().episodeWillBeOver() ||
            this->timeManager().willBeFinished();
    }

    /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     */
    void addOutputVtkFields()
        {
        //Here we calculate the analytical solution
            typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;
            unsigned numDofs = this->model().numDofs();

            //create required scalar fields
            ScalarField *temperatureExact = this->resultWriter().allocateManagedBuffer(numDofs);

            FVElementGeometry fvGeometry;
            VolumeVariables volVars;

            ElementIterator eIt = this->gridView().template begin<0>();
            ElementIterator eEndIt = this->gridView().template end<0>();
            fvGeometry.update(this->gridView(), *eIt);
            PrimaryVariables initialPriVars(0);
            GlobalPosition globalPos(0);
            initial_(initialPriVars, globalPos);

            //update the constant volume variables
            volVars.update(initialPriVars,
                           *this,
                           *eIt,
                           fvGeometry,
                           0,
                           false);

            Scalar porosity = this->spatialParams().porosity(*eIt, fvGeometry, 0);
            Scalar densityW = volVars.density();
            Scalar heatCapacityW = FluidSystem::heatCapacity(volVars.fluidState(), 0);
            Scalar effectiveHeatCapacityS = this->spatialParams().heatCapacity(*eIt, fvGeometry, 0);
            Scalar storage = densityW*heatCapacityW*porosity + effectiveHeatCapacityS;
            Scalar thermalConductivityS = this->spatialParams().thermalConductivitySolid(*eIt, fvGeometry, 0);
            Scalar effectiveThermalConductivity = ThermalConductivityModel::effectiveThermalConductivity(volVars, thermalConductivityS, porosity);
            Scalar time = std::max(this->timeManager().time() + this->timeManager().timeStepSize(), 1e-10);


            for (; eIt != eEndIt; ++eIt)
            {
                fvGeometry.update(this->gridView(), *eIt);
                for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
                {
                    int globalIdx = this->model().dofMapper().map(*eIt, scvIdx, dofCodim);
                    if (isBox)
                        globalPos = eIt->geometry().corner(scvIdx);
                    else
                        globalPos = eIt->geometry().center();

                    (*temperatureExact)[globalIdx] = temperatureHigh_ + (initialPriVars[temperatureIdx] - temperatureHigh_)
                                                     *std::erf(0.5*std::sqrt(globalPos[0]*globalPos[0]*storage/time/effectiveThermalConductivity));
                }
            }
            this->resultWriter().attachDofData(*temperatureExact, "temperatureExact", isBox);

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
        return name_.c_str();
    }

#if !NONISOTHERMAL
    /*!
     * \brief Returns the temperature within the domain [K].
     *
     * This problem assumes a temperature of 20 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 20; }; // in [K]

#endif

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param values The boundary types for the conservation equations
     * \param globalPos The position for which the bc type should be evaluated
     */
    void boundaryTypesAtPos(BoundaryTypes &values, 
                            const GlobalPosition &globalPos) const
    {
        if(globalPos[0] < eps_ || globalPos[0] > this->bBoxMax()[0] - eps_)
        {
            values.setAllDirichlet();
        }
        else
        {
            values.setAllNeumann();
        }
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The position for which the bc type should be evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);

//        condition for the N2 molefraction at left boundary
        if (globalPos[0] < eps_)
            {
#if NONISOTHERMAL
                values[temperatureIdx] = temperatureHigh_;
#endif
            }


    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann
     *        boundary segment.
     *
     * For this method, the \a priVars parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     *
     * The units must be according to either using mole or mass fractions. (mole/(m^2*s) or kg/(m^2*s))
     */
    void neumann(PrimaryVariables &priVars,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const Intersection &intersection,
                 const int scvIdx,
                 const int boundaryFaceIdx) const
    {
        priVars = 0;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * For this method, the \a priVars parameter stores the rate mass
     * of a component is generated or annihilate per volume
     * unit. Positive values mean that mass is created, negative ones
     * mean that it vanishes.
     *
     * The units must be according to either using mole or mass fractions. (mole/(m^3*s) or kg/(m^3*s))
     */
    void sourceAtPos(PrimaryVariables &priVars,
                     const GlobalPosition &globalPos) const
    {
        priVars = Scalar(0.0);
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The initial values for the primary variables
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    void initialAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
    }

    // \}

private:
    // the internal method for the initial condition
    void initial_(PrimaryVariables &priVars,
                  const GlobalPosition &globalPos) const
    {
        priVars[pressureIdx] = 1e5; // initial condition for the pressure
        priVars[massOrMoleFracIdx] = 1e-5;  // initial condition for the N2 molefraction
#if NONISOTHERMAL
        priVars[temperatureIdx] = 290.;
#endif
    }

    Scalar temperatureHigh_;
    const Scalar eps_;
    std::string name_;
    int outputInterval_;
};

} //end namespace
#endif
