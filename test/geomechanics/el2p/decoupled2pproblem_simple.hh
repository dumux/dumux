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
/*!F
 * \file
 *
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 */

#ifndef DUMUX_2PPROBLEM_SIMPLE_HH
#define DUMUX_2PPROBLEM_SIMPLE_HH

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/dnapl.hh>
#include <dumux/geomechanics/el2p/decoupled/2pmodel.hh>
#include <dumux/geomechanics/el2p/decoupled/2pfluxvariables.hh>
#include <dumux/geomechanics/el2p/decoupled/2plocalresidual.hh>
#include <dumux/geomechanics/el2p/decoupled/2pvolumevariables.hh>
#include <dumux/geomechanics/el2p/decoupled/2pelementvolumevariables.hh>
#include <dumux/geomechanics/el2p/decoupled/2pccelementvolumevariables.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/implicit/cellcentered/propertydefaults.hh>
#include <dumux/porousmediumflow/2p/implicit/gridadaptindicator.hh>
#include <dumux/implicit/adaptive/gridadaptinitializationindicator.hh>

#include <dumux/material/fluidsystems/brineco2.hh>
#include <test/geomechanics/el2p/co2tables.hh>

#include "decoupled2pspatialparams_simple.hh"
// #include "decoupled2pspatialparams.hh"

#include <dumux/porousmediumflow/implicit/cpdarcyfluxvariables.hh>

namespace Dumux
{

template <class TypeTag>
class TwoP_TestProblem_Simple;

// initial conditions for mass balance equations
template<class GridView, class Scalar>
class InitialPressSat;

//////////
// Specify the properties for the lens problem
//////////s
namespace Properties
{

#if PROBLEM_IS_CC==1
NEW_TYPE_TAG(TwoP_TestProblem_Simple, INHERITS_FROM(CCTwoP, TwoPSpatialParams_Simple));
#else
NEW_TYPE_TAG(TwoP_TestProblem_Simple, INHERITS_FROM(BoxTwoP, TwoPSpatialParams_Simple));
#endif
NEW_PROP_TAG(InitialPressSat);

NEW_PROP_TAG(BaseProblem);
SET_TYPE_PROP(TwoP_TestProblem_Simple, BaseProblem, ImplicitPorousMediaProblem<TypeTag>);

SET_TYPE_PROP(TwoP_TestProblem_Simple, Grid, Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>);

// Set the problem property
SET_TYPE_PROP(TwoP_TestProblem_Simple, Problem, TwoP_TestProblem_Simple<TypeTag>);

// Set fluid configuration
SET_PROP(TwoP_TestProblem_Simple, FluidSystem)
{
    typedef Dumux::BrineCO2FluidSystem<TypeTag> type;
};

// Set the CO2 table to be used; in this case not the the default table
// SET_TYPE_PROP(TwoP_TestProblem_Simple, CO2Table, Dumux::ViscoEl2P::CO2TablesAntonio);
SET_TYPE_PROP(TwoP_TestProblem_Simple, CO2Table, CO2Tables);
// Set the salinity mass fraction of the brine in the reservoir
SET_SCALAR_PROP(TwoP_TestProblem_Simple, ProblemSalinity, 0);

// Set the wetting phase
// SET_PROP(TwoP_TestProblem_Simple, WettingPhase)
// {
// private:
//     typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
// public:
//     typedef FluidSystems::LiquidPhase<Scalar, SimpleH2O<Scalar> > type;
// };
//
// // Set the non-wetting phase
// SET_PROP(TwoP_TestProblem_Simple, NonwettingPhase)
// {
// private:
//     typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
// public:
//     typedef FluidSystems::LiquidPhase<Scalar, DNAPL<Scalar> > type;
// };


// Set the initial pressure and saturation function
SET_PROP(TwoP_TestProblem_Simple, InitialPressSat)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
public:
    typedef Dumux::InitialPressSat<GridView, Scalar> type;
};

// Linear solver settings
SET_INT_PROP(TwoP_TestProblem_Simple, LinearSolverVerbosity, 0);

SET_SCALAR_PROP(TwoP_TestProblem_Simple, NewtonMaxRelativeShift, 1e-5);
// SET_BOOL_PROP(TwoP_TestProblem_Simple, NewtonWriteConvergence, true);
SET_BOOL_PROP(TwoP_TestProblem_Simple, NewtonUseLineSearch, true);
// SET_BOOL_PROP(TwoP_TestProblem_Simple, NewtonUseDampedUpdate, true);
SET_INT_PROP(TwoP_TestProblem_Simple, NewtonMaxSteps, 30);

SET_BOOL_PROP(TwoP_TestProblem_Simple, EvalGradientsAtSCVCenter, true);

// disable jacobian matrix recycling
SET_BOOL_PROP(TwoP_TestProblem_Simple, ImplicitEnableJacobianRecycling, false);
// disable partial reassembling
SET_BOOL_PROP(TwoP_TestProblem_Simple, ImplicitEnablePartialReassemble, false);
// Enable gravity
SET_BOOL_PROP(TwoP_TestProblem_Simple, ProblemEnableGravity, true);

// use the algebraic multigrid
// SET_TYPE_PROP(TwoP_TestProblem_Simple, LinearSolver, Dumux::AMGBackend<TypeTag> );
SET_TYPE_PROP(TwoP_TestProblem_Simple, LinearSolver, SuperLUBackend<TypeTag> );

// for better TPFA for CC use  CpDarcyFluxVariable
// Box uses DecoupledTwoPFluxVariables, which also includes the option of Keff
#if PROBLEM_IS_CC==1
SET_TYPE_PROP(TwoP_TestProblem_Simple, FluxVariables, CpDarcyFluxVariables<TypeTag>);
#else
SET_TYPE_PROP(TwoP_TestProblem_Simple, FluxVariables, DecoupledTwoPFluxVariables<TypeTag>);
#endif


// use the decoupled 2p model
SET_TYPE_PROP(TwoP_TestProblem_Simple, Model, DecoupledTwoPModel<TypeTag>);

// use the decoupled 2p volume variables
SET_TYPE_PROP(TwoP_TestProblem_Simple, VolumeVariables, DecoupledTwoPVolumeVariables<TypeTag>);

// use the decoupled 2p element volume variables
#if PROBLEM_IS_CC==1
SET_TYPE_PROP(TwoP_TestProblem_Simple, ElementVolumeVariables, DecoupledTwoPCCElementVolumeVariables<TypeTag>);
#else
SET_TYPE_PROP(TwoP_TestProblem_Simple, ElementVolumeVariables, DecoupledTwoPElementVolumeVariables<TypeTag>);
#endif

// Use the modified decoupled 2p local jacobian operator for the 2p model
SET_TYPE_PROP(TwoP_TestProblem_Simple,
              LocalResidual,
              DecoupledTwoPLocalResidual<TypeTag>);

// central differences to calculate the jacobian by default
SET_INT_PROP(TwoP_TestProblem_Simple, ImplicitNumericDifferenceMethod, 0);

// write the stress and displacement output according to rock mechanics
// sign convention (compressive stresses > 0)
// SET_BOOL_PROP(TwoP_TestProblem_Simple, VtkRockMechanicsSignConvention, true);
}

/*!
 * \ingroup TwoPModel
 * \ingroup ImplicitTestProblems
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 *
 * The domain is sized 6m times 4m and features a rectangular lens
 * with low permeablility which spans from (1 m , 2 m) to (4 m, 3 m)
 * and is surrounded by a medium with higher permability. Note that
 * this problem is discretized using only two dimensions, so from the
 * point of view of the two-phase model, the depth of the domain
 * implicitly is 1 m everywhere.
 *
 * On the top and the bottom of the domain neumann boundary conditions
 * are used, while dirichlet conditions apply on the left and right
 * boundaries.
 *
 * DNAPL is injected at the top boundary from 3m to 4m at a rate of
 * 0.04 kg/(s m^2), the remaining neumann boundaries are no-flow
 * boundaries.
 *
 * The dirichlet boundaries on the left boundary is the hydrostatic
 * pressure scaled by a factor of 1.125, while on the right side it is
 * just the hydrostatic pressure. The DNAPL saturation on both sides
 * is zero.
 *
 * This problem uses the \ref TwoPModel.
 *
 * This problem should typically be simulated until \f$t_{\text{end}}
 * \approx 20\,000\;s\f$ is reached. A good choice for the initial time step
 * size is \f$t_{\text{inital}} = 250\;s\f$.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box2p -parameterFile test_box2p.input</tt> or
 * <tt>./test_cc2p -parameterFile test_cc2p.input</tt>
 */
template <class TypeTag >
class TwoP_TestProblem_Simple : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Indices)) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    enum {
        // indices of the primary variables
            pressureIdx = Indices::pwIdx,
            saturationIdx = Indices::snIdx
    };
    enum {
        // indices of the equations+
            contiWEqIdx = Indices::contiWEqIdx,
            contiNEqIdx = Indices::contiNEqIdx
    };

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };


    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexMapper)) VertexMapper;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename Grid::ctype CoordScalar;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::BlockVector<GlobalPosition> InitialStressField;
    typedef Dune::FieldVector<Scalar, dim> DimVector;
    typedef Dune::BlockVector<Dune::FieldVector<double, dim> > VectorField;

//     typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalFEMSpace)) LocalFEMSpace;

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;//for inj_volume

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(CO2Table)) CO2Table;
    typedef Dumux::CO2<Scalar, CO2Table> CO2;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    TwoP_TestProblem_Simple(TimeManager &timeManager,
                const GridView &gridView)
    : ParentType(timeManager, gridView),
    gridView_(gridView),
    vertexMapper_(gridView),
    evalOriginalRhs_(false),
    evalGlobalResidual_(false)
    {
        GridCreator::grid().globalRefine(GET_RUNTIME_PARAM(TypeTag, Scalar,Grid.Refine));

        std::cout << "TwoP_TestProblem_Simple: Initializing the fluid system for the 2p model\n";

        // initialize the tables of the fluid system
//         FluidSystem::init(/*Tmin=*/273,
//                           /*Tmax=*/400,
//                           /*nT=*/120,
//                           /*pmin=*/1e5,
//                           /*pmax=*/1e8,
//                           /*np=*/200);

//         // resize the pressure field vector with the number of vertices
//         pInit_.resize(gridView.size(dim));
//         // fill the pressure field vector with zeros
//         std::fill( pInit_.begin(), pInit_.end(), 0.0 );

        // variable which determines if output should be written (initially set to false)
        output_ = false;
        // define if current run is initialization run
        // (initially set to true, will be set to false if initialization is over)
        initializationRun_ = true;
        // defines if feedback from geomechanics on flow is taken into account or not
        // (usually the coupling is switched off for the initialization run)
        coupled_ = false;

        // transfer the episode index to spatial parameters
        // (during intialization episode hydraulic different parameters might be applied)
        this->spatialParams().setEpisode(this->timeManager().episodeIndex());

        depthBOR_ = GET_RUNTIME_PARAM(TypeTag, Scalar, Injection.DepthBOR);

        episodeLength_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, EpisodeLengthInit);


        if (GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Output, PlotFluidMatrixInteractions))
            this->spatialParams().plotMaterialLaw();

//         eps_ = 3e-6;
//         brineDensity_ = 1000;
    }

    void init()
    {
        if (this->timeManager().time() < 1e-8)
        {
            // set the initial approximated hydrostatic pressure distribution
            // based on an averaged brine density
            // or based on a pressure polynomial
//             this->initializePressure();
            // output is written
//             this->setOutput(true);
            std::cout << "I am in the init() function" << std::endl;
        }

        ParentType::init();
    }

//     void initializePressure()
//     {
//         VertexIterator vIt = gridView_.template begin<dim>();
//         VertexIterator vEndIt = gridView_.template end<dim>();
//         for(; vIt != vEndIt; ++vIt)
//         {
//             int globalIdx = vertexMapper_.index(*vIt);
//             GlobalPosition globalPos = (*vIt).geometry().corner(0);
//
//             // initial approximate pressure distribution at start of initialization run
//             pInit_[globalIdx] = -(1.0e5 + (depthBOR_ - globalPos[1]) * brineDensity_ * 9.81);
//             std::cout.precision(15);
// //             std::cout << "pInit[" << globalIdx << "] = " << pInit_[globalIdx] << std::endl;
//         }
//     }

    // allows to change the coupled_ variable which defines if geomechanical feedback on flow is taken
    // into account
    void setCoupled(bool coupled)
    {
        coupled_ = coupled;
        std::cout << "coupled_ set to " << coupled_ << std::endl;
    }

    // returns the coupled_ variable which defines if geomechanical feedback on flow is taken
    // into account
    bool coupled() const
    {
        return coupled_;
    }

    // allows to change the output_ variable which defines if output is written
    void setOutput(bool output)
    {
        output_ = output;
    }

    // returns the initializationRun_ variable which defines if this is an initialization run
    bool initializationRun()
    {
        return initializationRun_;
    }

    // allows to set the initializationRun_ variable which defines if this is an initialization run
    void setInitializationRun(bool initializationRun)
    {
        initializationRun_ = initializationRun;
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const char *name() const
    {
        return "2pproblem";
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius at the ground surface
     * and a geothermal gradient of 0.03 K/m.
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    {
        Scalar T;
//         T = 308.15 /*+ (depthBOR_ - globalPos[dim-1]) * 0.025*/;
        T = 283.15 + (depthBOR_ - globalPos[dimWorld-1]) * 0.03;

        return T;
    };

    // returns the bottom of reservoir value (depth in m)
    const Scalar depthBOR() const
    {
        return depthBOR_;
    }

    // note: pInit is < 0 (just due to geomechanics sign convention applied here)
    // function which returns the initialized pressure at an arbitrary location within the element
    // called from finite element method (viscoel2plocaloperator.hh) and evaluated at Gauss points
//     Scalar pInit(const GlobalPosition& globalPos, const GlobalPosition& localPos, const Element& element) const
//     {
//         Scalar pValue = 0.0;
//
//         typename TwoP_TestProblem_Simple<TypeTag>::LocalFEMSpace feMap(this->gridView());
//         const typename LocalFEMSpace::Traits::FiniteElementType
//         &localFiniteElement = feMap.find(element.geometry().type());
//         typedef Dune::FieldVector<CoordScalar, 1> ShapeValue;
//         std::vector<ShapeValue> shapeVal;
//         localFiniteElement.localBasis().evaluateFunction(localPos, shapeVal);
//
// #if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
//         for (int i = 0; i < element.subEntities(dim); i++)
// #else
//         for (int i = 0; i < element.template count<dim>(); i++)
// #endif
//         {
//             int globalIdx = this->vertexMapper().subIndex(element, i, dim);
//             pValue += pInit_[globalIdx] * shapeVal[i];
//         }
//
//         return pValue;
//     }

    //function which returns brine density
    const Scalar brineDensity() const
    {
        return brineDensity_;
    }

//     // note: pInit is < 0
//     // function which returns initial pressure distribution
//     std::vector<Scalar> pInit()
//     {
//         return pInit_;
//     }

    // returns true if the current solution should be written to
    // disk (i.e. as a VTK file)
    // during initialization no output is written
    // during actual simulation output is written initially and
    // at episode/simulation end
    bool shouldWriteOutput()
    {
        return output_;
    }

    // returns true if the current solution should be written to
    // disk (i.e. as a drs file)
    bool shouldWriteRestartFile() const
    {
        return output_;
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param values Stores the value of the boundary type
     * \param globalPos The global position
     */
    void boundaryTypesAtPos(BoundaryTypes &values, const GlobalPosition& globalPos) const
    {
        values.setAllNeumann();

        // The solid displacement normal to the lateral boundaries is fixed.
        if(globalPos[0] < eps_ || globalPos[0] > this->bBoxMax()[0]-eps_)
        {
            if(initializationRun_ == false)
            {
                values.setDirichlet(pressureIdx);
                values.setDirichlet(saturationIdx);
            }
        }

        // for the initialization run the pressure and saturation
        // values are only given at the top boundary.
        if(globalPos[dimWorld-1] > this->bBoxMax()[dimWorld-1]-eps_)
        {
            values.setDirichlet(pressureIdx);
            values.setDirichlet(saturationIdx);
        }
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param vertex The vertex representing the "half volume on the boundary"
     *
     * For this method, the \a values parameter stores primary variables.
     */
//     void dirichlet(PrimaryVariables &values, const Intersection &intersection) const
//     {
//         const GlobalPosition globalPos = intersection.geometry().center();
//
//         values[pressureIdx] =  1.0e5 + (depthBOR_ - globalPos[1]) * brineDensity_ * 9.81;
//
//         std::cout << "values[pressureIdx] at" << globalPos << "] = " << values[pressureIdx] << std::endl;
//     }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * This function is called directly from dumux/geomechanics/viscoel2p/localoperator.hh
     * If it is renamed to dirichletAtPos it should be adjusted there as well.
     */
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
    {
        values[pressureIdx] =  1.013e5 + (depthBOR_ - globalPos[dimWorld-1]) * brineDensity_ * 9.81;
        values[saturationIdx] = 0.0;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param globalPos The position of the integration point of the boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    void neumannAtPos(PrimaryVariables &values,
                      const GlobalPosition &globalPos) const
    {
        values = 0.0;
    }
    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * \param values Stores the initial values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    void initialAtPos(PrimaryVariables &values,
                      const GlobalPosition &globalPos) const
    {
        typename GET_PROP_TYPE(TypeTag, FluidState) fluidState;

        // hydrostatic pressure
        values[pressureIdx] = 1.013e5 + (depthBOR_ - globalPos[dimWorld-1]) * brineDensity_ * 9.81;
        values[saturationIdx] = 0.0;
    }

    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * \param values Stores the initial values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    void source(PrimaryVariables &values,
            const Element &element,
            const FVElementGeometry &fvGeometry,
            int scvIdx) const
    {
        GlobalPosition globalPos = element.geometry().center();

        if (isBox)
            globalPos = element.geometry().corner(scvIdx);

        sourceAtPos(values, globalPos);
    }

    void sourceAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
    {
        values = 0.0;
        if(initializationRun_ == false){
            if(globalPos[0] > 499.99-eps_ && globalPos[0] < 500.01+eps_
                        && globalPos[1] > 499.99-eps_ && globalPos[1] < 500.01+eps_
                        && globalPos[dimWorld-1] > 499.99-eps_ && globalPos[dimWorld-1] < 500.01+eps_)
            {
                Scalar refinement = GET_RUNTIME_PARAM(TypeTag, int, Grid.Refinement);

                Scalar factor = std::pow(2, refinement) * std::pow(2, refinement);

                values[saturationIdx] = GET_RUNTIME_PARAM(TypeTag, Scalar, Injection.Rate)*factor; // injection
            }
        }
    }

    /*!
     * \brief Write mass balance information for both fluid phases
     */
    void postTimeStep()
    {
        PrimaryVariables mass;
        this->model().globalStorage(mass);
        double time = this->timeManager().time()+this->timeManager().timeStepSize();

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
//             std::cout.precision(17);
            std::cout << "TIME, MASS NPhase (kg), MASS WPhase (kg): \n"
            <<"mass: "
            <<time<< " , "
            <<mass[1] << " , "
            <<mass[0]
            <<"\n"
            <<"***************************************" <<std::endl;
        }
    }

    void advanceTimeLevel()
    {
//         // function in dumux/implicit/model:
//         // make the current solution the previous one.
//         uPrev_ = uCur_;
//         if (isBox)
//             prevHints_ = curHints_;
//
//         updatePrevHints();
        // is replaced by an empty function here to avoid uPrev_ during iterations
    }

    /*!
     * \brief Define length of next episode
     */
    void episodeEnd()
    {
        std::cout << "Episode control is excerted by the main problem coupling the subproblems" << std::endl;
    }

    void setEpisodeLength(Scalar episodeLength)
    {episodeLength_ = episodeLength; }

    /*!
       * \brief Get the effective porosity of an element in the transport problem.
    */
    Scalar getEffPorosity(const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx) const
    {
        int eIdx = this->model().elementMapper().index(element);
        return effPorosityVector_[eIdx][scvIdx];
    }

    /*!
     * \brief Set the effective porosity vector of the transport problem.
     */
    std::vector<std::vector<Scalar>> &setEffPorosity()
    { return effPorosityVector_; }

    /*!
       * \brief Get the old effective porosity of an element in the transport problem.
    */
    Scalar getEffPorosityOldTimestep(const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx) const
    {
        int eIdx = this->model().elementMapper().index(element);
        return effPorosityVectorOldTimestep_[eIdx][scvIdx];
    }

    /*!
     * \brief Set the old effective porosity vector of the transport problem.
     */
    std::vector<std::vector<Scalar>> &setEffPorosityOldTimestep()
    { return effPorosityVectorOldTimestep_; }

    /*!
       * \brief Get the old effective porosity of an element in the transport problem.
    */
    Scalar getEffPorosityOldIteration(const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx) const
    {
        int eIdx = this->model().elementMapper().index(element);
        return effPorosityVectorOldIteration_[eIdx];
    }

    /*!
     * \brief Set the effective Porosity of an element for the last iteration.
     */
    std::vector<std::vector<Scalar>> &setEffPorosityOldIteration()
    { return effPorosityVectorOldIteration_; }

    /*!
       * \brief Get the volumetricStrain of an element for the last iteration.
    */
    Scalar getDeltaVolumetricStrainOldIteration(const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx) const
    {
        int eIdx = this->model().elementMapper().index(element);
        return deltaVolumetricStrainOldIteration_[eIdx][scvIdx];
    }

    /*!
     * \brief Set the volumetricStrain of an element for the last iteration.
     */
    std::vector<std::vector<Scalar>> &setDeltaVolumetricStrainOldIteration()
    { return deltaVolumetricStrainOldIteration_; }

    /*!
       * \brief Get the volumetricStrain of an element for the last iteration.
    */
    Scalar getDeltaVolumetricStrainFullyCoupled(const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx) const
    {
        int eIdx = this->model().elementMapper().index(element);
        return deltaVolumetricStrainFullyCoupled_[eIdx][scvIdx];
    }

    /*!
     * \brief Set the volumetricStrain of an element for the last iteration.
     */
    std::vector<std::vector<Scalar>> &setDeltaVolumetricStrainFullyCoupled()
    { return deltaVolumetricStrainFullyCoupled_; }

/////////////////////////////////////////////////////////////////////////////////
    /*!
       * \brief Get pW of a scv for the previous iteration.
    */
    Scalar getpWInit_(const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx) const
    {
        int eIdx = this->model().elementMapper().index(element);
        return pWInit_[eIdx][scvIdx];
    }

    /*!
     * \brief Set pW of a scv for the previous iteration.
     */
    std::vector<std::vector<Scalar>> &setpWInit_()
    { return pWInit_; }

    /*!
       * \brief Get pN of a scv for the previous iteration.
    */
    Scalar getpNInit_(const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx) const
    {
        int eIdx = this->model().elementMapper().index(element);
        return pNInit_[eIdx][scvIdx];
    }

    /*!
     * \brief Set pN of a scv for the previous iteration.
     */
    std::vector<std::vector<Scalar>> &setpNInit_()
    { return pNInit_; }

    /*!
       * \brief Get sW of a scv for the previous iteration.
    */
    Scalar getsWInit_(const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx) const
    {
        int eIdx = this->model().elementMapper().index(element);
        return sWInit_[eIdx][scvIdx];
    }

    /*!
     * \brief Set sW of a scv for the previous iteration.
     */
    std::vector<std::vector<Scalar>> &setsWInit_()
    { return sWInit_; }
    /*!
       * \brief Get sN of a scv for the previous iteration.
    */
    Scalar getsNInit_(const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx) const
    {
        int eIdx = this->model().elementMapper().index(element);
        return sNInit_[eIdx][scvIdx];
    }

    /*!
     * \brief Set sN of a scv for the previous iteration.
     */
    std::vector<std::vector<Scalar>> &setsNInit_()
    { return sNInit_; }
/////////////////////////////////////////////////////////////////////////////////
    /*!
       * \brief Get pW of a scv for the previous iteration.
    */
    Scalar getpWOldIteration_(const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx) const
    {
        int eIdx = this->model().elementMapper().index(element);
        return pWOldIteration_[eIdx][scvIdx];
    }

    /*!
     * \brief Set pW of a scv for the previous iteration.
     */
    std::vector<std::vector<Scalar>> &setpWOldIteration_()
    { return pWOldIteration_; }

    /*!
       * \brief Get pN of a scv for the previous iteration.
    */
    Scalar getpNOldIteration_(const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx) const
    {
        int eIdx = this->model().elementMapper().index(element);
        return pNOldIteration_[eIdx][scvIdx];
    }

    /*!
     * \brief Set pN of a scv for the previous iteration.
     */
    std::vector<std::vector<Scalar>> &setpNOldIteration_()
    { return pNOldIteration_; }

    /*!
       * \brief Get sW of a scv for the previous iteration.
    */
    Scalar getsWOldIteration_(const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx) const
    {
        int eIdx = this->model().elementMapper().index(element);
        return sWOldIteration_[eIdx][scvIdx];
    }

    /*!
     * \brief Set sW of a scv for the previous iteration.
     */
    std::vector<std::vector<Scalar>> &setsWOldIteration_()
    { return sWOldIteration_; }
    /*!
       * \brief Get sN of a scv for the previous iteration.
    */
    Scalar getsNOldIteration_(const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx) const
    {
        int eIdx = this->model().elementMapper().index(element);
        return sNOldIteration_[eIdx][scvIdx];
    }

    /*!
     * \brief Set sN of a scv for the previous iteration.
     */
    std::vector<std::vector<Scalar>> &setsNOldIteration_()
    { return sNOldIteration_; }

    /*!
       * \brief Get the effective Porosity of an element for the last iteration.
    */
    Scalar getEffPermeability(const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx) const
    {
        int eIdx = this->model().elementMapper().index(element);
        return effPermeabilityVector_[eIdx][scvIdx];
    }

    /*!
     * \brief Set the permeablility vector of the transport problem.
     */
    std::vector<std::vector<Scalar>> &setEffPermeability()
    { return effPermeabilityVector_; }

   /*!
    * \brief Get the iteration number
    */
    int getIteration() const
    {
        return iteration_;
    }

    /*!
     * \brief Set the iteration number
     */
    int &setIteration()
    { return iteration_; }

    void setEvalOriginalRhs(bool evalOriginalRhs)
    {
        evalOriginalRhs_ = evalOriginalRhs;
    }

    bool evalOriginalRhs() const
    {
        return evalOriginalRhs_;
    }

    void setEvalGlobalResidual(bool evalGlobalResidual)
    {
        evalGlobalResidual_ = evalGlobalResidual;
    }

    bool evalGlobalResidual() const
    {
        return evalGlobalResidual_;
    }

private:
    static constexpr Scalar eps_ = 1e-10;
    Scalar depthBOR_;
    static constexpr Scalar brineDensity_ = 1059;
//     Scalar brineDensity_;
    Scalar episodeLength_;

//     std::vector<Scalar> pInit_;
    std::vector<std::vector<Scalar>> pInitPerScv_;
    GridView gridView_;
    VertexMapper vertexMapper_;
    std::vector<std::vector<Scalar>> effPorosityVector_, effPorosityVectorOldTimestep_, effPorosityVectorOldIteration_, effPermeabilityVector_, deltaVolumetricStrainOldIteration_, deltaVolumetricStrainFullyCoupled_;
    std::vector<std::vector<Scalar>> pWOldIteration_, pNOldIteration_, sWOldIteration_,sNOldIteration_;

    std::vector<std::vector<Scalar>> pWInit_, pNInit_, sWInit_,sNInit_;

    std::vector<std::vector<DimVector>> dUVector_;
        std::vector<Scalar> BNodeWiseAveraged_;
    int iteration_;
    bool evalOriginalRhs_, evalGlobalResidual_;
public:
    bool initializationRun_, coupled_, output_;
    InitialStressField initialStressField_;
};

} //end namespace

#endif
