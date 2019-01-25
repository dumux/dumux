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
 *
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 */

#ifndef DUMUX_2PPROBLEM_HH
#define DUMUX_2PPROBLEM_HH

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

#include <dumux/material/fluidsystems/h2oair.hh>
#include <test/geomechanics/el2p/co2tables.hh>

#include "decoupled2pspatialparamsLu.hh"
//  #include "decoupled2pspatialparamsLu.hh"

#include <dumux/porousmediumflow/implicit/cpdarcyfluxvariables.hh>

namespace Dumux
{

template <class TypeTag>
class TwoP_TestProblem;

// initial conditions for mass balance equations
template<class GridView, class Scalar>
class InitialPressSat;

//////////
// Specify the properties for the lens problem
//////////s
namespace Properties
{

#if PROBLEM_IS_CC==1
NEW_TYPE_TAG(TwoP_TestProblem, INHERITS_FROM(CCTwoP, TwoPSpatialParams));
#else
NEW_TYPE_TAG(TwoP_TestProblem, INHERITS_FROM(BoxTwoP, TwoPSpatialParams));
#endif
NEW_PROP_TAG(InitialPressSat);

NEW_PROP_TAG(BaseProblem);
SET_TYPE_PROP(TwoP_TestProblem, BaseProblem, ImplicitPorousMediaProblem<TypeTag>);

SET_TYPE_PROP(TwoP_TestProblem, Grid, Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>); //Sh: elrichards: ugGrid

// Set the problem property
SET_TYPE_PROP(TwoP_TestProblem, Problem, TwoP_TestProblem<TypeTag>);

// Set fluid configuration
SET_PROP(TwoP_TestProblem, FluidSystem)
{
    typedef Dumux::H2OAirFluidSystem<TypeTag> type;
};


// Set the initial pressure and saturation function
SET_PROP(TwoP_TestProblem, InitialPressSat)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
public:
    typedef Dumux::InitialPressSat<GridView, Scalar> type;
};

// Linear solver settings
SET_INT_PROP(TwoP_TestProblem, LinearSolverVerbosity, 0);

SET_SCALAR_PROP(TwoP_TestProblem, NewtonMaxRelativeShift, 1e-5);
// SET_BOOL_PROP(TwoP_TestProblem, NewtonWriteConvergence, true);
SET_BOOL_PROP(TwoP_TestProblem, NewtonUseLineSearch, true);
// SET_BOOL_PROP(TwoP_TestProblem, NewtonUseDampedUpdate, true);
SET_INT_PROP(TwoP_TestProblem, NewtonMaxSteps, 30);
// SET_SCALAR_PROP(El2P_TestProblem, NewtonUseLineSearch, true);


SET_BOOL_PROP(TwoP_TestProblem, EvalGradientsAtSCVCenter, true);

// disable jacobian matrix recycling
SET_BOOL_PROP(TwoP_TestProblem, ImplicitEnableJacobianRecycling, false);
// disable partial reassembling
SET_BOOL_PROP(TwoP_TestProblem, ImplicitEnablePartialReassemble, false);
// Enable gravity
SET_BOOL_PROP(TwoP_TestProblem, ProblemEnableGravity, true);

// use the algebraic multigrid
SET_TYPE_PROP(TwoP_TestProblem, LinearSolver, Dumux::AMGBackend<TypeTag> );
// SET_TYPE_PROP(TwoP_TestProblem, LinearSolver, SuperLUBackend<TypeTag> );

// for better TPFA for CC use  CpDarcyFluxVariable
// Box uses DecoupledTwoPFluxVariables, which also includes the option of Keff
#if PROBLEM_IS_CC==1
SET_TYPE_PROP(TwoP_TestProblem, FluxVariables, CpDarcyFluxVariables<TypeTag>);
#else
SET_TYPE_PROP(TwoP_TestProblem, FluxVariables, DecoupledTwoPFluxVariables<TypeTag>);
#endif


// use the decoupled 2p model
SET_TYPE_PROP(TwoP_TestProblem, Model, DecoupledTwoPModel<TypeTag>);

// use the decoupled 2p volume variables
SET_TYPE_PROP(TwoP_TestProblem, VolumeVariables, DecoupledTwoPVolumeVariables<TypeTag>);

// use the decoupled 2p element volume variables
#if PROBLEM_IS_CC==1
SET_TYPE_PROP(TwoP_TestProblem, ElementVolumeVariables, DecoupledTwoPCCElementVolumeVariables<TypeTag>);
#else
SET_TYPE_PROP(TwoP_TestProblem, ElementVolumeVariables, DecoupledTwoPElementVolumeVariables<TypeTag>);
#endif

// Use the modified decoupled 2p local jacobian operator for the 2p model
SET_TYPE_PROP(TwoP_TestProblem,
              LocalResidual,
              DecoupledTwoPLocalResidual<TypeTag>);

// central differences to calculate the jacobian by default
SET_INT_PROP(TwoP_TestProblem, ImplicitNumericDifferenceMethod, 0);

// write the stress and displacement output according to rock mechanics
// sign convention (compressive stresses > 0)
// SET_BOOL_PROP(TwoP_TestProblem, VtkRockMechanicsSignConvention, true);
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
class TwoP_TestProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Indices)) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;//khodam
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;//khodam
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

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
    enum {
        // indices of the phases
            wPhaseIdx = Indices::wPhaseIdx,
            nPhaseIdx = Indices::nPhaseIdx
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


public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    TwoP_TestProblem(TimeManager &timeManager,
                const GridView &gridView)
    : ParentType(timeManager, gridView),
    gridView_(gridView),
    vertexMapper_(gridView),
    evalOriginalRhs_(false),
    evalGlobalResidual_(false)
    {
        GridCreator::grid().globalRefine(GET_RUNTIME_PARAM(TypeTag, Scalar,Grid.Refine));

        std::cout << "TwoP_TestProblem: Initializing the fluid system for the 2p model\n";

        // initialize the tables of the fluid system
         FluidSystem::init(/*Tmin=*/273,
                           /*Tmax=*/300,
                           /*nT=*/5,
                           /*pmin=*/0,
                           /*pmax=*/1e6,
                           /*np=*/10);

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

        episodeLength_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, EpisodeLengthInit);


        if (GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Output, PlotFluidMatrixInteractions))
            this->spatialParams().plotMaterialLaw();

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

     //instead of bBoxMax which was constant in the whole domain, we define the max elevation or the ground surface at each point
    Scalar zMax(const GlobalPosition &globalPos) const
    {
      if (globalPos[0]<6.3+eps_)
          return 15.;
      else if (globalPos[0]>6.3-eps_ && globalPos[0]<23.6+eps_)
          return (5.-15.)/(23.6-6.3)*(globalPos[0])+ (322./17.3);
      else if (globalPos[0]>23.6-eps_)
          return 5.;
    }

        Scalar zMin(const GlobalPosition &globalPos) const
    {
          return -25.0;
    }

    Scalar depth(const GlobalPosition &globalPos) const
    {
        return  zMax(globalPos) - globalPos[1];//distance from the surface
    }

//     Scalar waterTable = -15;

    Scalar wTdepth(const GlobalPosition &globalPos) const
    {
       return (0. - globalPos[1]);
    }


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
        //T = 308.15 /*+ (depthBOR_ - globalPos[dim-1]) * 0.025*/;
        T=283.15;

        return T;
    };


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

//         if(onInlet_(globalPos))
//           {
//             if (initializationRun_ == true)
//             {
//                 values.setDirichlet(pressureIdx, contiWEqIdx);
//                 values.setDirichlet(saturationIdx, contiNEqIdx);
//             }
//
//            }

        if ( rightBoundaryU_(globalPos)){
            values.setDirichlet(pressureIdx, contiWEqIdx);
            values.setDirichlet(saturationIdx, contiNEqIdx);//this part is aplied all the time
           }
     }

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
        values = 0.0;
//
//         if ( (onInlet_(globalPos)) ){
//           const auto& materialLawParams = this->spatialParams().materialLawParams(globalPos);
//           const Scalar swr = materialLawParams.swr();
//           const Scalar snr = materialLawParams.snr();
//
//
//           const Scalar meterUeberGW = globalPos[1] +0;
//           const Scalar pc = std::max(0.0, 9.81*1000.0*meterUeberGW);
//           const Scalar sw = std::min(1.0-snr, std::max(swr, invertPcGW_(pc, materialLawParams)));
// //           std::cout<< "invertPcGW_=" <<invertPcGW_(pc, materialLawParams) << std::endl;
//           values[pressureIdx] = (1.0e5 + 1000. * 9.81 * wTdepth(globalPos));
//           values[saturationIdx] = 1.-sw;
//         } else
//         {
          values[pressureIdx] = (1.0e5 + 1000. * 9.81 * wTdepth(globalPos));
          values[saturationIdx] = 0.0;
//         }
    }

     void solDependentNeumann(PrimaryVariables &values,
                      const Element &element,
                      const FVElementGeometry &fvGeometry,
                      const Intersection &intersection,
                      const int scvIdx,
                      const int boundaryFaceIdx,
                      const ElementVolumeVariables &elemVolVars) const
    {
        values = 0.0;
        GlobalPosition globalPos = fvGeometry.boundaryFace[boundaryFaceIdx].ipGlobal;

        Scalar avgRain_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.avgRain);//[m/s]
        Scalar rb_ = 2.0E-4; //[1/s] =K/B conductance
        Scalar ks_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.ks);//1.39e-6;

        const Scalar pW = elemVolVars[scvIdx].pressure(pressureIdx);
        const Scalar pN = elemVolVars[scvIdx].pressure(nPhaseIdx);
        const Scalar satW = elemVolVars[scvIdx].saturation(wPhaseIdx);
        const Scalar satN = elemVolVars[scvIdx].saturation(nPhaseIdx);

        if (onInlet_(globalPos))
        {
               if (satW > 1. -eps_){
                   //std::cout << "saturation above 1 \n" ;
                    values[contiWEqIdx] = rb_ * pW/(9.81); // [kg/(m2*s)]
                    if (pN>1.e5) values[contiNEqIdx] = satN * (pN - 1.e5) * rb_;
               }
                else {//deactivate to match comsol
                   //std::cout << "saturation below 1, infiltration \n" ;
                    values[contiWEqIdx] = -avgRain_*1000.;
                    if (pN>1.e5) {
                    values[contiNEqIdx] = satN * (pN - 1.e5) * rb_;
                    }
//                 }
        }

     }
    }

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

        const auto& materialLawParams = this->spatialParams().materialLawParams(globalPos);
        const Scalar swr = materialLawParams.swr();
        const Scalar snr = materialLawParams.snr();

        if (globalPos[1] > 0.)
        {
            Scalar n_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.n);
            Scalar m_ = 1.0 - (1.0 / n_);
            Scalar alpha_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.alpha);
          const Scalar meterUeberGW_ = globalPos[1] + 0;
          const Scalar pc = std::max(0.0, 9.81*1000.0*meterUeberGW_);
//        const Scalar sw = std::min(1.0-snr, std::max(swr, invertPcGW_(pc, materialLawParams)));// use the numerical solution for Pc
           const Scalar sw = std::min(1.0 - snr, std::max(swr, pow(pow(alpha_*pc, n_) + 1, -m_)));//khodam use the analytical solution for sw

          values[pressureIdx] = (1.0e5 + 1000. * 9.81 * wTdepth(globalPos));
          values[saturationIdx] = 1 - sw;
        } else
        {
          values[pressureIdx] = (1.0e5 + 1000. * 9.81 * wTdepth(globalPos));
          values[saturationIdx] = 0.0;
        }
    }
//numerical solution for Pc
//     static Scalar invertPcGW_(const Scalar pcIn,
//                               const MaterialLawParams &pcParams)
//     {
//         Scalar lower(0.0);
//         Scalar upper(1.0);
//         const unsigned int maxIterations = 25;
//         const Scalar bisLimit = 1.0;
//
//         Scalar sw, pcGW;
//         for (unsigned int k = 1; k <= maxIterations; k++)
//         {
//             sw = 0.5*(upper + lower);
//             pcGW = MaterialLaw::pc(pcParams, sw);
//             const Scalar delta = std::abs(pcGW - pcIn);
//             if (delta < bisLimit)
//                 return sw;
//
//             if (k == maxIterations)
//                 return sw;
//
//             if (pcGW > pcIn)
//                 lower = sw;
//             else
//                 upper = sw;
//         }
//         return sw;
//     }


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

    }
        void preTimeStep()
    {
        this->spatialParams().setEpisode(this->timeManager().episodeIndex());
        this->timeManager().startNextEpisode(episodeLength_);
        std::cout << "2p: episodeLength_ is " << episodeLength_ << "\n";
        this->timeManager().setTimeStepSize(episodeLength_);
        std::cout << "2p: TimeStepSize_ " << this->timeManager().timeStepSize() << "\n";
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
    static constexpr Scalar eps_ = 3e-6;
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

    bool onInlet_(const GlobalPosition &globalPos) const
    {
        return (globalPos[1]>zMax(globalPos)-eps_);
    }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return  (globalPos[1] < zMin(globalPos) + eps_);
    }

    bool leftBoundaryO_(const GlobalPosition &globalPos) const
    {
        return (globalPos[0] < -32. + eps_ && globalPos[1] > 0. - eps_);
    }

    bool leftBoundaryU_(const GlobalPosition &globalPos) const
    {
        return (globalPos[0] < -32. + eps_ && globalPos[1] < 0. + eps_);
    }

    bool rightBoundaryO_(const GlobalPosition &globalPos) const
    {
        return (globalPos[0] > 35. - eps_ && globalPos[1] > 0. - eps_);
    }

    bool rightBoundaryU_(const GlobalPosition &globalPos) const
    {
        return (globalPos[0] > 35. - eps_ && globalPos[1] < 0. + eps_);
    }
public:
    bool initializationRun_, coupled_, output_;
    InitialStressField initialStressField_;
};

} //end namespace

#endif
