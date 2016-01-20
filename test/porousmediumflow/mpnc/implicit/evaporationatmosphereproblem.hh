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
 * \file evaporationatmosphereproblem.hh
 * \ingroup MpNcBoxproblems
 *
 * \brief Problem showcasing the capabilities of the kinetic model.
 *
 *        The whole domain is porous medium, but the upper half has properties mimicing the ones of a free-flow domain.
 *        This way a poor man's coupling approach is accomplished: Without the complications of coupling,
 *        the main characteristics a porous and a free-flow domain are depicted.
 *
 *        The porous domain is bypassed with dry air. This way the equilibration process on top of the porous domain can be studied.
 *
 *        The Problem is written, such that the kinetic consideration for mass and energy can
 *        be switched of by merely setting kinetic, kineticenergy respectivly to false.
 *        Boundary and initial conditions are specified for all cases.
 *
 * \author Philipp Nuske
 */
#ifndef DUMUX_EVAPORATION_ATMOSPHERE_PROBLEM_HH
#define DUMUX_EVAPORATION_ATMOSPHERE_PROBLEM_HH

// set this to one for looking at a different concept of mass transfer (see mpnclocalresidualmasskinetic)
#define FUNKYMASSTRANSFER 0

// this sets that the relation using pc_max is used.
// i.e. - only parameters for awn, ans are given,
//      - the fit for ans involves the maximum value for pc, where Sw, awn are zero.
// setting it here, because it impacts volume variables and spatialparameters
#define USE_PCMAX 1

#include <dune/common/parametertreeparser.hh>

#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/implicit/mpnc/mpncmodelkinetic.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/implicit/mpnc/velomodelnewtoncontroller.hh>

#include <dumux/material/fluidsystems/h2on2fluidsystemkinetic.hh>
#include <dumux/io/gnuplotinterface.hh>
#include <dumux/io/plotoverline2d.hh>

#include <dumux/material/fluidstates/nonequilibriumfluidstate.hh>
//#include <dumux/material/fluidstates/nonequilibriumenergyfluidstate.hh>
//#include <dumux/material/fluidstates/nonequilibriummassfluidstate.hh>
//#include <dumux/material/fluidstates/compositionalfluidstate.hh>

#include "evaporationatmospherespatialparams.hh"


namespace Dumux
{

template <class TypeTag>
class EvaporationAtmosphereProblem;

namespace Properties
{
NEW_TYPE_TAG(EvaporationAtmosphereProblem,
             INHERITS_FROM(BoxMPNCKinetic, EvaporationAtmosphereSpatialParams));

// Set the grid type
SET_TYPE_PROP(EvaporationAtmosphereProblem, Grid, Dune::YaspGrid<2, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);

// Set the problem property
SET_TYPE_PROP(EvaporationAtmosphereProblem,
              Problem,
              Dumux::EvaporationAtmosphereProblem<TTAG(EvaporationAtmosphereProblem)>);

// Set fluid configuration
SET_PROP(EvaporationAtmosphereProblem, FluidSystem)
{
private: typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:  typedef Dumux::FluidSystems::H2ON2Kinetic<Scalar, /*useComplexRelations=*/false> type;
};

// Set the newton controller
SET_TYPE_PROP(EvaporationAtmosphereProblem,
              NewtonController,
              Dumux::VeloModelNewtonController<TypeTag>);


//! Set the default pressure formulation: either pw first or pn first
SET_INT_PROP(EvaporationAtmosphereProblem,
             PressureFormulation,
             MpNcPressureFormulation::leastWettingFirst);

// Set the type used for scalar values
SET_TYPE_PROP(EvaporationAtmosphereProblem, Scalar, double);


// Specify whether diffusion is enabled
SET_BOOL_PROP(EvaporationAtmosphereProblem, EnableDiffusion, true);

// do not use a chopped newton method in the beginning
SET_BOOL_PROP(EvaporationAtmosphereProblem, NewtonEnableChop, false);

// Enable the re-use of the jacobian matrix whenever possible?
SET_BOOL_PROP(EvaporationAtmosphereProblem, ImplicitEnableJacobianRecycling, true);

//#################
// Which Code to compile
// Specify whether there is any energy equation
SET_BOOL_PROP(EvaporationAtmosphereProblem, EnableEnergy, true );
// Specify whether the kinetic energy module is used
SET_INT_PROP(EvaporationAtmosphereProblem, NumEnergyEquations, 3);
// Specify whether the kinetic mass module is use
SET_BOOL_PROP(EvaporationAtmosphereProblem, EnableKinetic, true);
//#################

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This has to be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 */
SET_PROP(EvaporationAtmosphereProblem, FluidState){
    private:    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    private:    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
//    public: typedef Dumux::NonEquilibriumEnergyFluidState<TypeTag> type;
//    public: typedef Dumux::NonEquilibriumMassFluidState<TypeTag> type;
    public: typedef Dumux::NonEquilibriumFluidState<Scalar, FluidSystem> type;
//    public: typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> type;
};

SET_BOOL_PROP(EvaporationAtmosphereProblem, UseMaxwellDiffusion, false);

// Output variables
SET_BOOL_PROP(EvaporationAtmosphereProblem, VtkAddVelocities, true);
SET_BOOL_PROP(EvaporationAtmosphereProblem, VtkAddReynolds, true);
SET_BOOL_PROP(EvaporationAtmosphereProblem, VtkAddPrandtl, true);
SET_BOOL_PROP(EvaporationAtmosphereProblem, VtkAddNusselt, true);
SET_BOOL_PROP(EvaporationAtmosphereProblem, VtkAddBoundaryTypes, true);
SET_BOOL_PROP(EvaporationAtmosphereProblem, VtkAddInterfacialArea, true);
SET_BOOL_PROP(EvaporationAtmosphereProblem, VtkAddTemperatures, true);
SET_BOOL_PROP(EvaporationAtmosphereProblem, VtkAddxEquil, true);
SET_BOOL_PROP(EvaporationAtmosphereProblem, VtkAddDeltaP, true);

SET_BOOL_PROP(EvaporationAtmosphereProblem, VelocityAveragingInModel, true);
}

/*!
 * \ingroup MpNcBoxproblems
 *
 * \brief Problem that simulates the coupled heat and mass transfer processes resulting form the evaporation of liquid water from
 *        a porous medium sub-domain into a gas filled "quasi-freeflow" sub-domain.
 */
template <class TypeTag>
class EvaporationAtmosphereProblem
    : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    /*!
     * \brief The fluid state which is used by the volume variables to
     *        store the thermodynamic state. This should be chosen
     *        appropriately for the model ((non-)isothermal, equilibrium, ...).
     *        This can be done in the problem.
     */
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename FluidSystem::ParameterCache ParameterCache;
    enum { dim = GridView::dimension}; // Grid and world dimension
    enum { dimWorld = GridView::dimensionworld};
    enum { numPhases       = GET_PROP_VALUE(TypeTag, NumPhases)};
    enum { numComponents   = GET_PROP_VALUE(TypeTag, NumComponents)};
    enum { S0Idx = Indices::s0Idx};
    enum { p0Idx = Indices::p0Idx};
    enum { conti00EqIdx    = Indices::conti0EqIdx };
    enum { energyEq0Idx    = Indices::energyEqIdx};
    enum { numEnergyEqs    = Indices::numPrimaryEnergyVars};
    enum { wPhaseIdx       = FluidSystem::wPhaseIdx};
    enum { nPhaseIdx       = FluidSystem::nPhaseIdx};
    enum { sPhaseIdx       = FluidSystem::sPhaseIdx};
    enum { wCompIdx        = FluidSystem::H2OIdx};
    enum { nCompIdx        = FluidSystem::N2Idx};
    enum {  enableKinetic       = GET_PROP_VALUE(TypeTag, EnableKinetic)};
    enum { enableEnergy        = GET_PROP_VALUE(TypeTag, EnableEnergy)};
    enum { numEnergyEquations = GET_PROP_VALUE(TypeTag, NumEnergyEquations)};
    enum { velocityOutput             = GET_PROP_VALUE(TypeTag, VtkAddVelocities)};
    enum { velocityAveragingInModel   = GET_PROP_VALUE(TypeTag, VelocityAveragingInModel)};

    // formulations
    enum {
        pressureFormulation = GET_PROP_VALUE(TypeTag, PressureFormulation),
        mostWettingFirst    = MpNcPressureFormulation::mostWettingFirst,
        leastWettingFirst   = MpNcPressureFormulation::leastWettingFirst
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef Dune::FieldVector<typename GridView::Grid::ctype, dimWorld> GlobalPosition;

    typedef std::vector<Dune::FieldVector<Scalar, 1> >  ScalarBuffer;
    typedef std::array<ScalarBuffer, numPhases>         PhaseBuffer;
    typedef Dune::FieldVector<Scalar, dim>              DimVector;
    typedef Dune::BlockVector<DimVector>           VelocityField;
    typedef std::array<VelocityField, numPhases>        PhaseVelocityField;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    EvaporationAtmosphereProblem(TimeManager &timeManager,
            const GridView &gridView)
        : ParentType(timeManager, gridView), gnuplot_(false)
    {
        this->timeManager().startNextEpisode(24.* 3600.);
    }

    void episodeEnd()
    {
        // each day is one episode, result file at the end of each day: comparison
        this->timeManager().startNextEpisode(24.* 3600.);
    }

    void init()
    {
            eps_                    = 1e-6;
            percentOfEquil_         = GET_RUNTIME_PARAM(TypeTag, Scalar,BoundaryConditions.percentOfEquil);
            nTemperature_           = GET_RUNTIME_PARAM(TypeTag, int, FluidSystem.nTemperature);
            nPressure_              = GET_RUNTIME_PARAM(TypeTag, int, FluidSystem.nPressure);
            outputName_             = GET_RUNTIME_PARAM(TypeTag, std::string, Constants.outputName);
            nRestart_               = GET_RUNTIME_PARAM(TypeTag, int, Constants.nRestart);
            TInitial_               = GET_RUNTIME_PARAM(TypeTag, Scalar, InitialConditions.TInitial);
            SwPMInitial_            = GET_RUNTIME_PARAM(TypeTag, Scalar, InitialConditions.SwPMInitial);
            SwFFInitial_            = GET_RUNTIME_PARAM(TypeTag, Scalar, InitialConditions.SwFFInitial);
            pnInitial_              = GET_RUNTIME_PARAM(TypeTag, Scalar,InitialConditions.pnInitial);
            pnInjection_            = GET_RUNTIME_PARAM(TypeTag, Scalar,InitialConditions.pnInjection);
            TInject_                = GET_RUNTIME_PARAM(TypeTag, Scalar,BoundaryConditions.TInject);
            massFluxInjectedPhase_  = GET_RUNTIME_PARAM(TypeTag, Scalar,BoundaryConditions.massFluxInjectedPhase);

        this->spatialParams().setInputInitialize();

        initFluidSystem_();

        ParentType::init();
        this->model().initVelocityStuff();
    }

    void initFluidSystem_()
    {
        // initialize the tables of the fluid system
        FluidSystem::init(TInitial_ - 15.0, 453.15, nTemperature_, // T_min, T_max, n_T
                          0.75*pnInitial_, 2.25*pnInitial_, nPressure_); // p_min, p_max, n_p
    }

    /*!
     * \brief User defined output after the time integration
     *
     * Will be called diretly after the time integration.
     */
    void postTimeStep()
    {
        // write two plot Over Lines to text files
        // each output gets it's writer object in order to allow for different header files
        PlotOverLine2D<TypeTag> writerInterface;
        PlotOverLine2D<TypeTag> writerDepth;

        // writer points: on the pm--ff Interface to file
        GlobalPosition pointInterfaceOne ;
        pointInterfaceOne[0] = this->bBoxMin()[0] ;
        pointInterfaceOne[1] = (this->bBoxMax()[1]) /2 ;

        GlobalPosition pointInterfaceTwo ;
        pointInterfaceTwo[0] = this->bBoxMax()[0] ;
        pointInterfaceTwo[1] = (this->bBoxMax()[1]) /2 ;

        writerInterface.write(*this,
                              pointInterfaceOne,
                              pointInterfaceTwo,
                              "plotOverLineInterface");

        // writer points: cut through the domain
        GlobalPosition pointDepthOne ;
        pointDepthOne[0] = (this->bBoxMax()[0])/2 ;
        pointDepthOne[1] = this->bBoxMin()[1] ;

        GlobalPosition pointDepthTwo ;
        pointDepthTwo[0] = (this->bBoxMax()[0])/2 ;
        pointDepthTwo[1] = (this->bBoxMax()[1])/2 ;

        writerDepth.write(*this,
                          pointDepthOne,
                          pointDepthTwo,
                          false, /* do not append data*/
                          "plotOverLineIntoDepth");

        // Calculate storage terms of the individual phases
        for (int phaseIdx = 0; phaseIdx < numEnergyEqs; ++phaseIdx) {
            PrimaryVariables phaseStorage;
            /* How much conserved quantity is in the system (in the unit of the balance equation)
             * Summed up storage term for the whole system.
             */
            this->model().globalPhaseStorage(phaseStorage, phaseIdx);

            if (this->gridView().comm().rank() == 0) {
                std::cout
                    <<"Storage in  "
                    << FluidSystem::phaseName(phaseIdx)
                    << "Phase: ["
                    << phaseStorage
                    << "]"
                    << "\n";
            }
        }

        // Calculate total storage terms
        PrimaryVariables storage;
        this->model().globalStorage(storage);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout
                <<"Storage total: [" << storage << "]"
                << "\n";
        }

        // use gnuplot for plotting the line data
        gnuplot_.setInteraction(true);
        gnuplot_.reset();
        gnuplot_.setXlabel("xN2w [-]");
        gnuplot_.setYlabel("y [m]");
        std::ostringstream stream;
        stream << this->timeManager().time();
        gnuplot_.setOption("set label 'at time " + stream.str() + "' at screen 0.15,0.90 left");
        gnuplot_.setDatafileSeparator(' ');
        std::string fileName = outputName_ + "plotOverLineIntoDepth.dat";
        gnuplot_.addFileToPlot(fileName, fileName, " u 11:4 w l");
        gnuplot_.plot("plot");
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
    const std::string name() const
    { return outputName_ ; }

    /*!
     * \brief Returns the temperature \f$ K \f$
     */
    Scalar boxTemperature(const Element &element,
                          const FVElementGeometry &fvGeometry,
                          const unsigned int scvIdx) const
    { return TInitial_; }

    /*!
     * \brief Write a restart file?
     *        yes of course:
     * - every nRestart_'th steps
     */
    bool shouldWriteRestartFile() const
    {
        return this->timeManager().timeStepIndex() > 0 and
            (this->timeManager().timeStepIndex() % nRestart_  == 0);
    }

    /*!
     * \brief Write a output file?
     */
    bool shouldWriteOutput() const
    { return true; }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param bTypes The boundarytypes types for the conservation equations
     * \param globalPos The global position
     */
    void boundaryTypesAtPos(BoundaryTypes &bTypes,
                            const GlobalPosition &globalPos) const
    {
        // Default: Neumann
        bTypes.setAllNeumann();

        // To the right: let out what wants out
        if(onRightBoundary_(globalPos) and this->spatialParams().inFF_(globalPos) )
        {
            bTypes.setAllOutflow();
        }

        // Put a dirichlet somewhere: we need this for convergence
        if(onUpperBoundary_(globalPos) )
        {
            bTypes.setAllDirichlet();
        }

        // In the porous part the *temperature* is fixed on the boundary.
        // Mass however, is not allowed to pass (default neumann=0)
        if(( onLeftBoundary_(globalPos) and this->spatialParams().inPM_(globalPos) )
            or ( onRightBoundary_(globalPos) and this->spatialParams().inPM_(globalPos) )
            or ( onLowerBoundary_(globalPos) and this->spatialParams().inPM_(globalPos) ) )
        {
            for (int energyEqIdx=0; energyEqIdx< numEnergyEqs; ++energyEqIdx)
                bTypes.setDirichlet(energyEq0Idx+energyEqIdx);
        }
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param priVars Stores the Dirichlet values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} ] \f$
     * \param globalPos The global position
     *
     */
    void dirichletAtPos(PrimaryVariables &priVars,
                        const GlobalPosition &globalPos) const
    {
        initial_(priVars, globalPos);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann
     *        boundary segment in dependency on the current solution.
     *
     * \param priVars Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param intersection The intersection between element and boundary
     * \param scvIdx The local index of the sub-control volume
     * \param boundaryFaceIdx The index of the boundary face
     *
     * This method is used for cases, when the Neumann condition depends on the
     * solution and requires some quantities that are specific to the fully-implicit method.
     * The \a values store the mass flux of each phase normal to the boundary.
     * Negative values indicate an inflow.
     */
    void neumann(PrimaryVariables & priVars,
                 const Element & element,
                 const FVElementGeometry & fvGeometry,
                 const Intersection & intersection,
                 const unsigned int scvIdx,
                 const unsigned int boundaryFaceIdx) const
    {
        // this is the coordinate of the vertex
        // distinction via vertex works better in this case, because this is also how the
        // permeabilities & co are defined. This way there is only injection in the free flow and
        // not also in the last porous medium node.
        const GlobalPosition & globalPos =  fvGeometry.subContVol[scvIdx].global;

        priVars = 0.0;

        const Scalar massFluxInjectedPhase = massFluxInjectedPhase_ ;

        ParameterCache dummyCache;
        FluidState fluidState;

        for(int phaseIdx=0; phaseIdx<numPhases; phaseIdx++)
            fluidState.setPressure(phaseIdx, pnInitial_);

        if(numEnergyEquations){
            fluidState.setTemperature(nPhaseIdx, TInject_ );
            fluidState.setTemperature(wPhaseIdx, TInitial_ ); // this value is a good one, TInject does not work
        }
        else
            fluidState.setTemperature(TInject_ );

        // This solves the system of equations defining x=x(p,T)
        FluidSystem::calculateEquilibriumMoleFractions(fluidState, dummyCache) ;

        // Now let's make the air phase less than fully saturated with water
        fluidState.setMoleFraction(nPhaseIdx, wCompIdx, fluidState.moleFraction(nPhaseIdx, wCompIdx) * percentOfEquil_ ) ;
        fluidState.setMoleFraction(nPhaseIdx, nCompIdx, 1.-fluidState.moleFraction(nPhaseIdx, wCompIdx) ) ;

        // compute density of injection phase
        const Scalar density = FluidSystem::density(fluidState,
                                                     dummyCache,
                                                     nPhaseIdx);
        fluidState.setDensity(nPhaseIdx, density);

        if(numEnergyEquations){
            for(int phaseIdx=0; phaseIdx<numPhases; phaseIdx++){
                const Scalar h = FluidSystem::enthalpy(fluidState,
                                                       dummyCache,
                                                       phaseIdx);
                fluidState.setEnthalpy(phaseIdx, h);
            }
        }
        else{
            const Scalar h = FluidSystem::enthalpy(fluidState,
                                                   dummyCache,
                                                   nPhaseIdx);
            fluidState.setEnthalpy(nPhaseIdx, h);
        }

        const Scalar molarFlux = massFluxInjectedPhase / fluidState.averageMolarMass(nPhaseIdx);

        // actually setting the fluxes
        if(onLeftBoundary_(globalPos) and this->spatialParams().inFF_(globalPos)){
            priVars[conti00EqIdx + nPhaseIdx * numComponents + wCompIdx]
             = -molarFlux * fluidState.moleFraction(nPhaseIdx, wCompIdx);
            priVars[conti00EqIdx + nPhaseIdx * numComponents + nCompIdx]
             = -molarFlux * fluidState.moleFraction(nPhaseIdx, nCompIdx);
            // energy equations are specified mass specifically
            if(numEnergyEquations){
                priVars[energyEq0Idx + nPhaseIdx] = - massFluxInjectedPhase
                                                        * fluidState.enthalpy(nPhaseIdx) ;
            }
            else if(enableEnergy)
                priVars[energyEq0Idx] = - massFluxInjectedPhase
                                         * fluidState.enthalpy(nPhaseIdx) ;
        }
    }

    // \}


    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param priVars Stores the initial solution for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    void initialAtPos(PrimaryVariables &priVars,
                      const GlobalPosition &globalPos) const
    {
        initial_(priVars, globalPos);
    }

    /*!
     * \brief Evaluate the source term for all balance equations within a given
     *        sub-control-volume.
     *
     * \param priVars Stores the solution for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} / (m^\textrm{dim} \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     *
     *      Positive values mean that mass is created, negative ones mean that it vanishes.
     */
    void source(PrimaryVariables & priVars,
                const Element & element,
                const FVElementGeometry & fvGeometry,
                const unsigned int scvIdx) const
    {
        priVars = Scalar(0.0);
    }

    // \}


private:
    // the internal method for the initial condition
    void initial_(PrimaryVariables &priVars,
                  const GlobalPosition &globalPos) const
    {
        const Scalar T = TInitial_;
        Scalar S[numPhases];

        if (this->spatialParams().inPM_(globalPos)){
            S[wPhaseIdx]    = SwPMInitial_;
            S[nPhaseIdx]    = 1. - S[wPhaseIdx] ;
        }
        else if (this->spatialParams().inFF_(globalPos)){
            S[wPhaseIdx]    = SwFFInitial_;
            S[nPhaseIdx]    = 1. - S[wPhaseIdx] ;
        }
        else
            DUNE_THROW(Dune::InvalidStateException,
                       "You should not be here: x=" << globalPos[0] << " y= "<< globalPos[dimWorld-1]);

        for (int i = 0; i < numPhases - 1; ++i)
            priVars[S0Idx + i] = S[i];

        // capillary pressure Params
        FluidState equilibriumFluidState;

        //set saturation to inital values, this needs to be done in order for the fluidState to tell me pc
        for (int phaseIdx = 0; phaseIdx < numPhases ; ++phaseIdx) {
            equilibriumFluidState.setSaturation(phaseIdx, S[phaseIdx]);

            if(numEnergyEquations){
                equilibriumFluidState.setTemperature(phaseIdx, TInitial_ );
            }
            else
                equilibriumFluidState.setTemperature(TInject_ );
        }

        const MaterialLawParams &materialParams =
            this->spatialParams().materialLawParamsAtPos(globalPos);
        Scalar capPress[numPhases];

        //obtain pc according to saturation
        MaterialLaw::capillaryPressures(capPress, materialParams, equilibriumFluidState);

        Scalar p[numPhases];
        if (this->spatialParams().inPM_(globalPos)){
            // Use homogenous pressure in the domain and let the newton find the pressure distribution
            p[wPhaseIdx] = pnInitial_  - std::abs(capPress[wPhaseIdx]);
            p[nPhaseIdx] = p[wPhaseIdx] + std::abs(capPress[wPhaseIdx]);
        }
        else if (this->spatialParams().inFF_(globalPos)){
            p[nPhaseIdx] = pnInitial_ ;
            p[wPhaseIdx] = pnInitial_  ;
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "You should not be here: x=" << globalPos[0] << " y= "<< globalPos[dimWorld-1]);

        if(pressureFormulation == mostWettingFirst){
            // This means that the pressures are sorted from the most wetting to the least wetting-1 in the primary variables vector.
            // For two phases this means that there is one pressure as primary variable: pw
            priVars[p0Idx] = p[wPhaseIdx];
        }
        else if(pressureFormulation == leastWettingFirst){
            // This means that the pressures are sorted from the least wetting to the most wetting-1 in the primary variables vector.
            // For two phases this means that there is one pressure as primary variable: pn
            priVars[p0Idx] = p[nPhaseIdx];
        }
        else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << pressureFormulation << " is invalid.");

        // temperature
        if(enableEnergy or numEnergyEquations)
            for (int energyEqIdx=0; energyEqIdx< numEnergyEqs; ++energyEqIdx)
                priVars[energyEq0Idx + energyEqIdx] = T;

        for (int phaseIdx=0; phaseIdx<numPhases; phaseIdx++)
             equilibriumFluidState.setPressure(phaseIdx, p[phaseIdx]);

         // This solves the system of equations defining x=x(p,T)
        ParameterCache dummyCache;
        FluidSystem::calculateEquilibriumMoleFractions(equilibriumFluidState, dummyCache) ;

        FluidState dryFluidState(equilibriumFluidState);
        // Now let's make the air phase less than fully saturated with vapor
        dryFluidState.setMoleFraction(nPhaseIdx, wCompIdx, dryFluidState.moleFraction(nPhaseIdx, wCompIdx) * percentOfEquil_ ) ;
        dryFluidState.setMoleFraction(nPhaseIdx, nCompIdx, 1.-dryFluidState.moleFraction(nPhaseIdx, wCompIdx) ) ;

        /* Difference between kinetic and MPNC:
         * number of component related primVar and how they are calculated (mole fraction, fugacities, resp.)
         */
        if(enableKinetic){
            for(int phaseIdx=0; phaseIdx < numPhases; ++ phaseIdx)
            {
                for(int compIdx=0; compIdx <numComponents; ++compIdx){
                    int offset = compIdx + phaseIdx * numComponents  ;

                    if (this->spatialParams().inPM_(globalPos)){
                        priVars[conti00EqIdx + offset] = equilibriumFluidState.moleFraction(phaseIdx,compIdx) ;
                    }
                    else if (this->spatialParams().inFF_(globalPos)){
                        priVars[conti00EqIdx + offset] = dryFluidState.moleFraction(phaseIdx,compIdx) ;
                    }
                    else
                        DUNE_THROW(Dune::InvalidStateException, "You should not be here: x=" << globalPos[0] << " y= "<< globalPos[dimWorld-1]);
                }
            }
        }
        else{ //enableKinetic
            // in the case I am using the "standard" mpnc model, the variables to be set are the "fugacities"
            const Scalar fugH2O = FluidSystem::H2O::vaporPressure(T) ;
            const Scalar fugN2 = p[nPhaseIdx] - fugH2O ;

            priVars[conti00EqIdx + FluidSystem::N2Idx] = fugN2 ;
            priVars[conti00EqIdx + FluidSystem::H2OIdx] = fugH2O ;

            Scalar xl[numComponents];
            Scalar beta[numComponents];

            const Scalar Henry              = Dumux::BinaryCoeff::H2O_N2::henry(TInitial_);
            const Scalar satVapPressure     = FluidSystem::H2O::vaporPressure(TInitial_);
            xl[FluidSystem::H2OIdx]         = x_[wPhaseIdx][wCompIdx];
            xl[FluidSystem::N2Idx]          = x_[wPhaseIdx][nCompIdx];
            beta[FluidSystem::H2OIdx]       = satVapPressure ;
            beta[FluidSystem::N2Idx]        = Henry ;

            for (int i = 0; i < numComponents; ++i)
                priVars[conti00EqIdx + i] = xl[i]*beta[i]; // this should be really fug0Idx but the compiler only knows one or the other
        }
    }

    /*!
     * \brief Give back whether the tested position (input) is a specific region (left) in the domain
     */
    bool onLeftBoundary_(const GlobalPosition & globalPos) const
    {       return globalPos[0] < this->bBoxMin()[0] + eps_;   }

    /*!
     * \brief Give back whether the tested position (input) is a specific region (right) in the domain
     */
    bool onRightBoundary_(const GlobalPosition & globalPos) const
    {        return globalPos[0] > this->bBoxMax()[0] - eps_;    }

    /*!
     * \brief Give back whether the tested position (input) is a specific region (down, (gravityDir)) in the domain
     */
    bool onLowerBoundary_(const GlobalPosition & globalPos) const
    {        return globalPos[dimWorld-1] < this->bBoxMin()[dimWorld-1] + eps_;    }

    /*!
     * \brief Give back whether the tested position (input) is a specific region (up, (gravityDir)) in the domain
     */
    bool onUpperBoundary_(const GlobalPosition & globalPos) const
    {        return globalPos[dimWorld-1] > this->bBoxMax()[dimWorld-1] - eps_;    }

private:
    Scalar eps_;
    Scalar percentOfEquil_ ;
    int nTemperature_;
    int nPressure_;
    std::string outputName_;
    int nRestart_ ;
    Scalar heatIntoSolid_;
    Scalar TInitial_ ;
    Scalar SwPMInitial_ ;
    Scalar SwFFInitial_ ;
    Scalar SnInitial_;
    Scalar pnInitial_;
    Scalar pnInjection_;
    Dune::ParameterTree inputParameters_;
    Scalar x_[numPhases][numComponents] ;
    Dumux::GnuplotInterface<Scalar> gnuplot_;

    Scalar TInject_;

    Scalar massFluxInjectedPhase_ ;

    PhaseVelocityField  volumeDarcyVelocity_;
    PhaseBuffer         volumeDarcyMagVelocity_ ;
    ScalarBuffer        boxSurface_;

public:

    Dune::ParameterTree getInputParameters() const
    { return inputParameters_; }
};

} //end namespace

#endif
