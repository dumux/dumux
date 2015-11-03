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
 * \file combustionproblem1c.hh
 *
 * \brief Problem where hot, pure liquid water is injected from the left hand side into a initially
 *        isotherm domain. The water is fully evaporated by a strong heat source.
 *        A local thermal non-equilibrium model is used: i.e. two different (fluid, solid)
 *        temperatures are primary variables.
 *
 * \ingroup MpNcBoxproblems
 *
 * \author Philipp Nuske
 */
#ifndef DUMUX_COMBUSTION_PROBLEM_ONE_COMPONENT_HH
#define DUMUX_COMBUSTION_PROBLEM_ONE_COMPONENT_HH

// st this to one for the Thermo-Gross Mass Transfer
#define FUNKYMASSTRANSFER 0

// this sets that the relation using pc_max is used.
// i.e. - only parameters for awn, ans are given,
//      - the fit for ans involves the maximum value for pc, where Sw, awn are zero.
// setting it here, because it impacts volume variables and spatialparameters
#define USE_PCMAX 0

#include <dune/common/parametertreeparser.hh>

#include <dumux/implicit/common/implicitporousmediaproblem.hh>

#include <dumux/implicit/mpnc/mpncmodelkinetic.hh>

#include "combustionspatialparams.hh"

#include <dumux/implicit/mpnc/velomodelnewtoncontroller.hh>

#include <dumux/material/fluidsystems/purewatersimplefluidsystem.hh>

#include <dumux/material/fluidstates/compositionalfluidstate.hh>

#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysimplefluidlumping.hh>

#include <dumux/material/constraintsolvers/computefromreferencephase.hh>

//#include <dumux/io/plotoverline1d.hh>

namespace Dumux {

template<class TypeTag>
class CombustionProblemOneComponent;

namespace Properties {
NEW_TYPE_TAG(CombustionProblemOneComponent,
        INHERITS_FROM(BoxMPNCKinetic, CombustionSpatialParams));

// Set the grid type
SET_TYPE_PROP(CombustionProblemOneComponent, Grid, Dune::OneDGrid);

#if HAVE_SUPERLU
SET_TYPE_PROP(CombustionProblemOneComponent, LinearSolver, SuperLUBackend<TypeTag>);
#endif

// Set the problem property
SET_TYPE_PROP(CombustionProblemOneComponent,
                Problem,
                Dumux::CombustionProblemOneComponent<TTAG(CombustionProblemOneComponent)>);

// Set fluid configuration
SET_PROP(CombustionProblemOneComponent, FluidSystem){
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::FluidSystems::PureWaterSimpleFluidSystem<Scalar, /*useComplexRelations=*/false> type;
};

// Set the newton controller
SET_TYPE_PROP(CombustionProblemOneComponent, NewtonController,
        Dumux::VeloModelNewtonController<TypeTag>);

//! Set the default pressure formulation: either pw first or pn first
SET_INT_PROP(CombustionProblemOneComponent,
        PressureFormulation,
        MpNcPressureFormulation::mostWettingFirst);

// Set the type used for scalar values
SET_TYPE_PROP(CombustionProblemOneComponent, Scalar, double );
// quad / double

// Specify whether diffusion is enabled
SET_BOOL_PROP(CombustionProblemOneComponent, EnableDiffusion, false);

// do not use a chopped newton method in the beginning
SET_BOOL_PROP(CombustionProblemOneComponent, NewtonEnableChop, false);

// Enable the re-use of the jacobian matrix whenever possible?
SET_BOOL_PROP(CombustionProblemOneComponent, ImplicitEnableJacobianRecycling, true);

// Specify whether the convergence rate ought to be written out by the
// newton method
SET_BOOL_PROP(CombustionProblemOneComponent, NewtonWriteConvergence, false);

//! Franz Lindners simple lumping
SET_TYPE_PROP(CombustionProblemOneComponent, ThermalConductivityModel, ThermalConductivitySimpleFluidLumping<TypeTag>);

//#################
// Which Code to compile
// Specify whether there is any energy equation
SET_BOOL_PROP(CombustionProblemOneComponent, EnableEnergy, true);
// Specify whether the kinetic energy module is used
SET_INT_PROP(CombustionProblemOneComponent, NumEnergyEquations, 2);
// Specify whether the kinetic mass module is use
SET_BOOL_PROP(CombustionProblemOneComponent, EnableKinetic, false);
//#################

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This has to be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 */
SET_PROP(CombustionProblemOneComponent, FluidState)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
public:
    typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> type;
};

SET_BOOL_PROP(CombustionProblemOneComponent, UseMaxwellDiffusion, false);

// Output variables
SET_BOOL_PROP(CombustionProblemOneComponent, VtkAddVelocities, true);
SET_BOOL_PROP(CombustionProblemOneComponent, VtkAddReynolds, true);
SET_BOOL_PROP(CombustionProblemOneComponent, VtkAddPrandtl, true);
SET_BOOL_PROP(CombustionProblemOneComponent, VtkAddNusselt, true);
SET_BOOL_PROP(CombustionProblemOneComponent, VtkAddBoundaryTypes, true);
SET_BOOL_PROP(CombustionProblemOneComponent, VtkAddInterfacialArea, true);
SET_BOOL_PROP(CombustionProblemOneComponent, VtkAddTemperatures, true);
SET_BOOL_PROP(CombustionProblemOneComponent, VtkAddxEquil, true);
SET_BOOL_PROP(CombustionProblemOneComponent, VtkAddDeltaP, true);

SET_BOOL_PROP(CombustionProblemOneComponent, VelocityAveragingInModel, true);
}
/*!
 * \ingroup MpNcBoxproblems
 *
 * \brief Problem where water is injected from the left hand side into a porous media filled domain,
 *        which is supplied with energy from the right hand side to evaporate the water.
 */
template<class TypeTag>
class CombustionProblemOneComponent: public ImplicitPorousMediaProblem<TypeTag> {

    typedef CombustionProblemOneComponent<TypeTag> ThisType;
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Indices)Indices;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    /*!
     * \brief The fluid state which is used by the volume variables to
     *        store the thermodynamic state. This should be chosen
     *        appropriately for the model ((non-)isothermal, equilibrium, ...).
     *        This can be done in the problem.
     */
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename FluidSystem::ParameterCache ParameterCache;
    enum {dim = GridView::dimension}; // Grid and world dimension
    enum {dimWorld = GridView::dimensionworld};
    enum {numPhases = GET_PROP_VALUE(TypeTag, NumPhases)};
    enum {numComponents = GET_PROP_VALUE(TypeTag, NumComponents)};
    enum {s0EqIdx = Indices::s0Idx};
    enum {p0EqIdx = Indices::p0Idx};
    enum {conti00EqIdx = Indices::conti0EqIdx};
    enum {temperature0Idx = Indices::temperatureIdx};
    enum {energyEq0Idx = Indices::energyEqIdx};
    enum {numEnergyEqs = Indices::numPrimaryEnergyVars};
    enum {energyEqSolidIdx = energyEq0Idx + numEnergyEqs - 1};
    enum {wPhaseIdx = FluidSystem::wPhaseIdx};
    enum {nPhaseIdx = FluidSystem::nPhaseIdx};
    enum {sPhaseIdx = FluidSystem::sPhaseIdx};
    enum {wCompIdx = FluidSystem::H2OIdx};
    enum {nCompIdx = FluidSystem::N2Idx};
    enum {enableKinetic = GET_PROP_VALUE(TypeTag, EnableKinetic)};
    enum {enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy)};
    enum {numEnergyEquations = GET_PROP_VALUE(TypeTag, NumEnergyEquations)};
    enum {velocityOutput = GET_PROP_VALUE(TypeTag, VtkAddVelocities)};
    enum {velocityAveragingInModel = GET_PROP_VALUE(TypeTag, VelocityAveragingInModel)};

    // formulations
    enum {
        pressureFormulation = GET_PROP_VALUE(TypeTag, PressureFormulation),
        mostWettingFirst = MpNcPressureFormulation::mostWettingFirst,
        leastWettingFirst = MpNcPressureFormulation::leastWettingFirst
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;

    typedef Dune::FieldVector<typename GridView::Grid::ctype, dim> DimVector;
    typedef Dune::FieldVector<typename GridView::Grid::ctype, dimWorld> GlobalPosition;

    typedef std::vector<Dune::FieldVector<Scalar, 1> > ScalarBuffer;
    typedef std::array<ScalarBuffer, numPhases> PhaseBuffer;
    typedef Dune::FieldVector<Scalar, dimWorld> VelocityVector;
    typedef Dune::BlockVector<VelocityVector> VelocityField;
    typedef std::array<VelocityField, numPhases> PhaseVelocityField;

public:
    CombustionProblemOneComponent(TimeManager &timeManager,
            const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
    }

    void init()
    {
            eps_ = 1e-6;
            outputName_ = GET_RUNTIME_PARAM(TypeTag, std::string, Constants.outputName);
            nRestart_ = GET_RUNTIME_PARAM(TypeTag, int, Constants.nRestart);
            TInitial_ = GET_RUNTIME_PARAM(TypeTag, Scalar, InitialConditions.TInitial);
            TRight_ = GET_RUNTIME_PARAM(TypeTag, Scalar, InitialConditions.TRight);
            pnInitial_ = GET_RUNTIME_PARAM(TypeTag, Scalar,InitialConditions.pnInitial);
            TBoundary_ = GET_RUNTIME_PARAM(TypeTag, Scalar,BoundaryConditions.TBoundary);
            SwBoundary_ = GET_RUNTIME_PARAM(TypeTag, Scalar,BoundaryConditions.SwBoundary);
            SwOneComponentSys_ = GET_RUNTIME_PARAM(TypeTag, Scalar,BoundaryConditions.SwOneComponentSys);
            massFluxInjectedPhase_ = GET_RUNTIME_PARAM(TypeTag, Scalar,BoundaryConditions.massFluxInjectedPhase);
            heatFluxFromRight_ = GET_RUNTIME_PARAM(TypeTag, Scalar,BoundaryConditions.heatFluxFromRight);
            coldTime_ = GET_RUNTIME_PARAM(TypeTag, Scalar,BoundaryConditions.coldTime);

        this->spatialParams().setInputInitialize();

        ParentType::init();
        this->model().initVelocityStuff();
    }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ for an isothermal problem.
     *
     * This is not specific to the discretization.
     */
    Scalar temperature() const
    {   return TInitial_;}

    /*!
     * \name Problem Params
     */
    // \{
    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string name() const
    {   return outputName_;}

    /*!
     * \brief Returns the temperature within the domain.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     *
     */
    Scalar boxTemperature(const Element &element,
            const FVElementGeometry &fvGeometry,
            const unsigned int scvIdx) const
    {   return TInitial_;}

    // \}

    /*!
     * \brief Called directly after the time integration.
     */
    void postTimeStep()
    {
//      if(this->timeManager().willBeFinished()){
//          // write plot over Line to text files
//          // each output gets it's writer object in order to allow for different header files
//          PlotOverLine1D<TypeTag> writer;
//
//          // writer points: in the porous medium, not the outflow
//          GlobalPosition pointOne ;
//          GlobalPosition pointTwo ;
//
//          pointOne[0] = this->bBoxMin()[0] ;
//          pointTwo[0] = (this->spatialParams().lengthPM()) ;
//
//          writer.write(*this,
//                        pointOne,
//                        pointTwo,
//                        "plotOverLine");
//      }

        // Calculate storage terms of the individual phases
        for (int phaseIdx = 0; phaseIdx < numEnergyEqs; ++phaseIdx) {
            PrimaryVariables phaseStorage;
            /* how much of the conserved quantity is in the system (in units of the balance eq.)
             * summed up storage term over the whole system
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
    }

    /*!
     * \brief Write a restart file?
     *yes of course:
     * - every nRestart_'th steps
     */
    bool shouldWriteRestartFile() const
    {
        return this->timeManager().timeStepIndex() > 0 and
        (this->timeManager().timeStepIndex() % nRestart_ == 0);
    }

    /*!
     * \brief Write a output file?
     */
    bool shouldWriteOutput() const
    {   return true;}

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * For this method, the \a values parameter stores the rate mass
     * of a component is generated or annihilate per volume
     * unit. Positive values mean that mass is created, negative ones
     * mean that it vanishes.
     *
     * \param priVars Stores the Dirichlet values for the conservation equations in
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
    void source(PrimaryVariables & priVars,
            const Element & element,
            const FVElementGeometry & fvGeometry,
            const unsigned int scvIdx) const
    {
        priVars = Scalar(0.0);
        const GlobalPosition & globalPos = fvGeometry.subContVol[scvIdx].global;

        const Scalar volume = fvGeometry.subContVol[scvIdx].volume;
        const unsigned int numScv = fvGeometry.numScv; // box: numSCV, cc:1

        if ((this->timeManager().time()) > coldTime_ )
        {
            if (onRightBoundaryPorousMedium_(globalPos))
            {
                if(numEnergyEquations>1) // for  2 and 3 energy equations
                {
                    if(numEnergyEquations==2) {
                        // Testing the location of a vertex, but function is called for each associated scv. Compensate for that
                        priVars[energyEqSolidIdx] = heatFluxFromRight_ / volume / (Scalar)numScv;
                    }
                    else
                    // Multiply by volume, because afterwards we divide by it -> make independet of grid resolution
                    priVars[energyEq0Idx + sPhaseIdx] = heatFluxFromRight_ / volume / (Scalar)numScv;
                }
                else if (enableEnergy) {
                    // Multiply by volume, because afterwards we divide by it -> make independet of grid resolution
                    priVars[energyEq0Idx] = heatFluxFromRight_ / volume / (Scalar)numScv;
                }
            }
        }
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param bTypes The bounentraldary types for the conservation equations
     * \param globalPos The global position

     */
    void boundaryTypesAtPos(BoundaryTypes &bTypes,
                            const GlobalPosition &globalPos) const
    {
        // Default: Neumann
        bTypes.setAllNeumann();

        if(onRightBoundary_(globalPos) ) {
            bTypes.setAllDirichlet();
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
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichletAtPos(PrimaryVariables &priVars,
                        const GlobalPosition &globalPos) const
    {
        initial_(priVars, globalPos);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     *
     * \param priVars Stores the Dirichlet values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} ] \f$
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param intersection The intersection between element and boundary
     * \param scvIdx The local index of the sub-control volume
     * \param boundaryFaceIdx The index of the boundary face
     * \param elemVolVars The element Volume Variables
     */

    void solDependentNeumann(PrimaryVariables &priVars,
            const Element &element,
            const FVElementGeometry &fvGeometry,
            const Intersection &intersection,
            const int scvIdx,
            const int boundaryFaceIdx,
            const ElementVolumeVariables &elemVolVars) const {

        const GlobalPosition & globalPos = fvGeometry.subContVol[scvIdx].global;
        const Scalar massFluxInjectedPhase = massFluxInjectedPhase_;

        FluidState fluidState;

        const Scalar pn = elemVolVars[scvIdx].fluidState().pressure(nPhaseIdx);
        const Scalar pw = elemVolVars[scvIdx].fluidState().pressure(wPhaseIdx);

        const Scalar Tn = elemVolVars[scvIdx].fluidState().temperature(nPhaseIdx);
        const Scalar Tw = elemVolVars[scvIdx].fluidState().temperature(wPhaseIdx);

        fluidState.setPressure(nPhaseIdx, pn);
        fluidState.setPressure(wPhaseIdx, pw);

        if(numEnergyEquations == 3 ) { // only in this case the fluidstate hase several temperatures
            fluidState.setTemperature(nPhaseIdx, Tn );
            fluidState.setTemperature(wPhaseIdx, Tw );
        }
        else
        fluidState.setTemperature(TBoundary_ );

        ParameterCache dummyCache;
        // This solves the system of equations defining x=x(p,T)
//           FluidSystem::calculateEquilibriumMoleFractions(fluidState, dummyCache) ;

        fluidState.setMoleFraction(wPhaseIdx, nCompIdx, 0.0);
        fluidState.setMoleFraction(wPhaseIdx, wCompIdx, 1.0);
        // compute density of injection phase
        const Scalar density = FluidSystem::density(fluidState,
                dummyCache,
                wPhaseIdx);
        fluidState.setDensity(wPhaseIdx, density);

        for(int phaseIdx=0; phaseIdx<numPhases; phaseIdx++) {
            const Scalar h = FluidSystem::enthalpy(fluidState,
                    dummyCache,
                    phaseIdx);
            fluidState.setEnthalpy(phaseIdx, h);
        }

        const Scalar molarFlux = massFluxInjectedPhase / fluidState.averageMolarMass(wPhaseIdx);

        // Default Neumann noFlow
        priVars = 0.0;

        if (onLeftBoundary_(globalPos))
        {
            if (enableKinetic) {
                priVars[conti00EqIdx + numComponents*wPhaseIdx+wCompIdx] = - molarFlux * fluidState.moleFraction(wPhaseIdx, wCompIdx);
                priVars[conti00EqIdx + numComponents*wPhaseIdx+nCompIdx] = - molarFlux * fluidState.moleFraction(wPhaseIdx, nCompIdx);
                if (enableEnergy) {
                    if (numEnergyEquations == 3)
                    priVars[energyEq0Idx + wPhaseIdx] = - massFluxInjectedPhase * fluidState.enthalpy(wPhaseIdx);
                    else if(numEnergyEquations == 2)
                    priVars[energyEq0Idx] = - massFluxInjectedPhase * fluidState.enthalpy(wPhaseIdx);
                    else if(numEnergyEquations == 1)
                    priVars[energyEq0Idx] = - massFluxInjectedPhase * fluidState.enthalpy(wPhaseIdx);
                    else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << pressureFormulation << " is invalid.");
                }
            }
            else {
                priVars[conti00EqIdx + wCompIdx] = - molarFlux * fluidState.moleFraction(wPhaseIdx, wCompIdx);;
                priVars[conti00EqIdx + nCompIdx] = - molarFlux * fluidState.moleFraction(wPhaseIdx, nCompIdx);;
                if (enableEnergy) {
                    if (numEnergyEquations == 3)
                    priVars[energyEq0Idx + wPhaseIdx] = - massFluxInjectedPhase * fluidState.enthalpy(wPhaseIdx);
                    else if(numEnergyEquations == 2)
                    priVars[energyEq0Idx] = - massFluxInjectedPhase * fluidState.enthalpy(wPhaseIdx);
                    else if(numEnergyEquations == 1)
                    priVars[energyEq0Idx] = - massFluxInjectedPhase * fluidState.enthalpy(wPhaseIdx);
                    else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << pressureFormulation << " is invalid.");
                }
            }
        }
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     *
     * \param values Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param globalPos The global position
     */
    void initialAtPos(PrimaryVariables &values,
                      const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
    }

private:
    // the internal method for the initial condition
    void initial_(PrimaryVariables &priVars,
            const GlobalPosition &globalPos) const
    {
        const Scalar curPos = globalPos[0];
        const Scalar slope = (SwBoundary_-SwOneComponentSys_) / (this->spatialParams().lengthPM());
        Scalar S[numPhases];
        const Scalar thisSaturation = SwOneComponentSys_ + curPos * slope;

        S[wPhaseIdx] = SwBoundary_;
        if (inPM_(globalPos) ) {
            S[wPhaseIdx] = thisSaturation;
        }

        S[nPhaseIdx] = 1. - S[wPhaseIdx];

        //////////////////////////////////////
        // Set saturation
        //////////////////////////////////////
        for (int i = 0; i < numPhases - 1; ++i) {
            priVars[s0EqIdx + i] = S[i];
        }

        FluidState fluidState;

        Scalar thisTemperature = TInitial_;
        if(onRightBoundary_(globalPos))
        thisTemperature = TRight_;

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            fluidState.setSaturation(phaseIdx, S[phaseIdx]);

            if(numEnergyEquations > 2) { // only in this case the fluidstate has more than one temperature
                fluidState.setTemperature(phaseIdx, thisTemperature );
            }
            else
            fluidState.setTemperature(thisTemperature );
        }

        //////////////////////////////////////
        // Set temperature
        //////////////////////////////////////
        if(enableEnergy) {
            for (int energyEqIdx=0; energyEqIdx< numEnergyEqs; ++energyEqIdx)
            priVars[energyEq0Idx + energyEqIdx] = thisTemperature;
        }

        const MaterialLawParams &materialParams =
        this->spatialParams().materialLawParamsAtPos(globalPos);
        Scalar capPress[numPhases];

        //obtain pc according to saturation
        MaterialLaw::capillaryPressures(capPress, materialParams, fluidState);

        Scalar p[numPhases];

        p[wPhaseIdx] = pnInitial_ - std::abs(capPress[wPhaseIdx]);
        p[nPhaseIdx] = p[wPhaseIdx] + std::abs(capPress[wPhaseIdx]);

        for (int phaseIdx=0; phaseIdx<numPhases; phaseIdx++)
        fluidState.setPressure(phaseIdx, p[phaseIdx]);

        //////////////////////////////////////
        // Set pressure
        //////////////////////////////////////
        if(pressureFormulation == mostWettingFirst) {
            // This means that the pressures are sorted from the most wetting to the least wetting-1 in the primary variables vector.
            // For two phases this means that there is one pressure as primary variable: pw
            priVars[p0EqIdx] = p[wPhaseIdx];
        }
        else if(pressureFormulation == leastWettingFirst) {
            // This means that the pressures are sorted from the least wetting to the most wetting-1 in the primary variables vector.
            // For two phases this means that there is one pressure as primary variable: pn
            priVars[p0EqIdx] = p[nPhaseIdx];
        }
        else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << pressureFormulation << " is invalid.");

        fluidState.setMoleFraction(wPhaseIdx, wCompIdx, 1.0);
        fluidState.setMoleFraction(wPhaseIdx, nCompIdx, 0.0);

        fluidState.setMoleFraction(nPhaseIdx, wCompIdx, 1.0);
        fluidState.setMoleFraction(nPhaseIdx, nCompIdx, 0.0);

        /* Difference between kinetic and MPNC:
         * number of component related primVar and how they are calculated (mole fraction, fugacities, resp.)          */
        if(enableKinetic) {
            for(int phaseIdx=0; phaseIdx < numPhases; ++ phaseIdx) {
                for(int compIdx=0; compIdx <numComponents; ++compIdx) {
                    int offset = compIdx + phaseIdx * numComponents;
                    priVars[conti00EqIdx + offset] = fluidState.moleFraction(phaseIdx,compIdx);
                }
            }
        }
        else { //enableKinetic
            int refPhaseIdx;

            // on right boundary: reference is gas
            refPhaseIdx = nPhaseIdx;

            if(inPM_(globalPos)) {
                refPhaseIdx = wPhaseIdx;
            }

            // obtain fugacities
            typedef Dumux::ComputeFromReferencePhase<Scalar, FluidSystem> ComputeFromReferencePhase;
            ParameterCache paramCache;
            ComputeFromReferencePhase::solve(fluidState,
                    paramCache,
                    refPhaseIdx,
                    /*setViscosity=*/false,
                    /*setEnthalpy=*/false);

            //////////////////////////////////////
            // Set fugacities
            //////////////////////////////////////
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                priVars[conti00EqIdx + compIdx] = fluidState.fugacity(refPhaseIdx,compIdx);
            }
        } // end not enableKinetic
    }

    /*!
     * \brief Give back whether the tested position (input) is a specific region (left) in the domain
     */
    bool onLeftBoundary_(const GlobalPosition & globalPos) const
    {   return globalPos[0] < this->bBoxMin()[0] + eps_;}

    /*!
     * \brief Give back whether the tested position (input) is a specific region (right) in the domain
     */
    bool onRightBoundary_(const GlobalPosition & globalPos) const
    {   return globalPos[0] > this->bBoxMax()[0] - eps_;}

    /*!
     * \brief Give back whether the tested position (input) is a specific region (right) in the domain
     *  \todo this needs to be more sophisticated in order to allow for meshes with nodes not directly on the boundary
     */
    bool onRightBoundaryPorousMedium_(const GlobalPosition & globalPos) const
    {   return ( std::fabs(globalPos[0] - (this->spatialParams().lengthPM())) < eps_ );}

    /*!
     * \brief Give back whether the tested position (input) is a specific region (right) in the domain
     */
    bool inPM_(const GlobalPosition & globalPos) const
    {   return
        not this->spatialParams().inOutFlow(globalPos);
    }

    /*!
     * \brief Give back whether the tested position (input) is a specific region (down, (gravityDir)) in the domain
     */
    bool onLowerBoundary_(const GlobalPosition & globalPos) const
    {   return globalPos[dimWorld-1] < this->bBoxMin()[dimWorld-1] + eps_;}

    /*!
     * \brief Give back whether the tested position (input) is a specific region (up, (gravityDir)) in the domain
     */
    bool onUpperBoundary_(const GlobalPosition & globalPos) const
    {   return globalPos[dimWorld-1] > this->bBoxMax()[dimWorld-1] - eps_;}

private:
    Scalar eps_;
    int nTemperature_;
    int nPressure_;
    std::string outputName_;
    int nRestart_;
    Scalar TInitial_;
    Scalar TRight_;

    Scalar pnInitial_;

    Dune::ParameterTree inputParameters_;
    Scalar x_[numPhases][numComponents];

    Scalar TBoundary_;
    Scalar SwBoundary_;
    Scalar SwOneComponentSys_;

    Scalar massFluxInjectedPhase_;
    Scalar heatFluxFromRight_;
    PhaseVelocityField volumeDarcyVelocity_;
    PhaseBuffer volumeDarcyMagVelocity_;
    ScalarBuffer boxSurface_;
    Scalar lengthPM_;
    Scalar coldTime_;

public:

    Dune::ParameterTree getInputParameters() const
    {   return inputParameters_;}
};

}
 //end namespace

#endif
