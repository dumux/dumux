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
 * \ingroup MPNCTests
 * \brief Problem where hot, pure liquid water is injected from the left hand side into a initially
 *        isotherm domain. The water is fully evaporated by a strong heat source.
 *        A local thermal non-equilibrium model is used: i.e. two different (fluid, solid)
 *        temperatures are primary variables.
 *
 * \author Philipp Nuske
 */
#ifndef DUMUX_COMBUSTION_PROBLEM_ONE_COMPONENT_HH
#define DUMUX_COMBUSTION_PROBLEM_ONE_COMPONENT_HH

#include <dumux/discretization/box/properties.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/mpnc/model.hh>

#include <dumux/material/fluidsystems/purewatersimple.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysimplefluidlumping.hh>
#include <dumux/material/constraintsolvers/computefromreferencephase.hh>

#include "combustionspatialparams.hh"

namespace Dumux
{

template<class TypeTag>
class CombustionProblemOneComponent;

namespace Properties
{
NEW_TYPE_TAG(CombustionProblemOneComponent, INHERITS_FROM(MPNCNonequil, CombustionSpatialParams));
NEW_TYPE_TAG(CombustionProblemOneComponentBoxProblem, INHERITS_FROM(BoxModel, CombustionProblemOneComponent));

// Set the grid type
SET_TYPE_PROP(CombustionProblemOneComponent, Grid, Dune::OneDGrid);

// Set the problem property
SET_TYPE_PROP(CombustionProblemOneComponent,
               Problem,
               CombustionProblemOneComponent<TypeTag>);

SET_TYPE_PROP(CombustionProblemOneComponent,
              FluidSystem,
              FluidSystems::PureWaterSimpleFluidSystem<typename GET_PROP_TYPE(TypeTag, Scalar), /*useComplexRelations=*/false>);

//! Set the default pressure formulation: either pw first or pn first
SET_INT_PROP(CombustionProblemOneComponent,
        PressureFormulation,
        MpNcPressureFormulation::mostWettingFirst);

// Set the type used for scalar values
SET_TYPE_PROP(CombustionProblemOneComponent, Scalar, double );
// quad / double

// Specify whether diffusion is enabled
SET_BOOL_PROP(CombustionProblemOneComponent, EnableMolecularDiffusion, false);

//! Franz Lindners simple lumping
SET_PROP(CombustionProblemOneComponent, ThermalConductivityModel)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
public:
    using type = ThermalConductivitySimpleFluidLumping<TypeTag, Scalar, Indices>;
};

SET_PROP(CombustionProblemOneComponent, FluidState)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    using type = CompositionalFluidState<Scalar, FluidSystem>;
};
//#################
//changes from the default settings which also assume chemical non-equilibrium
//set the number of energyequations we want to use
SET_INT_PROP(CombustionProblemOneComponent, NumEnergyEqFluid, 1);
SET_INT_PROP(CombustionProblemOneComponent, NumEnergyEqSolid, 1);

// by default chemical non equilibrium is enabled in the nonequil model, switch that off here
SET_BOOL_PROP(CombustionProblemOneComponent, EnableChemicalNonEquilibrium, false);
SET_INT_PROP(CombustionProblemOneComponent, NumEqBalance, GET_PROP_VALUE(TypeTag, NumPhases)+GET_PROP_VALUE(TypeTag, NumPhases));
//#################

}
/*!
 * \ingroup MPNCTests
 * \brief Problem where water is injected from the left hand side into a porous media filled domain,
 *        which is supplied with energy from the right hand side to evaporate the water.
 */
template<class TypeTag>
class CombustionProblemOneComponent: public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using ParameterCache = typename FluidSystem::ParameterCache;
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);

    enum {dimWorld = GridView::dimensionworld};
    enum {numPhases = GET_PROP_VALUE(TypeTag, NumPhases)};
    enum {numComponents = GET_PROP_VALUE(TypeTag, NumComponents)};
    enum {s0Idx = Indices::s0Idx};
    enum {p0Idx = Indices::p0Idx};
    enum {conti00EqIdx = Indices::conti0EqIdx};
    enum {energyEq0Idx = Indices::energyEqIdx};
    enum {numEnergyEqFluid = GET_PROP_VALUE(TypeTag, NumEnergyEqFluid)};
    enum {numEnergyEqSolid = GET_PROP_VALUE(TypeTag, NumEnergyEqSolid)};
    enum {energyEqSolidIdx = energyEq0Idx + numEnergyEqFluid + numEnergyEqSolid - 1};
    enum {wPhaseIdx = FluidSystem::wPhaseIdx};
    enum {nPhaseIdx = FluidSystem::nPhaseIdx};
    enum {wCompIdx = FluidSystem::H2OIdx};
    enum {nCompIdx = FluidSystem::N2Idx};

    // formulations
    enum {
        pressureFormulation = GET_PROP_VALUE(TypeTag, PressureFormulation),
        mostWettingFirst = MpNcPressureFormulation::mostWettingFirst,
        leastWettingFirst = MpNcPressureFormulation::leastWettingFirst
    };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    CombustionProblemOneComponent(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
        : ParentType(fvGridGeometry)
    {
            outputName_ = getParam<std::string>("Problem.Name");
            nRestart_ = getParam<Scalar>("Constants.nRestart");
            TInitial_ = getParam<Scalar>("InitialConditions.TInitial");
            TRight_ = getParam<Scalar>("InitialConditions.TRight");
            pnInitial_ = getParam<Scalar>("InitialConditions.pnInitial");
            TBoundary_ = getParam<Scalar>("BoundaryConditions.TBoundary");
            SwBoundary_ = getParam<Scalar>("BoundaryConditions.SwBoundary");
            SwOneComponentSys_= getParam<Scalar>("BoundaryConditions.SwOneComponentSys");
            massFluxInjectedPhase_ = getParam<Scalar>("BoundaryConditions.massFluxInjectedPhase");
            heatFluxFromRight_ = getParam<Scalar>("BoundaryConditions.heatFluxFromRight");
            coldTime_ =getParam<Scalar>("BoundaryConditions.coldTime");
            time_ = 0.0;
    }

    void setTime(Scalar time)
    { time_ = time; }

    void setGridVariables(std::shared_ptr<GridVariables> gridVariables)
    { gridVariables_ = gridVariables; }

     const GridVariables& gridVariables() const
    { return *gridVariables_; }

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

    // \}

    /*!
     * \brief Evaluate the source term for all balance equations within a given
     *        sub-control-volume.
     *
     * \param values Stores the solution for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} / (m^\textrm{dim} \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     *
     * Positive values mean that mass is created, negative ones mean that it vanishes.
     */
    //! \copydoc Dumux::ImplicitProblem::source()
    PrimaryVariables source(const Element &element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume &scv) const
    {
        PrimaryVariables priVars(0.0);

        const auto& globalPos = scv.dofPosition();

        const Scalar volume = scv.volume();
        const Scalar numScv = fvGeometry.numScv(); // box: numSCV, cc:1

        if (time_ > coldTime_ )
        {
            if (onRightBoundaryPorousMedium_(globalPos))
            {
                // Testing the location of a vertex, but function is called for each associated scv. Compensate for that
                priVars[energyEqSolidIdx] = heatFluxFromRight_ / volume / numScv;
            }
         }
        return priVars;
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param values The boundary types for the conservation equations
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        // Default: Neumann
         bcTypes.setAllNeumann();

         if(onRightBoundary_(globalPos) ) {
            bcTypes.setAllDirichlet();
         }
        return bcTypes;
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
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
       return  initial_(globalPos);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values The neumann values for the conservation equations
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the box scheme
     * \param intersection The intersection between element and boundary
     * \param scvIdx The local vertex index
     * \param boundaryFaceIdx The index of the boundary face
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    PrimaryVariables neumann(const Element &element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolumeFace& scvf) const
    {
        PrimaryVariables priVars(0.0);

        const auto& globalPos = fvGeometry.scv(scvf.insideScvIdx()).dofPosition();
        const auto& scvIdx = scvf.insideScvIdx();
        const Scalar massFluxInjectedPhase = massFluxInjectedPhase_;

        FluidState fluidState;

        const Scalar pn = elemVolVars[scvIdx].pressure(nPhaseIdx);
        const Scalar pw = elemVolVars[scvIdx].pressure(wPhaseIdx);

        fluidState.setPressure(nPhaseIdx, pn);
        fluidState.setPressure(wPhaseIdx, pw);

        fluidState.setTemperature(TBoundary_ );
        ParameterCache dummyCache;
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

        if (onLeftBoundary_(globalPos))
        {
            priVars[conti00EqIdx + wCompIdx] = - molarFlux * fluidState.moleFraction(wPhaseIdx, wCompIdx);;
            priVars[conti00EqIdx + nCompIdx] = - molarFlux * fluidState.moleFraction(wPhaseIdx, nCompIdx);;
            priVars[energyEq0Idx] = - massFluxInjectedPhase * fluidState.enthalpy(wPhaseIdx);
        }
        return priVars;
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
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
       return initial_(globalPos);
    }

private:
    // the internal method for the initial condition
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
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
            priVars[s0Idx + i] = S[i];
        }

        FluidState fluidState;

        Scalar thisTemperature = TInitial_;
        if(onRightBoundary_(globalPos))
        thisTemperature = TRight_;

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            fluidState.setSaturation(phaseIdx, S[phaseIdx]);

            fluidState.setTemperature(thisTemperature );

        }
        //////////////////////////////////////
        // Set temperature
        //////////////////////////////////////
        priVars[energyEq0Idx] = thisTemperature;
        priVars[energyEqSolidIdx] = thisTemperature;
        Scalar capPress[numPhases];

        //obtain pc according to saturation
        const auto &materialParams =
        this->spatialParams().materialLawParamsAtPos(globalPos);
         MaterialLaw::capillaryPressures(capPress, materialParams, fluidState);

        Scalar p[numPhases];

        using std::abs;
        p[wPhaseIdx] = pnInitial_ - abs(capPress[wPhaseIdx]);
        p[nPhaseIdx] = p[wPhaseIdx] + abs(capPress[wPhaseIdx]);

        for (int phaseIdx=0; phaseIdx<numPhases; phaseIdx++)
        fluidState.setPressure(phaseIdx, p[phaseIdx]);

        //////////////////////////////////////
        // Set pressure
        //////////////////////////////////////
        if(pressureFormulation == mostWettingFirst) {
            // This means that the pressures are sorted from the most wetting to the least wetting-1 in the primary variables vector.
            // For two phases this means that there is one pressure as primary variable: pw
            priVars[p0Idx] = p[wPhaseIdx];
        }
        else if(pressureFormulation == leastWettingFirst) {
            // This means that the pressures are sorted from the least wetting to the most wetting-1 in the primary variables vector.
            // For two phases this means that there is one pressure as primary variable: pn
            priVars[p0Idx] = p[nPhaseIdx];
        }
        else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << pressureFormulation << " is invalid.");

        fluidState.setMoleFraction(wPhaseIdx, wCompIdx, 1.0);
        fluidState.setMoleFraction(wPhaseIdx, nCompIdx, 0.0);

        fluidState.setMoleFraction(nPhaseIdx, wCompIdx, 1.0);
        fluidState.setMoleFraction(nPhaseIdx, nCompIdx, 0.0);

       int refPhaseIdx;

        // on right boundary: reference is gas
        refPhaseIdx = nPhaseIdx;

        if(inPM_(globalPos)) {
            refPhaseIdx = wPhaseIdx;
        }

        // obtain fugacities
        using ComputeFromReferencePhase = ComputeFromReferencePhase<Scalar, FluidSystem>;
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
        return priVars;
    }

    /*!
     * \brief Give back whether the tested position (input) is a specific region (left) in the domain
     */
    bool onLeftBoundary_(const GlobalPosition & globalPos) const
    {   return globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + eps_;}

    /*!
     * \brief Give back whether the tested position (input) is a specific region (right) in the domain
     */
    bool onRightBoundary_(const GlobalPosition & globalPos) const
    {   return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_;}

    /*!
     * \brief Give back whether the tested position (input) is a specific region (right) in the domain
     *  \todo this needs to be more sophisticated in order to allow for meshes with nodes not directly on the boundary
     */
    bool onRightBoundaryPorousMedium_(const GlobalPosition & globalPos) const
    {
        using std::abs;
        return (abs(globalPos[0] - (this->spatialParams().lengthPM())) < eps_ );
    }

    /*!
     * \brief Give back whether the tested position (input) is a specific region (right) in the domain
     */
    bool inPM_(const GlobalPosition & globalPos) const
    {   return
        not this->spatialParams().inOutFlow(globalPos);
    }

private:
    static constexpr Scalar eps_ = 1e-6;
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
    Scalar lengthPM_;
    Scalar coldTime_;

    Scalar time_;
    std::shared_ptr<GridVariables> gridVariables_;

};

} //end namespace

#endif
