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

#include <dune/grid/onedgrid.hh>

#include <dumux/discretization/box/properties.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/mpnc/model.hh>
#include <dumux/porousmediumflow/mpnc/pressureformulation.hh>

#include <dumux/material/solidstates/compositionalsolidstate.hh>
#include <dumux/material/solidsystems/compositionalsolidphase.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysimplefluidlumping.hh>
#include <dumux/material/constraintsolvers/computefromreferencephase.hh>

#include "spatialparams.hh"
#include "combustionfluidsystem.hh"
#include "combustionlocalresidual.hh"

namespace Dumux {

template<class TypeTag>
class CombustionProblemOneComponent;

//! Custom model traits to deactivate diffusion for this test
template<int numP, int numC, MpNcPressureFormulation formulation, bool useM>
struct CombustionModelTraits : public MPNCModelTraits<numP, numC, formulation, useM>
{
    static constexpr bool enableMolecularDiffusion() { return false; }
};

namespace Properties {
// Create new type tags
namespace TTag {
struct CombustionOneComponent { using InheritsFrom = std::tuple<MPNCNonequil>; };
struct CombustionOneComponentBox { using InheritsFrom = std::tuple<CombustionOneComponent, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::CombustionOneComponent> { using type = Dune::OneDGrid; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::CombustionOneComponent>
{ using type = CombustionProblemOneComponent<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::CombustionOneComponent>
{
    using FVGridGeometry = GetPropType<TypeTag, FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = CombustionSpatialParams<FVGridGeometry, Scalar>;
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::CombustionOneComponent>
{ using type = FluidSystems::CombustionFluidsystem<GetPropType<TypeTag, Properties::Scalar>>; };

//! Set the default pressure formulation: either pw first or pn first
template<class TypeTag>
struct PressureFormulation<TypeTag, TTag::CombustionOneComponent>
{
public:
    static const MpNcPressureFormulation value = MpNcPressureFormulation::mostWettingFirst;
};

// Set the type used for scalar values
template<class TypeTag>
struct Scalar<TypeTag, TTag::CombustionOneComponent> { using type = double ; };
// quad / double

// We use different model traits for the equilibrium part because we want to deactivate diffusion
template<class TypeTag>
struct EquilibriumModelTraits<TypeTag, TTag::CombustionOneComponent>
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
public:
    using type = CombustionModelTraits< FluidSystem::numPhases,
                                        FluidSystem::numComponents,
                                        getPropValue<TypeTag, Properties::PressureFormulation>(),
                                        getPropValue<TypeTag, Properties::UseMoles>() >;
};

template<class TypeTag>
struct FluidState<TypeTag, TTag::CombustionOneComponent>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
public:
    using type = CompositionalFluidState<Scalar, FluidSystem>;
};
//#################
//changes from the default settings which also assume chemical non-equilibrium
//set the number of energyequations we want to use
template<class TypeTag>
struct NumEnergyEqFluid<TypeTag, TTag::CombustionOneComponent> { static constexpr int value = 1; };
template<class TypeTag>
struct NumEnergyEqSolid<TypeTag, TTag::CombustionOneComponent> { static constexpr int value = 1; };

// by default chemical non equilibrium is enabled in the nonequil model, switch that off here
template<class TypeTag>
struct EnableChemicalNonEquilibrium<TypeTag, TTag::CombustionOneComponent> { static constexpr bool value = false; };
//#################

template<class TypeTag>
struct SolidSystem<TypeTag, TTag::CombustionOneComponent>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ComponentOne = Dumux::Components::Constant<1, Scalar>;
    using ComponentTwo = Dumux::Components::Constant<2, Scalar>;
    static constexpr int numInertComponents = 2;
    using type = SolidSystems::CompositionalSolidPhase<Scalar, ComponentOne, ComponentTwo, numInertComponents>;
};

template<class TypeTag>
struct SolidState<TypeTag, TTag::CombustionOneComponent>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
public:
    using type = CompositionalSolidState<Scalar, SolidSystem>;
};

template<class TypeTag>
struct EnergyLocalResidual<TypeTag, TTag::CombustionOneComponent>
{ using type = CombustionEnergyLocalResidual<TypeTag>; };
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
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;
    using ParameterCache = typename FluidSystem::ParameterCache;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    enum {dimWorld = GridView::dimensionworld};
    enum {numPhases = ModelTraits::numPhases()};
    enum {numComponents = ModelTraits::numComponents()};
    enum {s0Idx = Indices::s0Idx};
    enum {p0Idx = Indices::p0Idx};
    enum {conti00EqIdx = Indices::conti0EqIdx};
    enum {energyEq0Idx = Indices::energyEqIdx};
    enum {numEnergyEqFluid = ModelTraits::numEnergyEqFluid()};
    enum {numEnergyEqSolid = ModelTraits::numEnergyEqSolid()};
    enum {energyEqSolidIdx = energyEq0Idx + numEnergyEqFluid + numEnergyEqSolid - 1};
    enum {wPhaseIdx = FluidSystem::wPhaseIdx};
    enum {nPhaseIdx = FluidSystem::nPhaseIdx};
    enum {wCompIdx = FluidSystem::H2OIdx};
    enum {nCompIdx = FluidSystem::N2Idx};

    // formulations
    static constexpr auto pressureFormulation = ModelTraits::pressureFormulation();
    static constexpr auto mostWettingFirst = MpNcPressureFormulation::mostWettingFirst;
    static constexpr auto leastWettingFirst = MpNcPressureFormulation::leastWettingFirst;

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
    NumEqVector source(const Element &element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume &scv) const
    {
        NumEqVector values(0.0);

        const auto& globalPos = scv.dofPosition();

        const Scalar volume = scv.volume();
        const Scalar numScv = fvGeometry.numScv(); // box: numSCV, cc:1

        if (time_ > coldTime_ )
        {
            if (onRightBoundaryPorousMedium_(globalPos))
            {
                // Testing the location of a vertex, but function is called for each associated scv. Compensate for that
                values[energyEqSolidIdx] = heatFluxFromRight_ / volume / numScv;
            }
         }
        return values;
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
    NumEqVector neumann(const Element &element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        const auto& globalPos = fvGeometry.scv(scvf.insideScvIdx()).dofPosition();
        const auto& scvIdx = scvf.insideScvIdx();
        const Scalar massFluxInjectedPhase = massFluxInjectedPhase_;

        FluidState fluidState;

        const Scalar pn = elemVolVars[scvIdx].pressure(nPhaseIdx);
        const Scalar pw = elemVolVars[scvIdx].pressure(wPhaseIdx);

        fluidState.setPressure(nPhaseIdx, pn);
        fluidState.setPressure(wPhaseIdx, pw);

        fluidState.setTemperature(TBoundary_);
        ParameterCache dummyCache;
        fluidState.setMoleFraction(wPhaseIdx, nCompIdx, 0.0);
        fluidState.setMoleFraction(wPhaseIdx, wCompIdx, 1.0);
        // compute density of injection phase
        const Scalar density = FluidSystem::density(fluidState,
                                                    dummyCache,
                                                    wPhaseIdx);
        fluidState.setDensity(wPhaseIdx, density);
        const Scalar molarDensity = FluidSystem::molarDensity(fluidState,
                                                              dummyCache,
                                                              wPhaseIdx);
        fluidState.setMolarDensity(wPhaseIdx, molarDensity);

        for(int phaseIdx=0; phaseIdx<numPhases; phaseIdx++) {
            const Scalar h = FluidSystem::enthalpy(fluidState,
                                                   dummyCache,
                                                   phaseIdx);
            fluidState.setEnthalpy(phaseIdx, h);
        }

        const Scalar molarFlux = massFluxInjectedPhase / fluidState.averageMolarMass(wPhaseIdx);

        if (onLeftBoundary_(globalPos))
        {
            values[conti00EqIdx + wCompIdx] = - molarFlux * fluidState.moleFraction(wPhaseIdx, wCompIdx);
            values[conti00EqIdx + nCompIdx] = - molarFlux * fluidState.moleFraction(wPhaseIdx, nCompIdx);
            values[energyEq0Idx] = - massFluxInjectedPhase * fluidState.enthalpy(wPhaseIdx);
        }
        return values;
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
        std::vector<Scalar> capPress(numPhases);

        //obtain pc according to saturation
        const auto &materialParams =
        this->spatialParams().materialLawParamsAtPos(globalPos);
        using MaterialLaw = typename ParentType::SpatialParams::MaterialLaw;
        using MPAdapter = MPAdapter<MaterialLaw, numPhases>;
        MPAdapter::capillaryPressures(capPress, materialParams, fluidState, wPhaseIdx);

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
        else
            DUNE_THROW(Dune::InvalidStateException, "CombustionProblemOneComponent does not support the chosen pressure formulation.");

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
                                         refPhaseIdx);

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
    { return !this->spatialParams().inOutFlow(globalPos); }

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
