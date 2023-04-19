//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MPNCTests
 * \brief Problem where hot, pure liquid water is injected from the left hand side into an initially
 *        isothermal domain.
 *
 * The water is fully evaporated by a strong heat source.
 * A local thermal non-equilibrium model is used: i.e. two different (fluid, solid)
 * temperatures are primary variables.
 *
 * \author Philipp Nuske
 */

#ifndef DUMUX_COMBUSTION_PROBLEM_ONE_COMPONENT_HH
#define DUMUX_COMBUSTION_PROBLEM_ONE_COMPONENT_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/mpnc/initialconditionhelper.hh>

namespace Dumux {

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
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;

    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;
    using ParameterCache = typename FluidSystem::ParameterCache;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr int numComponents = ModelTraits::numFluidComponents();
    static constexpr int s0Idx = Indices::s0Idx;
    static constexpr int p0Idx = Indices::p0Idx;
    static constexpr int conti00EqIdx = Indices::conti0EqIdx;
    static constexpr int energyEq0Idx = Indices::energyEqIdx;
    static constexpr int numEnergyEqFluid = ModelTraits::numEnergyEqFluid();
    static constexpr int numEnergyEqSolid = ModelTraits::numEnergyEqSolid();
    static constexpr int energyEqSolidIdx = energyEq0Idx + numEnergyEqFluid + numEnergyEqSolid - 1;
    static constexpr int wPhaseIdx = FluidSystem::wPhaseIdx;
    static constexpr int nPhaseIdx = FluidSystem::nPhaseIdx;
    static constexpr int wCompIdx = FluidSystem::H2OIdx;
    static constexpr int nCompIdx = FluidSystem::N2Idx;

    static constexpr auto numPhases = ModelTraits::numFluidPhases();

public:
    CombustionProblemOneComponent(std::shared_ptr<const GridGeometry> gridGeometry)
        : ParentType(gridGeometry)
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
     * \brief Evaluates the source term for all balance equations within a given
     *        sub-control volume.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param elemVolVars The volume variables of the element
     * \param scv The local index of the sub-control volume
     *
     * Positive values mean that mass is created, negative ones mean that it vanishes.
     */
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
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     *
     * \param globalPos The global position
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
       return  initial_(globalPos);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the box scheme
     * \param elemVolVars The volume variables of the element
     * \param elemFluxVarsCache Flux variables caches for all faces in stencil
     * \param scvf The sub-control volume face
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    NumEqVector neumann(const Element &element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
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
     * \brief Evaluates the initial value for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     *
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
        const Scalar curPos = globalPos[0];
        const Scalar slope = (SwBoundary_ - SwOneComponentSys_)/this->spatialParams().lengthPM();
        Scalar S[numPhases];
        const Scalar thisSaturation = SwOneComponentSys_ + curPos * slope;

        S[wPhaseIdx] = SwBoundary_;
        if (inPM_(globalPos) )
            S[wPhaseIdx] = thisSaturation;
        S[nPhaseIdx] = 1. - S[wPhaseIdx];

        Scalar thisTemperature = TInitial_;
        if(onRightBoundary_(globalPos))
            thisTemperature = TRight_;

        FluidState fs;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            fs.setSaturation(phaseIdx, S[phaseIdx]);
            fs.setTemperature(thisTemperature);
        }

        //obtain pc according to saturation
        const int wettingPhaseIdx = this->spatialParams().template wettingPhaseAtPos<FluidSystem>(globalPos);
        const auto& fm = this->spatialParams().fluidMatrixInteractionAtPos(globalPos);
        const auto capPress = fm.capillaryPressures(fs, wettingPhaseIdx);
        Scalar p[numPhases];

        using std::abs;
        p[wPhaseIdx] = pnInitial_ - abs(capPress[wPhaseIdx]);
        p[nPhaseIdx] = p[wPhaseIdx] + abs(capPress[wPhaseIdx]);
        for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
            fs.setPressure(phaseIdx, p[phaseIdx]);

        fs.setMoleFraction(wPhaseIdx, wCompIdx, 1.0);
        fs.setMoleFraction(wPhaseIdx, nCompIdx, 0.0);
        fs.setMoleFraction(nPhaseIdx, wCompIdx, 1.0);
        fs.setMoleFraction(nPhaseIdx, nCompIdx, 0.0);

        const int refPhaseIdx = inPM_(globalPos) ? wPhaseIdx : nPhaseIdx;

        MPNCInitialConditionHelper<PrimaryVariables, FluidSystem, ModelTraits> helper;
        helper.solve(fs, MPNCInitialConditions::NotAllPhasesPresent{.refPhaseIdx = refPhaseIdx});
        auto priVars = helper.getPrimaryVariables(fs);

        // additionally set the temperature for thermal non-equilibrium for the phases
        priVars[energyEq0Idx] = thisTemperature;
        priVars[energyEqSolidIdx] = thisTemperature;
        return priVars;
    }

    /*!
     * \brief Returns whether the tested position is on the left boundary of the domain.
     */
    bool onLeftBoundary_(const GlobalPosition & globalPos) const
    {   return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_;}

    /*!
    * \brief Returns whether the tested position is on the right boundary of the domain.
     */
    bool onRightBoundary_(const GlobalPosition & globalPos) const
    {   return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_;}

    /*!
     * \brief Returns whether the tested position is in a specific region (right) in the domain
     *  \todo this needs to be more sophisticated in order to allow for meshes with nodes not directly on the boundary
     */
    bool onRightBoundaryPorousMedium_(const GlobalPosition & globalPos) const
    {
        using std::abs;
        return (abs(globalPos[0] - (this->spatialParams().lengthPM())) < eps_ );
    }

    /*!
     * \brief Returns whether the tested position is in the porous medium.
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

} // end namespace Dumux

#endif
