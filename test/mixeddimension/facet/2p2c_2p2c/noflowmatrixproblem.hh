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
#ifndef DUMUX_2P2C_NOFLOW_MATRIX_PROBLEM_HH
#define DUMUX_2P2C_NOFLOW_MATRIX_PROBLEM_HH

#include <dumux/mixeddimension/facet/mpfa/properties.hh>
#include <dumux/mixeddimension/subproblemproperties.hh>

#include <dumux/porousmediumflow/2p2c/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>

#include <dumux/material/fluidsystems/h2on2.hh>
#include "noflowmatrixspatialparams.hh"

namespace Dumux
{
template <class TypeTag>
class NoFlowDomainMatrixProblem;

namespace Properties
{
NEW_TYPE_TAG(NoFlowDomainMatrixProblem, INHERITS_FROM(FacetCouplingBulkMpfaModel, TwoPTwoCNI, MatrixSpatialParams));

// Set fluid configuration
SET_TYPE_PROP(NoFlowDomainMatrixProblem, FluidSystem, FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), true>);

// Set the problem property
SET_TYPE_PROP(NoFlowDomainMatrixProblem, Problem, NoFlowDomainMatrixProblem<TypeTag>);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(NoFlowDomainMatrixProblem, UseMoles, true);

// Set the spatial parameters
SET_TYPE_PROP(NoFlowDomainMatrixProblem, SpatialParams, TwoPMatrixSpatialParams<TypeTag>);

// Linear solver settings
SET_TYPE_PROP(NoFlowDomainMatrixProblem, LinearSolver, SuperLUBackend<TypeTag>);

// Enable gravity
SET_BOOL_PROP(NoFlowDomainMatrixProblem, ProblemEnableGravity, false);

// Solution-independent tensors
SET_BOOL_PROP(NoFlowDomainMatrixProblem, SolutionDependentAdvection, false);
SET_BOOL_PROP(NoFlowDomainMatrixProblem, EnableGlobalFVGeometryCache, true);
}

/*!
 * \ingroup OnePModel
 * \ingroup ImplicitTestProblems
 * \brief  Test problem for the one-phase model
 */

template <class TypeTag>
class NoFlowDomainMatrixProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    using ParentType = ImplicitPorousMediaProblem<TypeTag>;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    using GlobalProblemTypeTag = typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag);
    using CouplingManager = typename GET_PROP_TYPE(GlobalProblemTypeTag, CouplingManager);

    // copy some indices for convenience
    enum
    {
        // primary variable indices
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,

        // phase indices
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,

        // component indices
        wCompIdx = Indices::wCompIdx,
        nCompIdx = Indices::nCompIdx
    };

    //! property that defines whether mole or mass fractions are used
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

    static constexpr int numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;


public:
    NoFlowDomainMatrixProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        //initialize fluid system
        FluidSystem::init();

        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name) + "_matrix";
        eps_ = 1e-6;

        // stating in the console whether mole or mass fractions are used
        if(useMoles)
            std::cout<<"problem uses mole fractions"<<std::endl;
        else
            std::cout<<"problem uses mass fractions"<<std::endl;

        // write caption into output file
        std::ofstream outputFile;
        std::string fileName = name() + "timedata_a_" + GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, SpatialParams, FractureAperture) + ".log";
        outputFile.open(fileName, std::ios::out);
        outputFile << "time  \t | \t  "
                   << "K*gradP_n*n into fracture [kg/s] \t | \t "
                   << "phase transfer [kg/s] into fracture \t | \t "
                   << "advective N2 transfer in water into fracture \t | \t "
                   << "diffusive N2 transfer in water into fracture \t | \t "
                   << "advective water transfer in N2 into fracture \t | \t "
                   << "diffusive water transfer in N2 into fracture \t | \t "
                   << "heat transfer into fracture \t | \t "
                   << "K*gradP*n from fracture [kg/s] \t | \t "
                   << "phase transfer [kg/s] from fracture \t | \t "
                   << "advective N2 transfer in water from fracture \t | \t "
                   << "diffusive N2 transfer in water from fracture \t | \t "
                   << "advective water transfer in N2 from fracture \t | \t "
                   << "diffusive water transfer in N2 from fracture \t | \t "
                   << "heat transfer from fracture \n\n";
        outputFile.close();
    }

    void postTimeStep()
    {
        using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
        static const GlobalPosition fractureOrigin = {-0.5, 0.0};
        const auto apertureString = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, SpatialParams, FractureAperture);

        if (!this->timeManager().episodeWillBeFinished())
            return;

        Scalar kgradpInFracture = 0.0;
        Scalar phaseFluxInFracture = 0.0;
        Scalar advectiveN2FluxInFracture = 0.0;
        Scalar diffusiveN2FluxInFracture = 0.0;
        Scalar advectiveH2OFluxInFracture = 0.0;
        Scalar diffusiveH2OFluxInFracture = 0.0;
        Scalar heatFluxInFracture = 0.0;

        Scalar kgradpFromFracture = 0.0;
        Scalar phaseFluxFromFracture = 0.0;
        Scalar advectiveN2FluxFromFracture = 0.0;
        Scalar diffusiveN2FluxFromFracture = 0.0;
        Scalar advectiveH2OFluxFromFracture = 0.0;
        Scalar diffusiveH2OFluxFromFracture = 0.0;
        Scalar heatFluxFromFracture = 0.0;

        // at a given time we want to plot -(K gradP)*n and - (D gradC)*n along the fracture sides
        std::ofstream lowerFluxPlotFile, upperFluxPlotFile;
        const Scalar time = this->timeManager().time() + this->timeManager().timeStepSize();
        lowerFluxPlotFile.open(name() + "_fluxplotdata_loweredge_a_" + apertureString + "_t_" + std::to_string(int(time)) + ".log", std::ios::out);
        upperFluxPlotFile.open(name() + "_fluxplotdata_upperedge_a_" + apertureString + "_t_" + std::to_string(int(time)) + ".log", std::ios::out);

        // calculate outflow on the right side
        for (const auto& element : elements(this->gridView()))
        {
            couplingManager().setCouplingContext(element);

            auto fvGeometry = localView(this->model().globalFvGeometry());
            fvGeometry.bind(element);

            auto elemVolVars = localView(this->model().curGlobalVolVars());
            elemVolVars.bind(element, fvGeometry, this->model().curSol());

            auto elemFluxVarsCache = localView(this->model().globalFluxVarsCache());
            elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

            // check for scvfs at the lower fracture interface
            for (const auto& scvf : scvfs(fvGeometry))
            {
                const auto pos = scvf.ipGlobal();
                const auto y = pos[1];
                const Scalar yFracture = 0.0;

                if (y > yFracture - eps_ && y < yFracture + eps_ && scvf.unitOuterNormal()[1] > 0.0)
                {
                    FluxVariables fluxVars;
                    fluxVars.init(*this, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

                    auto dummyUpwindTerm = [](const auto& volVars)
                    { return 1; };

                    auto N2PhaseUpwindTerm = [] (const auto& volVars)
                    { return volVars.molarDensity(nPhaseIdx)*volVars.mobility(nPhaseIdx); };

                    auto N2InH2OUpwindTerm = [](const auto& volVars)
                    { return volVars.molarDensity(wPhaseIdx)*volVars.moleFraction(wPhaseIdx, nCompIdx)*volVars.mobility(wPhaseIdx); };

                    auto H2OInN2UpwindTerm = [](const auto& volVars)
                    { return volVars.molarDensity(nPhaseIdx)*volVars.moleFraction(nPhaseIdx, wCompIdx)*volVars.mobility(nPhaseIdx); };

                    const auto KGradP = fluxVars.advectiveFlux(nPhaseIdx, dummyUpwindTerm);
                    const auto q_phase = fluxVars.advectiveFlux(nPhaseIdx, N2PhaseUpwindTerm);
                    const auto q_adv_component = fluxVars.advectiveFlux(wPhaseIdx, N2InH2OUpwindTerm);
                    const auto DGradC = fluxVars.molecularDiffusionFlux(wPhaseIdx, nCompIdx);
                    const auto LGradT = fluxVars.heatConductionFlux();

                    // water transfe
                    const auto q_adv_h2o = fluxVars.advectiveFlux(nPhaseIdx, H2OInN2UpwindTerm);
                    const auto q_diff_h2o = fluxVars.molecularDiffusionFlux(nPhaseIdx, wCompIdx);

                    kgradpInFracture += KGradP;
                    phaseFluxInFracture += q_phase;
                    advectiveN2FluxInFracture += q_adv_component;
                    diffusiveN2FluxInFracture += DGradC;
                    advectiveH2OFluxInFracture += q_adv_h2o;
                    diffusiveH2OFluxInFracture += q_diff_h2o;
                    heatFluxInFracture += LGradT;

                    const Scalar arc_length = (pos-fractureOrigin).two_norm();
                    lowerFluxPlotFile << std::setprecision(8)
                                      << arc_length << '\t'
                                      << KGradP/scvf.area() << '\t'
                                      << q_phase/scvf.area() << '\t'
                                      << q_adv_component/scvf.area() << '\t'
                                      << DGradC/scvf.area() << '\t'
                                      << LGradT/scvf.area() << '\n';
                }

                if (y > yFracture - eps_ && y < yFracture + eps_ && scvf.unitOuterNormal()[1] < 0.0)
                {
                    FluxVariables fluxVars;
                    fluxVars.init(*this, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

                    auto dummyUpwindTerm = [](const auto& volVars)
                    { return 1; };

                    auto N2PhaseUpwindTerm = [] (const auto& volVars)
                    { return volVars.molarDensity(nPhaseIdx)*volVars.mobility(nPhaseIdx); };

                    auto N2InH2OUpwindTerm = [](const auto& volVars)
                    { return volVars.molarDensity(wPhaseIdx)*volVars.moleFraction(wPhaseIdx, nCompIdx)*volVars.mobility(wPhaseIdx); };

                    auto H2OInN2UpwindTerm = [](const auto& volVars)
                    { return volVars.molarDensity(nPhaseIdx)*volVars.moleFraction(nPhaseIdx, wCompIdx)*volVars.mobility(nPhaseIdx); };

                    const auto KGradP = fluxVars.advectiveFlux(nPhaseIdx, dummyUpwindTerm);
                    const auto q_phase = fluxVars.advectiveFlux(nPhaseIdx, N2PhaseUpwindTerm);
                    const auto q_adv_component = fluxVars.advectiveFlux(wPhaseIdx, N2InH2OUpwindTerm);
                    const auto DGradC = fluxVars.molecularDiffusionFlux(wPhaseIdx, nCompIdx);
                    const auto LGradT = fluxVars.heatConductionFlux();

                    // water transfe
                    const auto q_adv_h2o = fluxVars.advectiveFlux(nPhaseIdx, H2OInN2UpwindTerm);
                    const auto q_diff_h2o = fluxVars.molecularDiffusionFlux(nPhaseIdx, wCompIdx);

                    kgradpFromFracture += KGradP;
                    phaseFluxFromFracture += q_phase;
                    advectiveN2FluxFromFracture += q_adv_component;
                    diffusiveN2FluxFromFracture += DGradC;
                    advectiveH2OFluxFromFracture += q_adv_h2o;
                    diffusiveH2OFluxFromFracture += q_diff_h2o;
                    heatFluxFromFracture += LGradT;

                    const Scalar arc_length = (pos-fractureOrigin).two_norm();
                    upperFluxPlotFile << std::setprecision(8)
                                      << arc_length << '\t'
                                      << KGradP/scvf.area() << '\t'
                                      << q_phase/scvf.area() << '\t'
                                      << q_adv_component/scvf.area() << '\t'
                                      << DGradC/scvf.area() << '\t'
                                      << LGradT/scvf.area() << '\n';
                }
            }
        }

        // write output to file
        std::ofstream outputFile;
        std::string fileName = name() + "timedata_a_" + GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, SpatialParams, FractureAperture) + ".log";
        outputFile.open(fileName, std::ios::out | std::ios::app);
        outputFile << std::setprecision(8)
                   << this->timeManager().time() + this->timeManager().timeStepSize() << "\t\t"
                   << kgradpInFracture << '\t'
                   << phaseFluxInFracture << '\t'
                   << advectiveN2FluxInFracture << '\t'
                   << diffusiveN2FluxInFracture << '\t'
                   << advectiveH2OFluxInFracture << '\t'
                   << diffusiveH2OFluxInFracture << '\t'
                   << heatFluxInFracture << '\t'
                   << kgradpFromFracture << '\t'
                   << phaseFluxFromFracture << '\t'
                   << advectiveN2FluxFromFracture << '\t'
                   << diffusiveN2FluxFromFracture << '\t'
                   << advectiveH2OFluxFromFracture << '\t'
                   << diffusiveH2OFluxFromFracture << '\t'
                   << heatFluxFromFracture << '\n';
        outputFile.close();
    }

    /*!
     * \brief Returns the initial phase state for a control volume.
     *
     * \param scv The sub-control volume
     */
    template<class Entity>
    int initialPhasePresence(const Entity& entity) const
    {
        if (isInLens(entity.geometry().center()))
            return Indices::bothPhases;
        else
            return Indices::wPhaseOnly;
    }

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }

    /*!
     * \brief Return the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; }

    /*!
     * \brief Return the sources within the domain.
     */
    PrimaryVariables sourceAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();

        // set Dirichlet at outlet
        if (isOnInlet(globalPos) || isOnOutlet(globalPos))
            values.setAllDirichlet();

        return values;
    }

    /*!
     * \brief Specifies if a given intersection is on an interior boundary
     */
    bool isInteriorBoundary(const Element& element, const Intersection& is) const
    { return couplingManager().isInteriorBoundary(element, is); }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume (isothermal case).
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        static const Scalar deltaP = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, OverPressure);

        auto values = initialAtPos(globalPos);
        if (isOnInlet(globalPos))
            values[pressureIdx] += deltaP;

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     */
    PrimaryVariables neumann(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolumeFace& scvf) const
    { return PrimaryVariables(0.0); }

    /*!
     * \brief Evaluate the initial value for a control volume (isothermal case)
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        static const Scalar initLensSaturation = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, InitialLensSaturation);
        static const Scalar initLensTemperature = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, InitialLensTemperature);

        PrimaryVariables priVars(0.0);
        priVars[pressureIdx] = 1e5; // initial condition for the pressure
        priVars[Indices::temperatureIdx] = temperature();
        if (isInLens(globalPos))
        {
            priVars[switchIdx] = initLensSaturation;
            priVars[Indices::temperatureIdx] = initLensTemperature;
        }
        return priVars;
    }

    //! Set the coupling manager
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    CouplingManager& couplingManager()
    { return *couplingManager_; }

    bool isInLens(const GlobalPosition& globalPos) const
    { return globalPos[0] < -0.2 + eps_ && globalPos[1] > -0.3 - eps_ && globalPos[1] < -0.1 + eps_; }

    bool isOnInlet(const GlobalPosition& globalPos) const
    { return globalPos[0] < this->bBoxMin()[0] + eps_ && globalPos[1] < -0.4; }

    bool isOnOutlet(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->bBoxMax()[0] - eps_ && globalPos[1] > 0.4; }

private:
    std::string name_;
    Scalar eps_;
    std::shared_ptr<CouplingManager> couplingManager_;
};
} //end namespace

#endif
