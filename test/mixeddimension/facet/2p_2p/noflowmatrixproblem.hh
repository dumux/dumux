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
#ifndef DUMUX_NOFLOW_MATRIX_PROBLEM_HH
#define DUMUX_NOFLOW_MATRIX_PROBLEM_HH

#include <dumux/mixeddimension/facet/mpfa/properties.hh>
#include <dumux/mixeddimension/subproblemproperties.hh>

#include <dumux/porousmediumflow/2p/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>

#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/dnapl.hh>
#include "matrixspatialparams.hh"

#define COMPLEXPHYSICS 0

namespace Dumux
{
template <class TypeTag>
class NoFlowDomainMatrixProblem;

namespace Properties
{
NEW_TYPE_TAG(NoFlowDomainMatrixProblem, INHERITS_FROM(FacetCouplingBulkMpfaModel, TwoP, MatrixSpatialParams));

// Set fluid configuration
#if COMPLEXPHYSICS
SET_TYPE_PROP(NoFlowDomainMatrixProblem, FluidSystem, FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), true>);
#else
SET_TYPE_PROP(NoFlowDomainMatrixProblem, FluidSystem, FluidSystems::TwoPImmiscible<typename GET_PROP_TYPE(TypeTag, Scalar),
                                                                                   FluidSystems::LiquidPhase<typename GET_PROP_TYPE(TypeTag, Scalar), SimpleH2O<typename GET_PROP_TYPE(TypeTag, Scalar)>>,
                                                                                   FluidSystems::LiquidPhase<typename GET_PROP_TYPE(TypeTag, Scalar), DNAPL<typename GET_PROP_TYPE(TypeTag, Scalar)>>>);
#endif

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
        pwIdx = Indices::pwIdx,
        snIdx = Indices::snIdx,

        // equation indices
        contiWEqIdx = Indices::contiWEqIdx,
        contiNEqIdx = Indices::contiNEqIdx,

        // phase indices
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
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
        std::string fileName = name() + "timedata_ld_a_" + GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, SpatialParams, FractureAperture) + ".log";
        outputFile.open(fileName, std::ios::out);
        outputFile << "t  \t | \t  advective phase transfer into fracture [kg/s] \t | \t "
                   << "advective component transfer [kg/s] into fracture \t | \t "
                   << "diffusive component transfer [kg/s] into fracture \t | \t "
                   << "concentration in producer\t | \t "
                   << "advective phase transfer from fracture [kg/s] \t | \t "
                   << "advective component transfer [kg/s] from fracture \t | \t "
                   << "diffusive component transfer [kg/s] from fracture \t | \t "
                   << "advective phase inflow [kg/s] through boundary \t | \t "
                   << "advective component inflow [kg/s] through boundary \t | \t "
                   << "diffusive component inflow [kg/s] through boundary \t | \t "
                   << "advective phase outflow [kg/s] through boundary \t | \t "
                   << "advective component outflow [kg/s] through boundary \t | \t "
                   << "diffusive component outflow [kg/s] through boundary \n\n";
        outputFile.close();
    }

    void postTimeStep() const
    {
        using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);

        Scalar phaseFluxInFracture = 0.0;
        Scalar advCompFluxInFracture = 0.0;
        Scalar diffCompFluxInFracture = 0.0;
        Scalar x_outlet;
        Scalar phaseFluxFromFracture = 0.0;
        Scalar advCompFluxFromFracture = 0.0;
        Scalar diffCompFluxFromFracture = 0.0;
        Scalar boundaryInfluxPhase = 0.0;
        Scalar boundaryInfluxAdvective = 0.0;
        Scalar boundaryInfluxDiffusive = 0.0;
        Scalar boundaryOutfluxPhase = 0.0;
        Scalar boundaryOutfluxAdvective = 0.0;
        Scalar boundaryOutfluxDiffusive = 0.0;

        // at a given time we want to plot -(K gradP)*n and - (D gradC)*n along the fracture sides
        std::ofstream lowerFluxPlotFile, upperFluxPlotFile;
        const Scalar dt = this->timeManager().timeStepSize();
        const Scalar time = this->timeManager().time();
        const bool plotFlux = this->timeManager().willBeFinished();

        if (plotFlux)
        {
            const auto a = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, SpatialParams, FractureAperture);
            const auto t = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, TimeForFractureEdgePlots);
            lowerFluxPlotFile.open(name() + "fluxplotdata_ld_loweredge_a_" + a + "_t_" + t + ".log", std::ios::out);
            upperFluxPlotFile.open(name() + "fluxplotdata_ld_upperedge_a_" + a + "_t_" + t + ".log", std::ios::out);
        }

        // calculate outflow on the right side
        for (const auto& element : elements(this->gridView()))
        {
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
                const auto yFracture = 0.0;

                if (y > yFracture - eps_ && y < yFracture + eps_ && scvf.unitOuterNormal()[1] > 0.0)
                {
                    FluxVariables fluxVars;
                    fluxVars.init(*this, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

                    auto dummyUpwindTerm = [](const auto& volVars)
                    { return 1; };

                    auto compUpwindTerm = [](const auto& volVars)
                    { return volVars.density(nPhaseIdx)*volVars.saturation(nPhaseIdx)*volVars.mobility(nPhaseIdx); };

                    const auto KGradP = fluxVars.advectiveFlux(nPhaseIdx, dummyUpwindTerm);
                    const auto DGradX = 0.0;

                    phaseFluxInFracture += KGradP;
                    advCompFluxInFracture += fluxVars.advectiveFlux(nPhaseIdx, compUpwindTerm);
                    diffCompFluxInFracture += 0.0;

                    if (plotFlux)
                    {
                        const GlobalPosition interceptFracture = [&] () { GlobalPosition tmp(0.0); tmp[0] = this->bBoxMin()[0]; return tmp; } ();
                        const Scalar arc_length = (pos-interceptFracture).two_norm();
                        lowerFluxPlotFile << std::setprecision(8)
                                          << arc_length << '\t'
                                          << KGradP/scvf.area() << '\t'
                                          << DGradX/scvf.area() << '\n';
                    }
                }

                if (y > yFracture - eps_ && y < yFracture + eps_ && scvf.unitOuterNormal()[1] < 0.0)
                {
                    FluxVariables fluxVars;
                    fluxVars.init(*this, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

                    auto dummyUpwindTerm = [](const auto& volVars)
                    { return 1; };

                    auto compUpwindTerm = [](const auto& volVars)
                    { return volVars.density(nPhaseIdx)*volVars.saturation(nPhaseIdx)*volVars.mobility(nPhaseIdx); };

                    const auto KGradP = fluxVars.advectiveFlux(nPhaseIdx, dummyUpwindTerm);
                    const auto DGradX = 0.0;

                    phaseFluxFromFracture += KGradP;
                    advCompFluxFromFracture += fluxVars.advectiveFlux(nPhaseIdx, compUpwindTerm);
                    diffCompFluxFromFracture += 0.0;

                    if (plotFlux)
                    {
                        const GlobalPosition interceptFracture = [&] () { GlobalPosition tmp(0.0); tmp[0] = this->bBoxMin()[0]; return tmp; } ();
                        const Scalar arc_length = (pos-interceptFracture).two_norm();
                        upperFluxPlotFile << std::setprecision(8)
                                          << arc_length << '\t'
                                          << KGradP/scvf.area() << '\t'
                                          << DGradX/scvf.area() << '\n';
                    }
                }

                if (isOnOutlet(pos))
                {
                    FluxVariables fluxVars;
                    fluxVars.init(*this, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

                    auto phaseUpwindTerm = [](const auto& volVars)
                    { return 1; /*volVars.density(0)*volVars.mobility(0);*/ };

                    auto compUpwindTerm = [](const auto& volVars)
                    { return volVars.density(nPhaseIdx)*volVars.saturation(nPhaseIdx)*volVars.mobility(nPhaseIdx); };

                    boundaryOutfluxPhase += fluxVars.advectiveFlux(nPhaseIdx, phaseUpwindTerm);
                    boundaryOutfluxAdvective += fluxVars.advectiveFlux(nPhaseIdx, compUpwindTerm);
                    boundaryOutfluxDiffusive += 0.0;

                    x_outlet = elemVolVars[scvf.insideScvIdx()].saturation(nPhaseIdx);
                }

                if (isOnInlet(pos))
                {
                    FluxVariables fluxVars;
                    fluxVars.init(*this, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

                    auto phaseUpwindTerm = [](const auto& volVars)
                    { return 1; /*volVars.density(0)*volVars.mobility(0);*/ };

                    auto compUpwindTerm = [](const auto& volVars)
                    { return volVars.density(nPhaseIdx)*volVars.saturation(nPhaseIdx)*volVars.mobility(nPhaseIdx); };

                    boundaryInfluxPhase += fluxVars.advectiveFlux(nPhaseIdx, phaseUpwindTerm);
                    boundaryInfluxAdvective += fluxVars.advectiveFlux(nPhaseIdx, compUpwindTerm);
                    boundaryInfluxDiffusive += 0.0;
                }
            }
        }

        // write output to file
        std::ofstream outputFile;
        std::string fileName = name() + "timedata_ld_a_" + GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, SpatialParams, FractureAperture) + ".log";
        outputFile.open(fileName, std::ios::out | std::ios::app);
        outputFile << std::setprecision(8)
                   << this->timeManager().time() + this->timeManager().timeStepSize() << "\t\t"
                   << phaseFluxInFracture << '\t'
                   << advCompFluxInFracture << '\t'
                   << diffCompFluxInFracture << '\t'
                   << x_outlet << '\t'
                   << phaseFluxFromFracture << '\t'
                   << advCompFluxFromFracture << '\t'
                   << diffCompFluxFromFracture << '\t'
                   << boundaryInfluxPhase << '\t'
                   << boundaryInfluxAdvective << '\t'
                   << boundaryInfluxDiffusive << '\t'
                   << boundaryOutfluxPhase << '\t'
                   << boundaryOutfluxAdvective << '\t'
                   << boundaryOutfluxDiffusive << '\n';
        outputFile.close();
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
        {
            values[pwIdx] += deltaP;
            values[snIdx] = 1.0;
        }

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
        PrimaryVariables priVars(0.0);
        priVars[pwIdx] = 1e5; // initial condition for the pressure
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

    bool isOnOutlet(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->bBoxMax()[0] - eps_ && globalPos[1] > 0.25; }

    bool isOnInlet(const GlobalPosition& globalPos) const
    { return globalPos[0] < this->bBoxMin()[0] + eps_ && globalPos[1] < -0.25; }

private:
    std::string name_;
    Scalar eps_;
    std::shared_ptr<CouplingManager> couplingManager_;
};
} //end namespace

#endif
