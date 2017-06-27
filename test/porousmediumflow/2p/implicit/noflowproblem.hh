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
#ifndef DUMUX_2P_NOFLOW_PROBLEM_HH
#define DUMUX_2P_NOFLOW_PROBLEM_HH

#include <dumux/implicit/box/properties.hh>
#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/implicit/cellcentered/mpfa/properties.hh>
#include <dumux/porousmediumflow/2p/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>

#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/dnapl.hh>
#include "noflowspatialparams.hh"

#define COMPLEXPHYSICS 0

namespace Dumux
{

//! Forward declaration of the problem class
template <class TypeTag>
class NoFlowProblem;

namespace Properties
{
NEW_TYPE_TAG(NoFlowProblem, INHERITS_FROM(TwoP, TwoPSpatialParams));

NEW_TYPE_TAG(NoFlowCCMpfaProblem, INHERITS_FROM(CCMpfaModel, NoFlowProblem));
SET_BOOL_PROP(NoFlowCCMpfaProblem, SolutionDependentAdvection, false);
SET_BOOL_PROP(NoFlowCCMpfaProblem, EnableGlobalFVGeometryCache, true);
SET_BOOL_PROP(NoFlowCCMpfaProblem, EnableGlobalVolumeVariablesCache, true);
SET_BOOL_PROP(NoFlowCCMpfaProblem, EnableGlobalFluxVariablesCache, true);

// Set the grid type
SET_TYPE_PROP(NoFlowProblem, Grid, Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);

// Set the problem property
SET_TYPE_PROP(NoFlowProblem, Problem, NoFlowProblem<TypeTag>);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(NoFlowProblem, UseMoles, true);

// Enable velocity output
SET_BOOL_PROP(NoFlowProblem, VtkAddVelocity, true);

// Disable gravity
SET_BOOL_PROP(NoFlowProblem, ProblemEnableGravity, false);

// Set fluid configuration
#if COMPLEXPHYSICS
SET_TYPE_PROP(NoFlowProblem, FluidSystem, FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), true>);
#else
SET_TYPE_PROP(NoFlowProblem, FluidSystem, FluidSystems::TwoPImmiscible<typename GET_PROP_TYPE(TypeTag, Scalar),
                                                                       FluidSystems::LiquidPhase<typename GET_PROP_TYPE(TypeTag, Scalar), SimpleH2O<typename GET_PROP_TYPE(TypeTag, Scalar)>>,
                                                                       FluidSystems::LiquidPhase<typename GET_PROP_TYPE(TypeTag, Scalar), DNAPL<typename GET_PROP_TYPE(TypeTag, Scalar)>>>);
#endif
}

template <class TypeTag>
class NoFlowProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    using ParentType = ImplicitPorousMediaProblem<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);

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

    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    NoFlowProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        //initialize fluid system
        FluidSystem::init();

        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);
        aperture_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, FractureAperture);

        // stating in the console whether mole or mass fractions are used
        if(useMoles)
            std::cout<<"problem uses mole fractions"<<std::endl;
        else
            std::cout<<"problem uses mass fractions"<<std::endl;

        // write caption into output file
        std::ofstream outputFile;
        const auto apertureString = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, SpatialParams, FractureAperture);
        outputFile.open(name() + "_timedata_ed_a_" + apertureString + ".log", std::ios::out);
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
                   << "diffusive component outflow [kg/s] through boundary \n\n";;
        outputFile.close();
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
    const std::string& name() const
    {
        return name_;
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
            lowerFluxPlotFile.open(name() + "_fluxplotdata_ed_loweredge_a_" + a + "_t_" + t + ".log", std::ios::out);
            upperFluxPlotFile.open(name() + "_fluxplotdata_ed_upperedge_a_" + a + "_t_" + t + ".log", std::ios::out);
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
                const auto x = pos[0];
                const auto y = pos[1];
                const auto yFractureLow = this->spatialParams().yLowerFracture(x);
                const auto yFractureUp = this->spatialParams().yUpperFracture(x);

                if (y > yFractureLow - eps_ && y < yFractureLow + eps_ && scvf.unitOuterNormal()[1] > 0.0)
                {
                    FluxVariables fluxVars;
                    fluxVars.init(*this, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

                    auto dummyUpwindTerm = [](const auto& volVars)
                    { return 1; };

                    auto compUpwindTerm = [](const auto& volVars)
                    { return volVars.density(nPhaseIdx)*volVars.saturation(nPhaseIdx)*volVars.mobility(nPhaseIdx); };

                    const auto KGradP = fluxVars.advectiveFlux(0, dummyUpwindTerm);
                    const auto DGradX = 0.0;

                    phaseFluxInFracture += KGradP;
                    advCompFluxInFracture += fluxVars.advectiveFlux(0, compUpwindTerm);
                    diffCompFluxInFracture += 0.0;

                    if (plotFlux)
                    {
                        const GlobalPosition interceptFracture = [&] () { GlobalPosition tmp(0.0); tmp[1] = this->spatialParams().yLowerFracture(0.0); return tmp; } ();
                        const Scalar arc_length = (pos-interceptFracture).two_norm();
                        lowerFluxPlotFile << std::setprecision(8)
                                          << arc_length << '\t'
                                          << KGradP/scvf.area() << '\t'
                                          << DGradX/scvf.area() << '\n';
                    }
                }

                if (y > yFractureUp - eps_ && y < yFractureUp + eps_ && scvf.unitOuterNormal()[1] < 0.0)
                {
                    FluxVariables fluxVars;
                    fluxVars.init(*this, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

                    auto dummyUpwindTerm = [](const auto& volVars)
                    { return 1; };

                    auto compUpwindTerm = [](const auto& volVars)
                    { return volVars.density(nPhaseIdx)*volVars.saturation(nPhaseIdx)*volVars.mobility(nPhaseIdx); };

                    const auto KGradP = fluxVars.advectiveFlux(0, dummyUpwindTerm);
                    const auto DGradX = 0.0;

                    phaseFluxFromFracture += KGradP;
                    advCompFluxFromFracture += fluxVars.advectiveFlux(nPhaseIdx, compUpwindTerm);
                    diffCompFluxFromFracture += 0.0;

                    if (plotFlux)
                    {
                        const GlobalPosition interceptFracture = [&] () { GlobalPosition tmp(0.0); tmp[1] = this->spatialParams().yUpperFracture(0.0); return tmp; } ();
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
        const auto apertureString = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, SpatialParams, FractureAperture);
        outputFile.open(name() + "_timedata_ed_a_" + apertureString + ".log", std::ios::out | std::ios::app);
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

#if !NONISOTHERMAL
    /*!
     * \brief Returns the temperature within the domain [K].
     *
     * This problem assumes a temperature of 20 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 20; } // in [K]

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
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();

        // set Dirichlet at inlet and producer
        if (isOnInlet(globalPos) || isOnOutlet(globalPos))
            values.setAllDirichlet();

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        static const Scalar deltaP = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, OverPressure);
        auto values = initial_(globalPos);

        if (isOnInlet(globalPos))
        {
            values[pwIdx] += deltaP;
            values[snIdx] = 1.0;
        }
        return values;
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
    PrimaryVariables neumannAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

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
    PrimaryVariables sourceAtPos(const GlobalPosition &globalPos) const
    { return PrimaryVariables(0.0); }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return initial_(globalPos); }

    /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    template<class OutputModule>
    void addVtkOutputFields(OutputModule& outputModule) const
    {
        // create the required scalar fields
        auto& isFractureElement = outputModule.createScalarField("isFractureElement", 0);
        auto& touchesLowerFractureInterface = outputModule.createScalarField("touchesLowerFractureInterface", 0);

        for (const auto& element : elements(this->gridView()))
        {
            auto fvGeometry = localView(this->model().globalFvGeometry());
            fvGeometry.bind(element);

            const auto eIdxGlobal = this->elementMapper().index(element);
            isFractureElement[eIdxGlobal] = this->spatialParams().isFractureElement(element);

            for (const auto& scvf : scvfs(fvGeometry))
            {
                const auto pos = scvf.ipGlobal();
                const auto x = pos[0];
                const auto y = pos[1];
                const auto yFracture = this->spatialParams().yLowerFracture(x);

                if (y > yFracture - eps_ && y < yFracture + eps_ && scvf.unitOuterNormal()[1] > 0.0)
                    touchesLowerFractureInterface[eIdxGlobal] = true;
            }
        }
    }

    Scalar aperture() const
    { return aperture_; }

    bool isOnOutlet(const GlobalPosition& globalPos) const
    { return globalPos[1] > this->bBoxMax()[1] - eps_; }

    bool isOnInlet(const GlobalPosition& globalPos) const
    { return globalPos[1] < this->bBoxMin()[1] + eps_; }

    // \}

private:
    // the internal method for the initial condition
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars;
        priVars[pwIdx] = 1e5; // initial condition for the pressure
        priVars[snIdx] = 0.0;  // initial condition for the N2 molefraction
        return priVars;
    }

    static constexpr Scalar eps_ = 1e-6;
    Scalar aperture_;
    std::string name_;
};

} //end namespace Dumux

#endif
