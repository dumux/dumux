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
 * \brief Non-isothermal gas injection problem where a gas (e.g. air)
 *        is injected into a fully water saturated medium.
 */
#ifndef DUMUX_WINDTUNNEL_DARCY_SUBPROBLEM_HH
#define DUMUX_WINDTUNNEL_DARCY_SUBPROBLEM_HH

#include <dumux/io/gnuplotinterface.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysomerton.hh>
#include <dumux/porousmediumflow/2p2c/implicit/model.hh>

#include <appl/staggeredgrid/multidomain/navierstokes2ctdarcy2p2ct/properties.hh> // TODO

#include <dumux/io/readwritedatafile.hh> // TODO

#include "windtunnel_spatialparams.hh"

#define ISOTHERMAL 0

// TODO necessary? (used in darcytestproblem in stokesdarcy1p)
//#include <dumux/implicit/cellcentered/tpfa/properties.hh>
//#include <dumux/porousmediumflow/implicit/problem.hh>
//#include <dumux/linear/seqsolverbackend.hh>
//#include <dumux/material/components/simpleh2o.hh>
//#include <dumux/material/components/constant.hh>
//#include <dumux/material/fluidsystems/liquidphase.hh>
//// coupling-specific includes
//#include <dumux/multidomain/subproblemproperties.hh>

namespace Dumux
{
template <class TypeTag>
class WindtunnelDarcySubProblem;

namespace Properties
{
// TODO necessary? (stokesdarcy1p, replace OneP)
//NEW_TYPE_TAG(DarcyTestProblem, INHERITS_FROM(CCTpfaModel, OneP, WindTunnelSpatialParams));

// Set the problem property // TODO necessary? (stokesdarcy1p) replace --> WindtunnelDarcySubProblem?
//SET_TYPE_PROP(DarcySubProblem, Problem, Dumux::DarcySubProblem<TypeTag>);

// Set the spatial parameters // TODO necessary? (stokesdarcy1p)
//SET_TYPE_PROP(DarcyTestProblem, SpatialParams, OnePSpatialParams<TypeTag>);

// Enable gravity
SET_BOOL_PROP(DarcySubProblem, ProblemEnableGravity, true);

// choose pn and Sw as primary variables
SET_INT_PROP(DarcySubProblem, Formulation, TwoPTwoCFormulation::pnsw);

// the gas component balance (air) is replaced by the total mass balance
SET_INT_PROP(DarcySubProblem, ReplaceCompEqIdx, GET_PROP_TYPE(TypeTag, Indices)::contiNEqIdx);

// Use Kelvin equation to adapt the saturation vapor pressure
SET_BOOL_PROP(DarcySubProblem, UseKelvinEquation, true);

// Somerton is used as model to compute the effective thermal heat conductivity
SET_TYPE_PROP(DarcySubProblem, ThermalConductivityModel,
              ThermalConductivitySomerton<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// Depth in third dimension
NEW_PROP_TAG(GridExtrusionFactor);
SET_SCALAR_PROP(DarcySubProblem, GridExtrusionFactor, 1.0);

// Live plot of the evaporation rates
NEW_PROP_TAG(OutputPlotEvaporationRate);
SET_BOOL_PROP(DarcySubProblem, OutputPlotEvaporationRate, false);

// Set the output frequency
NEW_PROP_TAG(OutputFreqFluxOutput);
SET_INT_PROP(DarcySubProblem, OutputFreqFluxOutput, 100000);

// Set the output frequency
NEW_PROP_TAG(TimeManagerAbortIfEvapRateIsConstant);
SET_BOOL_PROP(DarcySubProblem, TimeManagerAbortIfEvapRateIsConstant, false);

// TODO necessary? (stokesdarcy1p)
SET_BOOL_PROP(DarcySubProblem, EnableGlobalFVGeometryCache, true);
NEW_PROP_TAG(GlobalProblemTypeTag);
NEW_PROP_TAG(CouplingManager);
}

/*!
 * \ingroup TwoPTwoCModel
 * \ingroup ImplicitTestProblems
 * \brief Non-isothermal gas injection problem where a gas (e.g. air)
 *        is injected into a fully water saturated medium. During
 *        buoyancy driven upward migration the gas passes a high
 *        temperature area.
 *
 */
template <class TypeTag >
class WindtunnelDarcySubProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    using ParentType = ImplicitPorousMediaProblem<TypeTag>;

//    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainTypeTag) MultiDomainTypeTag;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
//    typedef typename GET_PROP_TYPE(MultiDomainTypeTag, SubDomainGridView) SubDomainGridView;

    // copy some indices for convenience
//    typedef typename GET_PROP_TYPE(MultiDomainTypeTag, Indices) MultiDomainIndices;
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    // the equation indices
    enum { contiTotalMassIdx = Indices::contiNEqIdx,
           contiWEqIdx = Indices::contiWEqIdx,
           energyEqIdx = Indices::energyEqIdx };
    // the indices of the primary variables
    enum { pressureIdx = Indices::pressureIdx,
           switchIdx = Indices::switchIdx,
           temperatureIdx = Indices::temperatureIdx };
    // the indices for the phase presence
    enum { wPhaseOnly = Indices::wPhaseOnly,
           nPhaseOnly = Indices::nPhaseOnly,
           bothPhases = Indices::bothPhases };
    // the phase and component indices
    enum { wPhaseIdx = Indices::wPhaseIdx,
           nPhaseIdx = Indices::nPhaseIdx,
           wCompIdx = Indices::wCompIdx };
    // Grid and world dimension
    enum { dim = GridView::dimension,
           dimWorld = GridView::dimensionworld };

    enum { dofCodim = 0 };

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);

    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<dim>::Entity;
    using Intersection = typename GridView::Intersection;

    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry) ;
    using LocalJacobian = typename GET_PROP_TYPE(TypeTag, LocalJacobian);
    using LocalResidual = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using ThermalConductivityModel = typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel);

    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;

    // property that defines whether mole or mass fractions are used
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
//    static const unsigned int darcySubDomainIdx = MultiDomainIndices::darcySubDomainIdx;

    using GlobalTypeTag = typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag);
    using CouplingManager = typename GET_PROP_TYPE(GlobalTypeTag, CouplingManager);

public:
    /*!
     * \brief The constructor.
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    WindtunnelDarcySubProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
        , gravity_(0.0)
        , gridView_(gridView)
    {
        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
            gravity_[1] = -9.81;

        std::vector<Scalar> positions0 = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::vector<Scalar>, Grid, Positions0);
        std::vector<Scalar> positions1 = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::vector<Scalar>, Grid, Positions1);

        bBoxMin_[0] = std::max(positions0.front(),GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, NoDarcyX1));
        bBoxMax_[0] = std::min(positions0.back(),GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, NoDarcyX2));
#if DUMUX_MULTIDOMAIN_DIM > 2 // TODO
        bBoxMin_[2] = std::max(positions0.front(),GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, NoDarcyZ1));
        bBoxMax_[2] = std::min(positions0.back(),GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, NoDarcyZ2));
#endif
        extrusionFactor_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, ExtrusionFactor);
        if (extrusionFactor_ < 0.999999 || extrusionFactor_ > 1.000001)
            DUNE_THROW(Dune::NotImplemented, "An extrusionFactor_ != 1, has to be implemented for all accumulate terms in the staggered balance equations.");

        interfacePos_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, InterfacePosY);

        velocity_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefVelocity);

        pressure_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PorousMedium, RefPressure);
        switch_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PorousMedium, RefWaterSaturation);
        initialPhasePresence_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, PorousMedium, InitialPhasePresence);
        temperature_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PorousMedium, RefTemperature);
        solDependentEnergySource_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, PorousMedium, SolDependentEnergySource);
        solDependentEnergyWalls_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, PorousMedium, SolDependentEnergyWalls);

        try {
            temperatureDataFilePM_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, PorousMedium, TemperatureDataFile);
            useTemperatureDataFilePM_ = true;
            readData(temperatureDataFilePM_, temperatureDataPM_);
        }
        catch (...) {
            useTemperatureDataFilePM_ = false;
        }

        freqFluxOutput_ = GET_PARAM_FROM_GROUP(TypeTag, int, Output, FreqFluxOutput);

        storageLastTimestep_ = Scalar(0);
        storageChangeLastTimestep_ = Scalar(0);
        lastMassOutputTime_ = Scalar(0);

        std::string storageFile = name() + "-storage.csv";
        outfile.open(storageFile);
        outfile << "Time[s]" << ";"
                << "TotalMassChange[kg/(s*mDepth)]" << ";"
                << "WaterMassChange[kg/(s*mDepth))]" << ";"
                << "IntEnergyChange[J/(m^3*s*mDepth)]" << ";"
                << "WaterMass[kg/mDepth]" << ";"
                << "WaterMassLoss[kg/mDepth]" << ";"
                << "EvaporationRate[mm/s]" << ";"
                << "MoistureContent[-]"
                << std::endl;
    }

    //! \brief The destructor
    ~WindtunnelDarcySubProblem()
    {
        outfile.close();
    }

    /*!
     * \name Simulation steering
     */
    // \{

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string name() const
    { return GET_RUNTIME_PARAM(TypeTag, std::string, Output.Name) + "_staggered"; } // TODO necessary?

    /*!
     * \brief Called by the Dumux::TimeManager in order to
     *        initialize the problem.
     */
    void init()
    {
        ParentType::init();
        globalStorage(storageLastTimestep_);
        initialWaterContent_ = storageLastTimestep_[contiWEqIdx];
    }

    // suppress output from DuMuX
    bool shouldWriteOutput() const
    {
        return false;
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param values The boundary types for the conservation equations
     * \param globalPos The position for which the bc type should be evaluated
     */
    void boundaryTypesAtPos(BoundaryTypes &values,
                            const GlobalPosition &globalPos) const
    {
        values.setAllNeumann();
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The position for which the bc type should be evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        initialAtPos(values, globalPos);
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
     *
     * The units must be according to either using mole or mass fractions. (mole/(m^2*s) or kg/(m^2*s))
     */
    void neumannAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        values = 0;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{


    //! \copydoc Dumux::ImplicitProblem::solDependentSource()
    void solDependentSource(PrimaryVariables &values,
                            const Element &element,
                            const FVElementGeometry &fvGeometry,
                            const int scvIdx,
                            const ElementVolumeVariables &elemVolVars) const
    {
        values = 0.;

        if (!solDependentEnergyWalls_ && !solDependentEnergySource_)
            return;

        // cell center global
        const Dune::FieldVector<Scalar, dim>& insideCellCenterLocal =
            Dune::ReferenceElements<Scalar, dim>::general(element.geometry().type()).position(0, 0);
        Dune::FieldVector<Scalar, dim> insideCellCenterGlobal = element.geometry().global(insideCellCenterLocal);

        // assume thermal conduction through the plexiglass box
        Scalar plexiglassThermalConductivity = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PorousMedium, PlexiglassThermalConductivity);
        Scalar plexiglassThickness = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PorousMedium, PlexiglassThickness);
        Scalar thermalConductivityInside = ThermalConductivityModel::effectiveThermalConductivity(
                                              elemVolVars[scvIdx].saturation(wPhaseIdx),
                                              elemVolVars[scvIdx].fluidThermalConductivity(wPhaseIdx),
                                              elemVolVars[scvIdx].fluidThermalConductivity(nPhaseIdx),
                                              this->spatialParams().solidThermalConductivityAtPos(insideCellCenterGlobal),
                                              this->spatialParams().porosityAtPos(insideCellCenterGlobal),
                                              this->spatialParams().solidDensityAtPos(insideCellCenterGlobal));
        Scalar temperatureInside = elemVolVars[scvIdx].temperature();
        Scalar temperatureRef = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PorousMedium, SolDependentEnergyTemperature);
        if (useTemperatureDataFilePM_)
            temperatureRef = evaluateData(temperatureDataPM_, this->timeManager().time(), this->timeManager().time()+this->timeManager().timeStepSize());

        if (solDependentEnergyWalls_)
        {
            typedef typename GridView::IntersectionIterator IntersectionIterator;
            for (IntersectionIterator is = gridView_.ibegin(element);
                is != gridView_.iend(element); ++is)
            {
                // face center global
                const Dune::FieldVector<Scalar, dim-1>& faceCenterLocal =
                    Dune::ReferenceElements<Scalar, dim-1>::general(is->geometry().type()).position(0, 0);
                Dune::FieldVector<Scalar, dim> faceCenterGlobal = is->geometry().global(faceCenterLocal);

                Dune::FieldVector<Scalar, dim> distance(faceCenterGlobal);
                distance -= insideCellCenterGlobal;
                Scalar harmonicThermalConductivityFactor = plexiglassThermalConductivity * thermalConductivityInside
                                                           / (plexiglassThermalConductivity * distance.two_norm()
                                                              + thermalConductivityInside * plexiglassThickness);

                if (is->boundary()
                    || (gridView_.indexSet().contains(darcySubDomainIdx, is->inside())
                        != gridView_.indexSet().contains(darcySubDomainIdx, is->outside())))
                {
                    Scalar area = is->geometry().integrationElement(faceCenterLocal);
                    Scalar volume = is->inside().geometry().volume();

                    // NOTE: - sign switches compared to DUMUX Neumann boundary conditions
                    //       - account for additional thickness with distance.two_norm(), dofs are not located at the boundary
                    // NOTE: - discussed and confirmed with Martin S. 2017-09-20
                    values[energyEqIdx] -= harmonicThermalConductivityFactor
                                           * (temperatureInside - temperatureRef)
                                           * area
                                           / volume; // account for the fact that it is treated in the storage term
                }
            }
        }

        if(solDependentEnergySource_ && dim < 3)
        {
            DUNE_THROW(Dune::NotImplemented, "The solDependentEnergySource_ is not working, the extrusionFactor has to be implemented for all accumulate terms in the balance equations");
            Scalar extrusionFactor = this->extrusionFactorAtPos(insideCellCenterGlobal);
            // NOTE: - account for additional thickness with extrusionFactorForSourceTerm_ / 2.0, dofs are not located at the boundary
            // NOTE: - discussed and confirmed with Martin S. 2017-09-20
            Scalar harmonicThermalConductivityFactor = plexiglassThermalConductivity * thermalConductivityInside
                                                        / (plexiglassThermalConductivity * extrusionFactor / 2.0
                                                          + thermalConductivityInside * plexiglassThickness);
            values[energyEqIdx] -= 2.0 * harmonicThermalConductivityFactor
                                   * (temperatureInside - temperatureRef)
                                   / extrusionFactor;
        }
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The initial values for the primary variables
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    void initialAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        values[pressureIdx] = pressure_;
        values[switchIdx] = switch_;
        if (GET_PROP_VALUE(TypeTag, Formulation) == TwoPTwoCFormulation::pwsn
            && initialPhasePresenceAtPos(globalPos) == Indices::bothPhases)
        {
            values[switchIdx] = 1.0 - switch_;
        }
        values[temperatureIdx] = temperature_;
    }

    /*!
     * \brief Return the initial phase state inside a control volume.
     *
     * \param globalPos The global position
     */
    int initialPhasePresenceAtPos(const GlobalPosition &globalPos) const
    {
        if (std::strcmp(initialPhasePresence_.c_str(), "bothPhases") == 0)
        {
            return bothPhases;
        }
        else if (std::strcmp(initialPhasePresence_.c_str(), "nPhaseOnly") == 0)
        {
            return nPhaseOnly;
        }
        else
        {
            DUNE_THROW(Dune::NotImplemented, "Initial phases presence is unkown " << initialPhasePresence_.c_str());
        }
    }

    /*!
     * \brief Called by the time manager after the time integration to
     *        do some post processing on the solution.
     */
    void postTimeStep()
    {
        // Calculate masses
        PrimaryVariables storage;

        globalStorage(storage);
        const Scalar time = this->timeManager().time() +  this->timeManager().timeStepSize();

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0)
        {
            if (this->timeManager().timeStepIndex() % freqFluxOutput_ == 0
                || this->timeManager().episodeWillBeFinished()
                || this->timeManager().willBeFinished())
            {

                GlobalPosition globalPos(0.0);
                PrimaryVariables storageChange(0.);
                storageChange = storageLastTimestep_ - storage;

                assert(time - lastMassOutputTime_ != 0);
                storageChange /= (time - lastMassOutputTime_);

                Scalar evaprate = storageChange[contiWEqIdx] / (bBoxMax_[0]-bBoxMin_[0]) * 86400.0;
#if DUMUX_MULTIDOMAIN_DIM > 2 // TODO
                evaprate /= (bBoxMax_[2]-bBoxMin_[2]);
#else
                evaprate /= this->extrusionFactorAtPos(globalPos);
#endif
                std::cout << "Time[s]: " << time
                          << " TotalMass[kg]: " << storage[contiTotalMassIdx]
                          << " WaterMass[kg]: " << storage[contiWEqIdx]
                          << " IntEnergy[J/m^3]: " << storage[energyEqIdx]
                          << " WaterMassChange[kg/s]: " << storageChange[contiWEqIdx]
                          << " EvaporationRate[mm/d]: " << evaprate
                          << std::endl;
                if (this->timeManager().time() != 0.)
                    outfile << time << ";"
                            << storageChange[contiTotalMassIdx] << ";"
                            << storageChange[contiWEqIdx] << ";"
                            << storageChange[energyEqIdx] << ";"
                            << storage[contiWEqIdx] << ";"
                            << initialWaterContent_ - storage[contiWEqIdx] << ";"
                            << evaprate << ";"
                            << moistureContent()
                            << std::endl;
                static double yMax = 15.0;
                static std::vector<double> x;
                static std::vector<double> y;

                x.push_back(time / 86400.0); // d
                y.push_back(evaprate);
                yMax = std::min(15.0,std::max(yMax, evaprate));

                gnuplot_.resetPlot();
                gnuplot_.setDatafileSeparator(';');
                gnuplot_.setXRange(0, x[x.size()-1]);
                gnuplot_.setYRange(0, yMax);
                gnuplot_.setXlabel("time [d]");
                gnuplot_.setYlabel("evaporation rate [mm/d]");
                gnuplot_.addDataSetToPlot(x, y, name() + "_evaprate.csv");
                if (GET_PARAM_FROM_GROUP(TypeTag, bool, Output, PlotEvaporationRate))
                {
                    gnuplot_.plot(name() + "_evaprate");
                }

                if (GET_PARAM_FROM_GROUP(TypeTag, bool, TimeManager, AbortIfEvapRateIsConstant)
                    && (this->timeManager().timeStepIndex() % freqFluxOutput_ == 0
                        || this->timeManager().episodeWillBeFinished()
                        || this->timeManager().willBeFinished()))
                {
                    std::cout << "eNew/eOld[-]: " << storageChange[contiWEqIdx] / storageChangeLastTimestep_[contiWEqIdx]
                              << " de[mm/d]: " << evaprate - evaprateLastTimestep_
                              << " e[mm/d]: " << evaprate
                              << std::endl;

                    if (storageChangeLastTimestep_[contiWEqIdx] * 1.001 > storageChange[contiWEqIdx]
                        && storageChangeLastTimestep_[contiWEqIdx] * 0.999 < storageChange[contiWEqIdx]
                        && (time - lastMassOutputTime_) > GET_PARAM_FROM_GROUP(TypeTag, int, TimeManager, MaxTimeStepSize))
                    {
                        std::cout << "Evaporation rate is constant, simulation will be aborted." << std::endl;
                        std::cout << "Final steady state evaporation rate [mm/d]: " << evaprate
                                  << " velocity[m/s]: " << velocity_
                                  << " patchSize[m]: " << (bBoxMax_[0]-bBoxMin_[0])
                                  << " patchLocation[m]: " << (bBoxMax_[0]+bBoxMin_[0]) / 2.0
                                  << std::endl;
                        exit(0);
                    }
                }

                storageLastTimestep_ = storage;
                storageChangeLastTimestep_ = storageChange;
                evaprateLastTimestep_ = evaprate;
                lastMassOutputTime_ = time;
            }
        }
    }

    /*!
     * \brief Compute the integral over the domain of the storage
     *        terms of all conservation quantities.
     *
     * \param storage Stores the result
     */
    void globalStorage(PrimaryVariables &storage)
    {
        // for the staggered model the local jacobian has to be initialized
        localJacobian_.init(*this);
        storage = 0;

        for (const auto& element : elements(gridView_))
        {
            if(element.partitionType() == Dune::InteriorEntity
               && gridView_.indexSet().contains(darcySubDomainIdx, element))
            {
                localResidual().evalStorage(element);
                storage += localResidual().storageTerm()[0];
            }
        }

        if (gridView_.comm().size() > 1)
            storage = gridView_.comm().sum(storage);
    }

    /*!
     * \brief Compute the integral moisture content
     */
    Scalar moistureContent()
    {
        Scalar vaporContent = 0.;
        Scalar liquidWaterContent = 0.;
        Scalar solidContent = 0.;

        for (const auto& element : elements(gridView_))
        {
            if(element.partitionType() == Dune::InteriorEntity
               && gridView_.indexSet().contains(darcySubDomainIdx, element))
            {
                FVElementGeometry fvGeometry;
                fvGeometry.update(gridView_, element);

                ElementVolumeVariables elemVolVars;
                elemVolVars.update(*this,
                                   element,
                                   fvGeometry,
                                   false /* oldSol? */);

                for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
                {
                    Scalar porosity = this->spatialParams().porosityAtPos(element.geometry().center());
                    liquidWaterContent += porosity
                                          * elemVolVars[scvIdx].saturation(wPhaseIdx)
                                          * elemVolVars[scvIdx].density(wPhaseIdx)
                                          * fvGeometry.subContVol[scvIdx].volume;

                    vaporContent += porosity
                                    * elemVolVars[scvIdx].saturation(nPhaseIdx)
                                    * elemVolVars[scvIdx].density(nPhaseIdx)
                                    * elemVolVars[scvIdx].moleFraction(nPhaseIdx, wCompIdx)
                                    * fvGeometry.subContVol[scvIdx].volume;

                    solidContent += (1.0 - porosity)
                                    * this->spatialParams().solidDensityAtPos(element.geometry().center())
                                    * fvGeometry.subContVol[scvIdx].volume;
                }
            }
        }
        return (liquidWaterContent + vaporContent)
               / (liquidWaterContent + vaporContent + solidContent);
    }

    /*!
     * \brief The depth of the problem in third dimension
     */
    Scalar extrusionFactorAtPos(const GlobalPosition &globalPos) const
    { return extrusionFactor_; }

    /*!
     * \brief Reference to the local residual object
     */
    LocalResidual &localResidual()
    { return localJacobian_.localResidual(); }

    //! \brief Returns the acceleration due to gravity
    const GlobalPosition &gravity() const
    { return gravity_; }
    // \}

    /*!
     * \brief Set the coupling manager
     *
     * \param couplingManager The coupling manager for the global problem
     */
    void setCouplingManager(std::shared_ptr<CouplingManager> couplingManager)
    {
        couplingManager_ = couplingManager;
    }

    /*!
     * \brief Get the coupling manager
     */
    CouplingManager& couplingManager() const
    { return *couplingManager_; }

    /*!
     * \brief Check if on coupling interface
     *
     * \param globalPos The global position
     *
     * Returns true if globalPos is on coupling interface
     * (here: upper boundary of Darcy domain)
     */
    bool onCouplingInterface(const GlobalPosition &globalPos) const
    {return onUpperBoundary_(globalPos); }

private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < bBoxMin_[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > bBoxMax_[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < bBoxMin_[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > bBoxMax_[1] - eps_; }

    static constexpr Scalar eps_ = 1e-8;

    GlobalPosition bBoxMin_;
    GlobalPosition bBoxMax_;
    Scalar extrusionFactor_;
    Scalar interfacePos_;

    GlobalPosition gravity_;
    Scalar velocity_;
    Scalar pressure_;
    Scalar switch_;
    std::string initialPhasePresence_;
    Scalar temperature_;
    bool solDependentEnergySource_;
    bool solDependentEnergyWalls_;
    bool useTemperatureDataFilePM_;
    std::string temperatureDataFilePM_;
    std::vector<double> temperatureDataPM_[2];

    int freqFluxOutput_;
    PrimaryVariables storageLastTimestep_;
    PrimaryVariables storageChangeLastTimestep_;
    Scalar evaprateLastTimestep_;
    Scalar initialWaterContent_;
    Scalar lastMassOutputTime_;
    std::ofstream outfile;
    Dumux::GnuplotInterface<double> gnuplot_;

    LocalJacobian localJacobian_;
    GridView gridView_;

    std::shared_ptr<CouplingManager> couplingManager_;
};

} //end namespace

#endif // DUMUX_WINDTUNNEL_DARCY_SUBPROBLEM_HH
