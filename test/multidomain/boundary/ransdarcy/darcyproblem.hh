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
 * \brief The porous medium sub problem
 */
#ifndef DUMUX_DARCY2P2C_SUBPROBLEM_HH
#define DUMUX_DARCY2P2C_SUBPROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/io/gnuplotinterface.hh>
#include <dumux/porousmediumflow/2p2c/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/material/fluidsystems/h2oair.hh>

#include "spatialparams.hh"

namespace Dumux
{
template <class TypeTag>
class DarcySubProblem;

namespace Properties
{
NEW_TYPE_TAG(DarcyTwoPTwoCTypeTag, INHERITS_FROM(CCTpfaModel, TwoPTwoCNI));

// Set the problem property
SET_TYPE_PROP(DarcyTwoPTwoCTypeTag, Problem, Dumux::DarcySubProblem<TypeTag>);

// the fluid system
SET_TYPE_PROP(DarcyTwoPTwoCTypeTag, FluidSystem, FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! Set the default formulation to pw-Sn: This can be over written in the problem.
SET_PROP(DarcyTwoPTwoCTypeTag, Formulation)
{ static constexpr auto value = TwoPFormulation::p1s0; };

// The gas component balance (air) is replaced by the total mass balance
SET_INT_PROP(DarcyTwoPTwoCTypeTag, ReplaceCompEqIdx, 3);

// Set the grid type
SET_TYPE_PROP(DarcyTwoPTwoCTypeTag, Grid, Dune::YaspGrid<2, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);

SET_BOOL_PROP(DarcyTwoPTwoCTypeTag, UseMoles, true);

SET_TYPE_PROP(DarcyTwoPTwoCTypeTag, SpatialParams, TwoPSpatialParams<TypeTag>);
}

/*!
 * \brief The porous medium sub problem
 */
template <class TypeTag>
class DarcySubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;

    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    // copy some indices for convenience
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    enum {
        // primary variable indices
        conti0EqIdx = Indices::conti0EqIdx,
        contiWEqIdx = Indices::conti0EqIdx + FluidSystem::H2OIdx,
        contiNEqIdx = Indices::conti0EqIdx + FluidSystem::AirIdx,
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);
    using TimeLoopPtr = std::shared_ptr<TimeLoop<Scalar>>;

    using DiffusionCoefficientAveragingType = typename StokesDarcyCouplingOptions::DiffusionCoefficientAveragingType;

public:
    DarcySubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                   std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(fvGridGeometry, "Darcy"), eps_(1e-7), couplingManager_(couplingManager)
    {
        pressure_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.Pressure");
        initialSw_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.Saturation");
        temperature_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.Temperature");
        initialPhasePresence_ = getParamFromGroup<int>(this->paramGroup(), "Problem.InitPhasePresence");

        diffCoeffAvgType_ = StokesDarcyCouplingOptions::stringToEnum(DiffusionCoefficientAveragingType{},
                                                                     getParamFromGroup<std::string>(this->paramGroup(), "Problem.InterfaceDiffusionCoefficientAvg"));

        // initialize output file
        plotFluxes_ = getParamFromGroup<bool>(this->paramGroup(), "Problem.PlotFluxes", false);
        plotStorage_ = getParamFromGroup<bool>(this->paramGroup(), "Problem.PlotStorage", false);
        storageFileName_ = "storage_" + getParam<std::string>("Problem.Name") + "_" + this->name() + ".csv";
        storageFile_.open(storageFileName_);
        storageFile_ << "#Time[s]" << ";"
                     << "WaterMass[kg]" << ";"
                     << "WaterMassLoss[kg]" << ";"
                     << "EvaporationRate[mm/d]"
                     << std::endl;
    }

    /*!
     * \name Simulation steering
     */
    // \{

    /*!
     * \brief Initialize the problem.
     */
    template<class SolutionVector, class GridVariables>
    void init(const SolutionVector& curSol,
              const GridVariables& gridVariables)
    {
        initialWaterContent_ = evaluateWaterMassStorageTerm(curSol, gridVariables);
        lastWaterMass_ = initialWaterContent_;
    }

    template<class SolutionVector, class GridVariables>
    void postTimeStep(const SolutionVector& curSol,
                      const GridVariables& gridVariables)

    {
        evaluateWaterMassStorageTerm(curSol, gridVariables);
        evaluateInterfaceFluxes(curSol, gridVariables);

        gnuplotStorage_.resetPlot();
        gnuplotStorage_.setDatafileSeparator(';');
        gnuplotStorage_.setXlabel("time [d]");
        gnuplotStorage_.setXRange(0.0, getParam<Scalar>("TimeLoop.TEnd") / 86400);
        gnuplotStorage_.setYlabel("evaporation rate [mm/d]");
        gnuplotStorage_.setOption("set yrange [0.0:15.0]");
        gnuplotStorage_.setOption("set y2label 'cumulative mass loss [kg]'");
        gnuplotStorage_.setOption("set y2range [0.0:15.0]");
        gnuplotStorage_.setOption("set ytics nomirror");
        gnuplotStorage_.setOption("set y2tics");

        gnuplotStorage_.addFileToPlot(storageFileName_, "using ($1/86400):4 with lines title 'evaporation rate'");
        gnuplotStorage_.addFileToPlot(storageFileName_, "using ($1/86400):3 axes x1y2 with lines title 'cumulative mass loss'");
        if (plotStorage_)
            gnuplotStorage_.plot("temp");
    }

    template<class SolutionVector, class GridVariables>
    Scalar evaluateWaterMassStorageTerm(const SolutionVector& curSol,
                                        const GridVariables& gridVariables)

    {
        // compute the mass in the entire domain
        Scalar waterMass = 0.0;

        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVariables.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, curSol);

            for (auto&& scv : scvs(fvGeometry))
            {
                const auto& volVars = elemVolVars[scv];
                for(int phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
                {
                    // insert calculation of the water mass here
                    waterMass += volVars.massFraction(phaseIdx, FluidSystem::H2OIdx) * volVars.density(phaseIdx)
                                 * volVars.saturation(phaseIdx) * volVars.porosity()
                                 * scv.volume() * volVars.extrusionFactor();
                }
            }
        }

        Scalar cumMassLoss = initialWaterContent_ - waterMass;
        Scalar evaporationRate = (lastWaterMass_ - waterMass) * 86400
                                 / (this->fvGridGeometry().bBoxMax()[0] - this->fvGridGeometry().bBoxMin()[0])
                                 / timeLoop_->timeStepSize();
        lastWaterMass_ = waterMass;

        storageFile_ << timeLoop_->time() << ";"
                     << waterMass << ";"
                     << cumMassLoss << ";"
                     << evaporationRate
                     << std::endl;

        std::cout << "Mass of water is: " << waterMass << std::endl;
        std::cout << "Evaporation rate is: " << evaporationRate << std::endl;

        return waterMass;
    }

    template<class SolutionVector, class GridVariables>
    void evaluateInterfaceFluxes(const SolutionVector& curSol,
                                 const GridVariables& gridVariables)

    {
        using std::max;
        using std::min;
        std::vector<Scalar> x;
        std::vector<Scalar> y;
        static Scalar maxFlux = -9e9;
        static Scalar minFlux = 9e9;

        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVariables.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, curSol);

            for (auto&& scvf : scvfs(fvGeometry))
            {
                if (!couplingManager().isCoupledEntity(CouplingManager::darcyIdx, scvf))
                    continue;

                // NOTE: binding the coupling context is necessary
                couplingManager_->bindCouplingContext(CouplingManager::darcyIdx, element);
                const auto massFlux = couplingManager().couplingData().massCouplingCondition(fvGeometry, elemVolVars, scvf, diffCoeffAvgType_);
                NumEqVector flux(0.0);
                for(int i = 0; i< massFlux.size(); ++i)
                    flux[i] = massFlux[i];

                x.push_back(scvf.center()[0]);
                y.push_back(flux[contiWEqIdx]);
                maxFlux = max(maxFlux,flux[contiWEqIdx]);
                minFlux = min(minFlux,flux[contiWEqIdx]);
            }
        }

        gnuplotInterfaceFluxes_.resetPlot();
        gnuplotInterfaceFluxes_.setXlabel("x-position [m]");
        gnuplotInterfaceFluxes_.setXRange(this->fvGridGeometry().bBoxMin()[0], this->fvGridGeometry().bBoxMax()[0]);
        gnuplotInterfaceFluxes_.setYlabel("flux [kg/(m^2 s)]");
        gnuplotInterfaceFluxes_.setYRange(minFlux, maxFlux);
        gnuplotInterfaceFluxes_.setOption("set label 'time: " + std::to_string(timeLoop_->time()/86400.) + "d' at graph 0.8,0.8 ");
        std::string fluxFileName = "flux_" + std::to_string(timeLoop_->timeStepIndex()) +
                                   "_" + getParam<std::string>("Problem.Name") + "_" + this->name() + ".csv";
        gnuplotInterfaceFluxes_.addDataSetToPlot(x, y, fluxFileName, "with lines title 'water mass flux'");
        if (plotFluxes_)
            gnuplotInterfaceFluxes_.plot("flux_" + std::to_string(timeLoop_->timeStepIndex()));
    }

    /*!
     * \brief Return the temperature within the domain in [K].
     */
    Scalar temperature() const
    { return temperature_; }
    // \}

     /*!
     * \name Boundary conditions
     */
    // \{
    /*!
      * \brief Specifies which kind of boundary condition should be
      *        used for which equation on a given boundary control volume.
      *
      * \param element The element
      * \param scvf The boundary sub control volume face
      */
    BoundaryTypes boundaryTypes(const Element &element, const SubControlVolumeFace &scvf) const
    {
        BoundaryTypes values;
        values.setAllNeumann();

        if (couplingManager().isCoupledEntity(CouplingManager::darcyIdx, scvf))
            values.setAllCouplingNeumann();

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Dirichlet control volume.
     *
     * \param element The element for which the Dirichlet boundary condition is set
     * \param scvf The boundary subcontrolvolumeface
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichlet(const Element &element, const SubControlVolumeFace &scvf) const
    {
        PrimaryVariables values(0.0);
        values = initialAtPos(scvf.center());

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann
     *        control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeomentry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param scvf The boundary sub control volume face
     *
     */
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        if (couplingManager().isCoupledEntity(CouplingManager::darcyIdx, scvf))
        {
            const auto massFlux = couplingManager().couplingData().massCouplingCondition(fvGeometry, elemVolVars, scvf, diffCoeffAvgType_);

            for(int i = 0; i< massFlux.size(); ++i)
                values[i] = massFlux[i];

            values[Indices::energyEqIdx] = couplingManager().couplingData().energyCouplingCondition(fvGeometry, elemVolVars, scvf, diffCoeffAvgType_);
        }

        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{
    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * \param element The element for which the source term is set
     * \param fvGeomentry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param scv The subcontrolvolume
     *
     * For this method, the \a values variable stores the rate mass
     * of a component is generated or annihilate per volume
     * unit. Positive values mean that mass is created, negative ones
     * mean that it vanishes.
     */
    NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume &scv) const
    { return NumEqVector(0.0); }

    // \}

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values.setState(initialPhasePresence_);

        values[pressureIdx] = pressure_ + 1. * this->gravity()[1] * (globalPos[1] - this->fvGridGeometry().bBoxMax()[1]);
        values[switchIdx] = initialSw_;
        values[Indices::temperatureIdx] = temperature_;

        return values;
    }

    // \}

    /*!
     * \brief Set the coupling manager
     */
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    /*!
     * \brief Get the coupling manager
     */
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    void setTimeLoop(TimeLoopPtr timeLoop)
    { timeLoop_ = timeLoop; }

private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->fvGridGeometry().bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->fvGridGeometry().bBoxMax()[1] - eps_; }

    Scalar pressure_;
    Scalar initialSw_;
    Scalar temperature_;
    int initialPhasePresence_;

    TimeLoopPtr timeLoop_;

    Scalar eps_;

    std::shared_ptr<CouplingManager> couplingManager_;
    DiffusionCoefficientAveragingType diffCoeffAvgType_;

    std::string storageFileName_;
    std::ofstream storageFile_;
    bool plotFluxes_;
    bool plotStorage_;
    Scalar initialWaterContent_ = 0.0;
    Scalar lastWaterMass_ = 0.0;
    Dumux::GnuplotInterface<Scalar> gnuplotInterfaceFluxes_;
    Dumux::GnuplotInterface<Scalar> gnuplotStorage_;
};
} //end namespace

#endif //DUMUX_DARCY2P2C_SUBPROBLEM_HH
