// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup BoundaryTests
 * \brief A Darcy Subproblem (cell-centered finite volume method).
 */

#ifndef DUMUX_DARCY2P2C_SUBPROBLEM_HH
#define DUMUX_DARCY2P2C_SUBPROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/multidomain/boundary/stokesdarcy/couplingdata.hh>

#include <dumux/porousmediumflow/2p2c/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivityconstanttortuosity.hh>

#include "spatialparams.hh"

namespace Dumux {
template <class TypeTag>
class DarcySubProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct DarcyTwoPTwoCNI { using InheritsFrom = std::tuple<TwoPTwoCNI, CCTpfaModel>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DarcyTwoPTwoCNI> { using type = Dumux::DarcySubProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DarcyTwoPTwoCNI> { using type = FluidSystems::H2OAir<GetPropType<TypeTag, Properties::Scalar>>; };

template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::DarcyTwoPTwoCNI> { static constexpr int value = 10; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DarcyTwoPTwoCNI>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::DarcyTwoPTwoCNI>
{
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = TwoPTwoCSpatialParams<FVGridGeometry, Scalar>;
};

template<class TypeTag>
struct UseMoles<TypeTag, TTag::DarcyTwoPTwoCNI> { static constexpr bool value = true; };

//! Set the default formulation to pw-Sn: This can be over written in the problem.
template<class TypeTag>
struct Formulation<TypeTag, TTag::DarcyTwoPTwoCNI>
{ static constexpr auto value = TwoPFormulation::p1s0; };

} // end namespace Properties

template <class TypeTag>
class DarcySubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;

    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    enum {
        // grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        conti0EqIdx = Indices::conti0EqIdx,
        contiWEqIdx = Indices::conti0EqIdx + FluidSystem::H2OIdx,
        contiNEqIdx = Indices::conti0EqIdx + FluidSystem::AirIdx,
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using TimeLoopPtr = std::shared_ptr<TimeLoop<Scalar>>;

    using DiffusionCoefficientAveragingType = typename StokesDarcyCouplingOptions::DiffusionCoefficientAveragingType;


public:
    DarcySubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                   std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(fvGridGeometry, "Darcy")
    , couplingManager_(couplingManager)
    , eps_(1e-7)

    {
        pressure_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.Pressure");
        initialSaturation_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InitialSaturation");
        temperature_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.Temperature");
        initialPhasePresence_ = getParamFromGroup<int>(this->paramGroup(), "Problem.InitPhasePresence");

        diffCoeffAvgType_ = StokesDarcyCouplingOptions::stringToEnum(DiffusionCoefficientAveragingType{},
                                                                     getParamFromGroup<std::string>(this->paramGroup(), "Problem.InterfaceDiffusionCoefficientAvg"));

        problemName_  =  getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
    }


    /*!
     * \name Mass Exchange Evaluations
     */
    // \{

    /*!
     * \brief Initialize the problem.
     */
    template<class SolutionVector, class GridVariables>
    void init(const SolutionVector& curSol,
              const GridVariables& gridVariables,
              const std::string ransName)
    {
        initialWaterContent_ = evaluateWaterMassStorageTerm(curSol, gridVariables);
        lastWaterMass_ = initialWaterContent_;

        storageFileName_ = "storage_" + ransName + ".csv";
        storageFile_.open(storageFileName_, std::ios::app);
        storageFile_ << "#Time[s]" << ";"
                     << "WaterMass[kg]" << ";"
                     << "WaterMassLoss[kg]" << ";"
                     << "EvaporationRate[mm/d]"
                     << std::endl;
        storageFile_.close();
    }

    /*!
     * \brief Evaluate the exchange processes after each timestep.
     */
    template<class SolutionVector, class GridVariables>
    void postTimeStep(const SolutionVector& curSol,
                      const GridVariables& gridVariables)

    {
        evaluateWaterMassStorageTerm(curSol, gridVariables);
    }

    /*!
     * \brief Calculate the total water mass storage exchange.
     */
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
        storageFile_.open(this->storageFileName_, std::ios::app);
        storageFile_ << timeLoop_->time() << ";"
                     << waterMass << ";"
                     << cumMassLoss << ";"
                     << evaporationRate
                     << "\n";
        storageFile_.close();

        std::cout << "Mass of water is: " << waterMass << "\n";
        std::cout << "Evaporation rate is: " << evaporationRate << "\n";
        std::cout << "\n";
        return waterMass;
    }

    // \}

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    {
        return problemName_;
    }


    /*!
     * \brief Returns the temperature within the domain in [K].
     *
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
     * \brief Evaluate the boundary conditions for a Neumann
     *        control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
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
            const auto massFlux = couplingManager().couplingData().massCouplingCondition(element, fvGeometry, elemVolVars, scvf, diffCoeffAvgType_);

            for(int i = 0; i< massFlux.size(); ++i)
                values[i] = massFlux[i];

            values[Indices::energyEqIdx] = couplingManager().couplingData().energyCouplingCondition(element, fvGeometry, elemVolVars, scvf, diffCoeffAvgType_);
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
     */
    template<class ElementVolumeVariables>
    NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume &scv) const
    { return NumEqVector(0.0); }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values.setState(Indices::bothPhases);

        FluidState fluidState;
        fluidState.setPressure(0, pressure_);
        fluidState.setTemperature(temperature());
        fluidState.setMoleFraction(0, 0, 1);
        Scalar density = FluidSystem::density(fluidState, 0);

        values[Indices::pressureIdx] = pressure_ - density * this->gravity()[dimWorld-1]* (this->fvGridGeometry().bBoxMax()[1] - globalPos[1]);
        values[switchIdx] = initialSaturation_;
        values[Indices::temperatureIdx] = temperature();

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
    Scalar temperature_;
    int initialPhasePresence_;
    Scalar initialSaturation_;

    TimeLoopPtr timeLoop_;
    std::shared_ptr<CouplingManager> couplingManager_;
    DiffusionCoefficientAveragingType diffCoeffAvgType_;
    Scalar eps_;

    std::string problemName_;
    std::string storageFileName_;
    std::ofstream storageFile_;

    bool plotFluxes_;
    bool plotStorage_;
    Scalar initialWaterContent_ = 0.0;
    Scalar lastWaterMass_ = 0.0;

};
} // end namespace Dumux

#endif //DUMUX_DARCY_SUBPROBLEM_HH
