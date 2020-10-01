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
 * \brief A simple Darcy test problem (cell-centered finite volume method).
 */

#ifndef DUMUX_DARCY2P2C_SUBPROBLEM_HH
#define DUMUX_DARCY2P2C_SUBPROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>

#include <dumux/common/boundarytypes.hh>

#include <dumux/material/fluidsystems/h2oair.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/2p2c/model.hh>

#include "spatialparams.hh"

namespace Dumux {
template <class TypeTag>
class DarcySubProblem;

namespace Properties {
// Create new type tags
namespace TTag {
#if !NONISOTHERMAL
struct DarcyTwoPTwoC { using InheritsFrom = std::tuple<TwoPTwoC, CCTpfaModel>; };
#else
struct DarcyTwoPTwoC { using InheritsFrom = std::tuple<TwoPTwoCNI, CCTpfaModel>; };
#endif
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DarcyTwoPTwoC> { using type = Dumux::DarcySubProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DarcyTwoPTwoC> { using type = FluidSystems::H2OAir<GetPropType<TypeTag, Properties::Scalar>>; };

//! Set the default formulation to pw-Sn: This can be over written in the problem.
template<class TypeTag>
struct Formulation<TypeTag, TTag::DarcyTwoPTwoC>
{ static constexpr auto value = TwoPFormulation::p1s0; };

//// The gas component balance (air) is replaced by the total mass balance
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::DarcyTwoPTwoC> { static constexpr int value = 3; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DarcyTwoPTwoC> { using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

template<class TypeTag>
struct UseMoles<TypeTag, TTag::DarcyTwoPTwoC> { static constexpr bool value = true; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::DarcyTwoPTwoC>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = TwoPTwoCSpatialParams<GridGeometry, Scalar>;
};

} // end namespace Properties

template <class TypeTag>
class DarcySubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    // copy some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
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

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    using DiffusionCoefficientAveragingType = typename StokesDarcyCouplingOptions::DiffusionCoefficientAveragingType;

public:
    DarcySubProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                   std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, "Darcy"), eps_(1e-7), couplingManager_(couplingManager)
    {
        pressure_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.Pressure");
        initialSw_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.Saturation");
        temperature_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.Temperature");
        initialPhasePresence_ = getParamFromGroup<int>(this->paramGroup(), "Problem.InitPhasePresence");

        diffCoeffAvgType_ = StokesDarcyCouplingOptions::stringToEnum(DiffusionCoefficientAveragingType{},
                                                                     getParamFromGroup<std::string>(this->paramGroup(), "Problem.InterfaceDiffusionCoefficientAvg"));
        problemName_  =  getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
    }

    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    {
        return problemName_;
    }

    template<class SolutionVector, class GridVariables>
    void printWaterMass(const SolutionVector& curSol,
                        const GridVariables& gridVariables,
                        const Scalar timeStepSize)

    {
        // compute the mass in the entire domain
        Scalar massWater = 0.0;

        // bulk elements
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVariables.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, curSol);

            for (auto&& scv : scvs(fvGeometry))
            {
                const auto& volVars = elemVolVars[scv];
                for(int phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
                {
                    massWater += volVars.massFraction(phaseIdx, FluidSystem::H2OIdx)*volVars.density(phaseIdx)
                    * scv.volume() * volVars.saturation(phaseIdx) * volVars.porosity() * volVars.extrusionFactor();
                }
            }
        }

        std::cout << std::setprecision(15) << "mass of water is: " << massWater << std::endl;
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the temperature within the domain in [K].
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
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param element The element for which the Dirichlet boundary condition is set
     * \param scvf The boundary sub control volume face
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
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFluxVarsCache Flux variables caches for all faces in stencil
     * \param scvf The boundary sub control volume face
     */
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        if (couplingManager().isCoupledEntity(CouplingManager::darcyIdx, scvf))
        {
#if !NONISOTHERMAL
            values = couplingManager().couplingData().massCouplingCondition(element, fvGeometry, elemVolVars, scvf, diffCoeffAvgType_);
#else
            const auto massFlux = couplingManager().couplingData().massCouplingCondition(element, fvGeometry, elemVolVars, scvf, diffCoeffAvgType_);

            for(int i = 0; i< massFlux.size(); ++i)
                values[i] = massFlux[i];

            values[Indices::energyEqIdx] = couplingManager().couplingData().energyCouplingCondition(element, fvGeometry, elemVolVars, scvf, diffCoeffAvgType_);
#endif
        }

        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{
    /*!
     * \brief Evaluates the source term for all phases within a given
     *        sub control volume.
     *
     * \param element The element for which the source term is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param scv The sub control volume
     *
     * For this method, the \a values variable stores the rate mass
     * of a component is generated or annihilated per volume
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
     * \brief Evaluates the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values.setState(initialPhasePresence_);

        values[pressureIdx] = pressure_ + 1000. * this->spatialParams().gravity(globalPos)[1] * (globalPos[1] - this->gridGeometry().bBoxMax()[1]);
        values[switchIdx] = initialSw_;

#if NONISOTHERMAL
        values[Indices::temperatureIdx] = temperature_;
#endif
        return values;
    }

    // \}

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_; }

    Scalar pressure_;
    Scalar initialSw_;
    Scalar temperature_;
    int initialPhasePresence_;
    std::string problemName_;
    Scalar eps_;

    std::shared_ptr<CouplingManager> couplingManager_;
    DiffusionCoefficientAveragingType diffCoeffAvgType_;
};
} // end namespace Dumux

#endif
