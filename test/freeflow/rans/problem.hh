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
 * \ingroup RANSTests
 * \brief Pipe flow test for the staggered grid RANS model
 *
 * This test simulates pipe flow experiments performed by John Laufer in 1954
 * \cite Laufer1954a.
 */
#ifndef DUMUX_PIPE_LAUFER_PROBLEM_HH
#define DUMUX_PIPE_LAUFER_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>
#include <dune/common/hybridutilities.hh>

#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/freeflow/turbulenceproperties.hh>
#include <dumux/freeflow/rans/problem.hh>
#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/components/air.hh>

#include <dumux/freeflow/rans/zeroeq/model.hh>
#include <dumux/freeflow/turbulencemodel.hh>
#include <dumux/freeflow/rans/oneeq/problem.hh>
#include <dumux/freeflow/rans/oneeq/model.hh>
#include <dumux/freeflow/rans/twoeq/lowrekepsilon/problem.hh>
#include <dumux/freeflow/rans/twoeq/lowrekepsilon/model.hh>
#include <dumux/freeflow/rans/twoeq/komega/problem.hh>
#include <dumux/freeflow/rans/twoeq/komega/model.hh>
#include <dumux/freeflow/rans/twoeq/kepsilon/problem.hh>
#include <dumux/freeflow/rans/twoeq/kepsilon/model.hh>

namespace Dumux {

template <class TypeTag>
class PipeLauferProblem;

namespace Properties {

// Create new type tags
namespace TTag {
// Base Typetag
struct RANSModel { using InheritsFrom = std::tuple<StaggeredFreeFlowModel>; };
// Isothermal Typetags
struct PipeLauferZeroEq { using InheritsFrom = std::tuple<RANSModel, ZeroEq>; };
struct PipeLauferOneEq { using InheritsFrom = std::tuple<RANSModel, OneEq>; };
struct PipeLauferKOmega { using InheritsFrom = std::tuple<RANSModel, KOmega>; };
struct PipeLauferLowReKEpsilon { using InheritsFrom = std::tuple<RANSModel, LowReKEpsilon>; };
struct PipeLauferKEpsilon { using InheritsFrom = std::tuple<RANSModel, KEpsilon>; };
// Non-Isothermal Typetags
struct PipeLauferNIZeroEq { using InheritsFrom = std::tuple<RANSModel, ZeroEqNI>; };
struct PipeLauferNIOneEq { using InheritsFrom = std::tuple<RANSModel, OneEqNI>; };
struct PipeLauferNIKOmega { using InheritsFrom = std::tuple<RANSModel, KOmegaNI>; };
struct PipeLauferNILowReKEpsilon { using InheritsFrom = std::tuple<RANSModel, LowReKEpsilonNI>; };
struct PipeLauferNIKEpsilon { using InheritsFrom = std::tuple<RANSModel, KEpsilonNI>; };
} // end namespace TTag


// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::RANSModel>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePGas<Scalar, Components::Air<Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::RANSModel>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::RANSModel>
{ using type = Dumux::PipeLauferProblem<TypeTag>; };

template<class TypeTag>
struct EnableFVGridGeometryCache<TypeTag, TTag::RANSModel> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::RANSModel> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::RANSModel> { static constexpr bool value = true; };
} // end namespace Properties

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the one-phase (Navier-) Stokes problem in a channel.
 *
 * This test simulates is based on pipe flow experiments by
 * John Laufers experiments in 1954 \cite Laufer1954a.
 */
template <class TypeTag>
class PipeLauferProblem : public RANSProblem<TypeTag>
{
    using ParentType = RANSProblem<TypeTag>;

    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Element = typename FVGridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using TimeLoopPtr = std::shared_ptr<CheckPointTimeLoop<Scalar>>;

    static constexpr auto dimWorld = FVGridGeometry::GridView::dimensionworld;

public:
    PipeLauferProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry), eps_(1e-6)
    {
        inletVelocity_ = getParam<Scalar>("Problem.InletVelocity");
        inletTemperature_ = getParam<Scalar>("Problem.InletTemperature", 283.15);
        wallTemperature_ = getParam<Scalar>("Problem.WallTemperature", 303.15);
        sandGrainRoughness_ = getParam<Scalar>("Problem.SandGrainRoughness", 0.0);

        FluidSystem::init();
        Dumux::TurbulenceProperties<Scalar, dimWorld, true> turbulenceProperties;
        FluidState fluidState;
        fluidState.setPressure(0, 1e5);
        fluidState.setTemperature(temperature());
        Scalar density = FluidSystem::density(fluidState, 0);
        Scalar kinematicViscosity = FluidSystem::viscosity(fluidState, 0) / density;
        Scalar diameter = this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1];

        // ideally the viscosityTilde parameter as inflow for the Spalart-Allmaras model should be zero
        viscosityTilde_ = 1e-3 * turbulenceProperties.viscosityTilde(inletVelocity_, diameter, kinematicViscosity);
        turbulentKineticEnergy_ = turbulenceProperties.turbulentKineticEnergy(inletVelocity_, diameter, kinematicViscosity);

        if (ModelTraits::turbulenceModel() == TurbulenceModel::komega)
            dissipation_ = turbulenceProperties.dissipationRate(inletVelocity_, diameter, kinematicViscosity);
        else
            dissipation_ = turbulenceProperties.dissipation(inletVelocity_, diameter, kinematicViscosity);

        if (ModelTraits::turbulenceModel() == TurbulenceModel::oneeq)
            initializationTime_ = getParam<Scalar>("TimeLoop.Initialization", 1.0);
        else
            initializationTime_ = getParam<Scalar>("TimeLoop.Initialization", -1.0);
        std::cout << std::endl;
    }

   /*!
     * \name Problem parameters
     */
    // \{

    bool isOnWallAtPos(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->fvGridGeometry().bBoxMin()[1] + eps_
               || globalPos[1] > this->fvGridGeometry().bBoxMax()[1] - eps_;
    }

    Scalar sandGrainRoughnessAtPos(const GlobalPosition &globalPos) const
    {
        return sandGrainRoughness_;
    }

    bool shouldWriteRestartFile() const
    {
        return false;
    }

   /*!
     * \brief Returns the temperature [K] within the domain for the isothermal model.
     */
    Scalar temperature() const
    { return inletTemperature_; }

   /*!
     * \brief Returns the sources within the domain.
     *
     * \param globalPos The global position
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        return NumEqVector(0.0);
    }
    // \}

   /*!
     * \name Boundary conditions
     */
    // \{

   /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        // common boundary types for all turbulence models
        if(isOutlet_(globalPos))
            values.setDirichlet(Indices::pressureIdx);
        else
        {
            // walls and inflow
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
        }

#if NONISOTHERMAL
        if(isOutlet_(globalPos))
            values.setOutflow(Indices::energyEqIdx);
        else
            values.setDirichlet(Indices::temperatureIdx);
#endif

        // turbulence model-specific boundary types
        static constexpr auto numEq = numTurbulenceEq(ModelTraits::turbulenceModel());
        setBcTypes_(values, globalPos, Dune::index_constant<numEq>{});

        return values;
    }

    /*!
     * \brief Returns whether a fixed Dirichlet value shall be used at a given cell.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param scv The sub control volume
     */
    template<class Element, class FVElementGeometry, class SubControlVolume>
    bool isDirichletCell(const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const SubControlVolume& scv,
                         int pvIdx) const
    {
        using IsKOmegaKEpsilon = std::integral_constant<bool, (ModelTraits::turbulenceModel() == TurbulenceModel::komega
                                                            || ModelTraits::turbulenceModel() == TurbulenceModel::kepsilon)>;
        return isDirichletCell_(element, fvGeometry, scv, pvIdx, IsKOmegaKEpsilon{});
    }

     /*!
      * \brief Evaluate the boundary conditions for a dirichlet values at the boundary.
      *
      * \param element The finite element
      * \param scvf the sub control volume face
      * \note used for cell-centered discretization schemes
      */
    PrimaryVariables dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    {
        const auto globalPos = scvf.ipGlobal();
        PrimaryVariables values(initialAtPos(globalPos));

#if NONISOTHERMAL
        values[Indices::temperatureIdx] = isOnWallAtPos(globalPos) ? wallTemperature_ : inletTemperature_;
#endif

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for fixed values at cell centers
     *
     * \param element The finite element
     * \param scv the sub control volume
     * \note used for cell-centered discretization schemes
     */
    template<bool enable = (ModelTraits::turbulenceModel() == TurbulenceModel::komega
                         || ModelTraits::turbulenceModel() == TurbulenceModel::kepsilon),
                         std::enable_if_t<!enable, int> = 0>
    PrimaryVariables dirichlet(const Element& element, const SubControlVolume& scv) const
    {
        const auto globalPos = scv.center();
        PrimaryVariables values(initialAtPos(globalPos));
        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for fixed values at cell centers
     *
     * \param element The finite element
     * \param scv the sub control volume
     * \note used for cell-centered discretization schemes
     */
    template<bool enable = (ModelTraits::turbulenceModel() == TurbulenceModel::komega
                         || ModelTraits::turbulenceModel() == TurbulenceModel::kepsilon),
                         std::enable_if_t<enable, int> = 0>
    PrimaryVariables dirichlet(const Element& element, const SubControlVolume& scv) const
    {
        using SetDirichletCellForBothTurbEq = std::integral_constant<bool, (ModelTraits::turbulenceModel() == TurbulenceModel::kepsilon)>;

        return dirichletTurbulentTwoEq_(element, scv, SetDirichletCellForBothTurbEq{});
    }

   /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        // common initial conditions for all turbulence models
        PrimaryVariables values(0.0);
        values[Indices::pressureIdx] = 1.0e+5;
        values[Indices::velocityXIdx] = time() > initializationTime_
                                        ? inletVelocity_
                                        : time() / initializationTime_ * inletVelocity_;
        if (isOnWallAtPos(globalPos))
            values[Indices::velocityXIdx] = 0.0;

#if NONISOTHERMAL
        values[Indices::temperatureIdx] = isOnWallAtPos(globalPos) ? wallTemperature_ : inletTemperature_;
#endif

        // turbulence model-specific initial conditions
        static constexpr auto numEq = numTurbulenceEq(ModelTraits::turbulenceModel());
        setInitialAtPos_(values, globalPos, Dune::index_constant<numEq>{});

        return values;
    }

    // \}

    void setTimeLoop(TimeLoopPtr timeLoop)
    {
        timeLoop_ = timeLoop;
    }

    Scalar time() const
    {
        return timeLoop_->time();
    }

private:
    bool isInlet_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + eps_;
    }

    bool isOutlet_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_;
    }

    //! Initial conditions for the zero-eq turbulence model (none)
    void setInitialAtPos_(PrimaryVariables& values, const GlobalPosition &globalPos, Dune::index_constant<0>) const {}

    //! Initial conditions for the one-eq turbulence model
    void setInitialAtPos_(PrimaryVariables& values, const GlobalPosition &globalPos, Dune::index_constant<1>) const
    {
        values[Indices::viscosityTildeIdx] = viscosityTilde_;
        if (isOnWallAtPos(globalPos))
            values[Indices::viscosityTildeIdx] = 0.0;
    }

    //! Initial conditions for the komega, kepsilon and lowrekepsilon turbulence models
    void setInitialAtPos_(PrimaryVariables& values, const GlobalPosition &globalPos, Dune::index_constant<2>) const
    {
        values[Indices::turbulentKineticEnergyIdx] = turbulentKineticEnergy_;
        values[Indices::dissipationIdx] = dissipation_;
        if (isOnWallAtPos(globalPos))
        {
            values[Indices::turbulentKineticEnergyIdx] = 0.0;
            values[Indices::dissipationIdx] = 0.0;
        }
    }

    //! Boundary condition types for the zero-eq turbulence model (none)
    void setBcTypes_(BoundaryTypes& values, const GlobalPosition& pos, Dune::index_constant<0>) const {}

    //! Boundary condition types for the one-eq turbulence model
    void setBcTypes_(BoundaryTypes& values, const GlobalPosition& pos, Dune::index_constant<1>) const
    {
        if(isOutlet_(pos))
            values.setOutflow(Indices::viscosityTildeIdx);
        else // walls and inflow
            values.setDirichlet(Indices::viscosityTildeIdx);
    }

    //! Boundary condition types for the komega, kepsilon and lowrekepsilon turbulence models
    void setBcTypes_(BoundaryTypes& values,const GlobalPosition& pos, Dune::index_constant<2>) const
    {
        if(isOutlet_(pos))
        {
            values.setOutflow(Indices::turbulentKineticEnergyEqIdx);
            values.setOutflow(Indices::dissipationEqIdx);
        }
        else
        {
            // walls and inflow
            values.setDirichlet(Indices::turbulentKineticEnergyIdx);
            values.setDirichlet(Indices::dissipationIdx);
        }
    }

    //! Forward to ParentType
    template<class Element, class FVElementGeometry, class SubControlVolume>
    bool isDirichletCell_(const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const SubControlVolume& scv,
                          int pvIdx,
                          std::false_type) const
    {
        return ParentType::isDirichletCell(element, fvGeometry, scv, pvIdx);
    }

    //! Specialization for the KOmega and KEpsilon Models
    template<class Element, class FVElementGeometry, class SubControlVolume>
    bool isDirichletCell_(const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const SubControlVolume& scv,
                          int pvIdx,
                          std::true_type) const
    {
        using SetDirichletCellForBothTurbEq = std::integral_constant<bool, (ModelTraits::turbulenceModel() == TurbulenceModel::kepsilon)>;
        return isDirichletCellTurbulentTwoEq_(element, fvGeometry, scv, pvIdx, SetDirichletCellForBothTurbEq{});
    }

    //! Specialization for the KEpsilon Model
    template<class Element>
    bool isDirichletCellTurbulentTwoEq_(const Element& element,
                                        const FVElementGeometry& fvGeometry,
                                        const SubControlVolume& scv,
                                        int pvIdx,
                                        std::true_type) const
    {
        const auto eIdx = this->fvGridGeometry().elementMapper().index(element);

        // set a fixed turbulent kinetic energy and dissipation near the wall
        if (this->inNearWallRegion(eIdx))
            return pvIdx == Indices::turbulentKineticEnergyEqIdx || pvIdx == Indices::dissipationEqIdx;

        // set a fixed dissipation at  the matching point
        if (this->isMatchingPoint(eIdx))
            return pvIdx == Indices::dissipationEqIdx;// set a fixed dissipation (omega) for all cells at the wall

        return false;
    }

    //! Specialization for the KOmega Model
    template<class Element>
    bool isDirichletCellTurbulentTwoEq_(const Element& element,
                                        const FVElementGeometry& fvGeometry,
                                        const SubControlVolume& scv,
                                        int pvIdx,
                                        std::false_type) const
    {
        // set a fixed dissipation (omega) for all cells at the wall
        for (const auto& scvf : scvfs(fvGeometry))
            if (isOnWallAtPos(scvf.center()) && pvIdx == Indices::dissipationIdx)
                return true;

        return false;

    }

    //! Specialization for the kepsilon
    template<class Element, class SubControlVolume>
    PrimaryVariables dirichletTurbulentTwoEq_(const Element& element,
                                              const SubControlVolume& scv,
                                              std::true_type) const
    {
        const auto globalPos = scv.center();
        PrimaryVariables values(initialAtPos(globalPos));
        unsigned int  elementIdx = this->fvGridGeometry().elementMapper().index(element);

        // fixed value for the turbulent kinetic energy
        values[Indices::turbulentKineticEnergyEqIdx] = this->turbulentKineticEnergyWallFunction(elementIdx);

        // fixed value for the dissipation
        values[Indices::dissipationEqIdx] = this->dissipationWallFunction(elementIdx);

        return values;
    }

    //! Specialization for the KOmega
    template<class Element, class SubControlVolume>
    PrimaryVariables dirichletTurbulentTwoEq_(const Element& element,
                                              const SubControlVolume& scv,
                                              std::false_type) const
    {
        const auto globalPos = scv.center();
        PrimaryVariables values(initialAtPos(globalPos));
        unsigned int  elementIdx = this->fvGridGeometry().elementMapper().index(element);

        const auto wallDistance = ParentType::wallDistance_[elementIdx];
        using std::pow;
        values[Indices::dissipationEqIdx] = 6.0 * ParentType::kinematicViscosity_[elementIdx]
                                                / (ParentType::betaOmega() * pow(wallDistance, 2));
        return values;
    }

    Scalar eps_;
    Scalar inletVelocity_;
    Scalar inletTemperature_;
    Scalar wallTemperature_;
    Scalar sandGrainRoughness_;
    Scalar initializationTime_;
    Scalar viscosityTilde_;
    Scalar turbulentKineticEnergy_;
    Scalar dissipation_;
    TimeLoopPtr timeLoop_;
};
} // end namespace Dumux

#endif
