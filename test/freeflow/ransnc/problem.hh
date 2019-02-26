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
 * \ingroup RANSNCTests
 * \brief Flat plate test for the multi-component staggered grid Reynolds-averaged Navier-Stokes model.
 */
#ifndef DUMUX_RANS_NC_TEST_PROBLEM_HH
#define DUMUX_RANS_NC_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/material/fluidsystems/1padapter.hh>
#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/freeflow/turbulenceproperties.hh>

#include <dumux/freeflow/rans/zeroeq/problem.hh>
#include <dumux/freeflow/compositional/zeroeqncmodel.hh>

#include <dumux/freeflow/rans/oneeq/problem.hh>
#include <dumux/freeflow/compositional/oneeqncmodel.hh>

#include <dumux/freeflow/rans/twoeq/komega/problem.hh>
#include <dumux/freeflow/compositional/komegancmodel.hh>

#include <dumux/freeflow/rans/twoeq/lowrekepsilon/problem.hh>
#include <dumux/freeflow/compositional/lowrekepsilonncmodel.hh>

#include <dumux/freeflow/rans/twoeq/kepsilon/problem.hh>
#include <dumux/freeflow/compositional/kepsilonncmodel.hh>


namespace Dumux {

template <class TypeTag>
class FlatPlateNCTestProblem;

namespace Properties {

// Create new type tags
namespace TTag {
// Base Typetag
struct RANSNCModel { using InheritsFrom = std::tuple<StaggeredFreeFlowModel>; };
// Isothermal Typetags
struct FlatPlateNCZeroEq { using InheritsFrom = std::tuple<RANSNCModel, ZeroEqNC>; };
struct FlatPlateNCOneEq { using InheritsFrom = std::tuple<RANSNCModel, OneEqNC>; };
struct FlatPlateNCKOmega { using InheritsFrom = std::tuple<RANSNCModel, KOmegaNC>; };
struct FlatPlateNCLowReKEpsilon { using InheritsFrom = std::tuple<RANSNCModel, LowReKEpsilonNC>; };
struct FlatPlateNCKEpsilon { using InheritsFrom = std::tuple<RANSNCModel, KEpsilonNC>; };
// Isothermal Typetags
struct FlatPlateNCNIZeroEq { using InheritsFrom = std::tuple<RANSNCModel, ZeroEqNCNI>; };
struct FlatPlateNCNIOneEq { using InheritsFrom = std::tuple<RANSNCModel, OneEqNCNI>; };
struct FlatPlateNCNIKOmega { using InheritsFrom = std::tuple<RANSNCModel, KOmegaNCNI>; };
struct FlatPlateNCNILowReKEpsilon { using InheritsFrom = std::tuple<RANSNCModel, LowReKEpsilonNCNI>; };
struct FlatPlateNCNIKEpsilon { using InheritsFrom = std::tuple<RANSNCModel, KEpsilonNCNI>; };
} // end namespace TTag

// The fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::RANSNCModel>
{
  using H2OAir = FluidSystems::H2OAir<GetPropType<TypeTag, Properties::Scalar>>;
  static constexpr auto phaseIdx = H2OAir::gasPhaseIdx; // simulate the air phase
  using type = FluidSystems::OnePAdapter<H2OAir, phaseIdx>;
};

// replace the main component balance eq with a total balance eq
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::RANSNCModel> { static constexpr int value = 0; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::RANSNCModel>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::RANSNCModel> { using type = Dumux::FlatPlateNCTestProblem<TypeTag> ; };

template<class TypeTag>
struct EnableFVGridGeometryCache<TypeTag, TTag::RANSNCModel> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::RANSNCModel> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::RANSNCModel> { static constexpr bool value = true; };

template<class TypeTag>
struct UseMoles<TypeTag, TTag::RANSNCModel> { static constexpr bool value = true; };
} // end namespace Properties

/*!
 * \ingroup RANSNCTests
 * \brief  Test problem for the one-phase model.
 *
 * Dry air is entering from the left side and flows above a 1-D a flat plate.
 * In the middle of the inlet, water vapor is injected, which spreads by turbulent diffusion.
 * For the non-isothermal model the bottom has a constant temperature
 * which is \f$ \unit[30]{K} \f$ higher than the initial and inlet temperature.
 */
template <class TypeTag>
class FlatPlateNCTestProblem : public RANSProblem<TypeTag>
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

    static constexpr auto dimWorld = GetPropType<TypeTag, Properties::GridView>::dimensionworld;
    static constexpr auto transportEqIdx = Indices::conti0EqIdx + 1;
    static constexpr auto transportCompIdx = Indices::conti0EqIdx + 1;

public:
    FlatPlateNCTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry), eps_(1e-6)
    {
        inletVelocity_ = getParam<Scalar>("Problem.InletVelocity", 0.1);
        inletTemperature_ = getParam<Scalar>("Problem.InletTemperature", 283.15);
        wallTemperature_ = getParam<Scalar>("Problem.WallTemperature", 313.15);
        inletMoleFraction_ = getParam<Scalar>("Problem.InletMoleFraction", 1e-3);

        FluidSystem::init();
        Dumux::TurbulenceProperties<Scalar, dimWorld, true> turbulenceProperties;
        FluidState fluidState;
        const auto phaseIdx = 0;
        fluidState.setPressure(phaseIdx, 1e5);
        fluidState.setTemperature(temperature());
        fluidState.setMassFraction(phaseIdx, phaseIdx, 1.0);
        Scalar density = FluidSystem::density(fluidState, phaseIdx);
        Scalar kinematicViscosity = FluidSystem::viscosity(fluidState, phaseIdx) / density;
        Scalar diameter = this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1];
        viscosityTilde_ = 1e-3 * turbulenceProperties.viscosityTilde(inletVelocity_, diameter, kinematicViscosity);
        turbulentKineticEnergy_ = turbulenceProperties.turbulentKineticEnergy(inletVelocity_, diameter, kinematicViscosity);
        if (ModelTraits::turbulenceModel() == TurbulenceModel::komega)
            dissipation_ = turbulenceProperties.dissipationRate(inletVelocity_, diameter, kinematicViscosity);
        else
            dissipation_ = turbulenceProperties.dissipation(inletVelocity_, diameter, kinematicViscosity);
        turbulenceModelName_ = turbulenceModelToString(ModelTraits::turbulenceModel());
        std::cout << "Using the "<< turbulenceModelName_ << " Turbulence Model. \n";
        std::cout << std::endl;
    }

   /*!
     * \name Problem parameters
     */
    // \{

    bool isOnWallAtPos(const GlobalPosition& globalPos) const
    {
        return globalPos[1] < eps_;
    }


    bool shouldWriteRestartFile() const
    {
        return false;
    }

   /*!
     * \brief Returns the temperature within the domain in [K].
     *
     * The isothermal problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10C


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

        // turbulence model-specific boundary types
        static constexpr auto numEq = numTurbulenceEq(ModelTraits::turbulenceModel());
        setBcTypes_(values, globalPos, Dune::index_constant<numEq>{});

        if(isInlet_(globalPos))
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
            values.setDirichlet(transportCompIdx);
#if NONISOTHERMAL
            values.setDirichlet(Indices::temperatureIdx);
#endif
        }
        else if(isOutlet_(globalPos))
        {
            values.setDirichlet(Indices::pressureIdx);
            values.setOutflow(transportEqIdx);
#if NONISOTHERMAL
            values.setOutflow(Indices::energyEqIdx);
#endif
        }
        else if(isOnWallAtPos(globalPos))
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
            values.setNeumann(transportEqIdx);
#if NONISOTHERMAL
            values.setDirichlet(Indices::temperatureIdx);
#endif
        }
        else
            values.setAllSymmetry();

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

        if (isInlet_(globalPos)
            && globalPos[1] > 0.4 * this->fvGridGeometry().bBoxMax()[1]
            && globalPos[1] < 0.6 * this->fvGridGeometry().bBoxMax()[1])
        {
            values[transportCompIdx] = (time() > 10.0) ? inletMoleFraction_ : 0.0;
        }

#if NONISOTHERMAL
        if (time() > 10.0 && isOnWallAtPos(globalPos))
        {
            values[Indices::temperatureIdx] = wallTemperature_;
        }
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
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values[Indices::pressureIdx] = 1.0e+5;
        values[transportCompIdx] = 0.0;
#if NONISOTHERMAL
        values[Indices::temperatureIdx] = temperature();
#endif

        // block velocity profile
        values[Indices::velocityXIdx] = 0.0;
        if (!isOnWallAtPos(globalPos))
            values[Indices::velocityXIdx] =  inletVelocity_;
        values[Indices::velocityYIdx] = 0.0;

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
        return globalPos[0] < eps_;
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

    const Scalar eps_;
    Scalar inletVelocity_;
    Scalar inletMoleFraction_;
    Scalar inletTemperature_;
    Scalar wallTemperature_;
    Scalar viscosityTilde_;
    Scalar turbulentKineticEnergy_;
    Scalar dissipation_;
    TimeLoopPtr timeLoop_;
    std::string turbulenceModelName_;
};
} // end namespace Dumux

#endif
