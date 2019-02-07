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
 * \brief The free-flow RANS subproblem (staggered grid model).
 */

#ifndef DUMUX_RANS_SUBPROBLEM_HH
#define DUMUX_RANS_SUBPROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/material/fluidsystems/1padapter.hh>
#include <dumux/material/fluidsystems/h2oair.hh>

#include <dumux/freeflow/turbulencemodel.hh>
#include <dumux/freeflow/turbulenceproperties.hh>
#include <dumux/discretization/staggered/freeflow/properties.hh>

#include <dumux/freeflow/rans/problem.hh>

#include <dumux/freeflow/compositional/zeroeqncmodel.hh>

#include <dumux/freeflow/rans/oneeq/problem.hh>
#include <dumux/freeflow/compositional/oneeqncmodel.hh>

#include <dumux/freeflow/rans/twoeq/komega/problem.hh>
#include <dumux/freeflow/compositional/komegancmodel.hh>

#include <dumux/freeflow/rans/twoeq/lowrekepsilon/problem.hh>
#include <dumux/freeflow/compositional/lowrekepsilonncmodel.hh>

#include <dumux/freeflow/rans/twoeq/kepsilon/problem.hh>
#include <dumux/freeflow/compositional/kepsilonncmodel.hh>

#include <dumux/multidomain/boundary/stokesdarcy/couplingdata.hh>

namespace Dumux {
template <class TypeTag>
class RANSSubProblem;

namespace Properties {

// Create new type tags
namespace TTag {

// Base Typetag
struct RANSSubModel { using InheritsFrom = std::tuple<StaggeredFreeFlowModel>; };

// Typetags
struct RANSSubProblemZeroEqNCNI { using InheritsFrom = std::tuple<RANSSubModel, ZeroEqNCNI>; };
struct RANSSubProblemOneEqNCNI { using InheritsFrom = std::tuple<RANSSubModel, OneEqNCNI>; };
struct RANSSubProblemKOmegaNCNI { using InheritsFrom = std::tuple<RANSSubModel, KOmegaNCNI>; };
struct RANSSubProblemLowReKEpsilonNCNI { using InheritsFrom = std::tuple<RANSSubModel, LowReKEpsilonNCNI>; };
struct RANSSubProblemKEpsilonNCNI { using InheritsFrom = std::tuple<RANSSubModel, KEpsilonNCNI>; };

} // end namespace TTag

// The fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::RANSSubModel>
{
  using H2OAir = FluidSystems::H2OAir<GetPropType<TypeTag, Properties::Scalar>>;
  static constexpr auto phaseIdx = H2OAir::gasPhaseIdx;
  using type = FluidSystems::OnePAdapter<H2OAir, phaseIdx>;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::RANSSubModel>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::RANSSubModel> { static constexpr int value = 10; };

// Use formulation based on mole fractions
template<class TypeTag>
struct UseMoles<TypeTag, TTag::RANSSubModel> { static constexpr bool value = true; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::RANSSubModel> { using type = Dumux::RANSSubProblem<TypeTag> ; };

template<class TypeTag>
struct EnableFVGridGeometryCache<TypeTag, TTag::RANSSubModel> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::RANSSubModel> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::RANSSubModel> { static constexpr bool value = true; };
} // end namespace Properties

/*!
 * \ingroup BoundaryTests
 * \brief Test problem for the one-phase Reynolds Averaged (Navier-) Stokes problem.
 *
 * Horizontal flow from left to right with a turbulent velocity profile.
 */
template <class TypeTag>
class RANSSubProblem : public RANSProblem<TypeTag>
{
    using ParentType = RANSProblem<TypeTag>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

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
    using GridView = GetPropType<TypeTag, Properties::GridView>;

    using TimeLoopPtr = std::shared_ptr<TimeLoop<Scalar>>;
    using DiffusionCoefficientAveragingType = typename StokesDarcyCouplingOptions::DiffusionCoefficientAveragingType;

    static constexpr auto transportEqIdx = Indices::conti0EqIdx + 1;
    static constexpr auto transportCompIdx = Indices::conti0EqIdx + 1;
    static constexpr auto dimWorld = FVGridGeometry::GridView::dimensionworld;

public:
    RANSSubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(fvGridGeometry, "RANS"), eps_(1e-6), couplingManager_(couplingManager)
    {
        pressure_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InletPressure");
        inletMassFrac_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InletMassFrac");
        temperature_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.Temperature");
        inletVelocity_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InletVelocity");

        Dumux::TurbulenceProperties<Scalar, FVGridGeometry::GridView::dimensionworld, true> turbulenceProperties;
        FluidState fluidState;
        const auto phaseIdx = 0;
        fluidState.setPressure(phaseIdx, 1e5);
        fluidState.setTemperature(temperature());
        fluidState.setMassFraction(phaseIdx, phaseIdx, inletMassFrac());
        Scalar density = FluidSystem::density(fluidState, 0);
        Scalar kinematicViscosity = FluidSystem::viscosity(fluidState, 0) / density;
        Scalar charLength = this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1];
        diffCoeffAvgType_ = StokesDarcyCouplingOptions::stringToEnum(DiffusionCoefficientAveragingType{},
                                                                     getParam<std::string>("Problem.InterfaceDiffusionCoefficientAvg"));

        std::cout << "The reynolds Number (l_char = channel diameter) is: " << inletVelocity_ * charLength / kinematicViscosity << "\n";


        // ideally the viscosityTilde parameter as inflow for the Spalart-Allmaras model should be zero
        viscosityTilde_ = 1e-3 * turbulenceProperties.viscosityTilde(inletVelocity_, charLength, kinematicViscosity);
        turbulentKineticEnergy_ = turbulenceProperties.turbulentKineticEnergy(inletVelocity_, charLength, kinematicViscosity);

        if (ModelTraits::turbulenceModel() == TurbulenceModel::komega)
            dissipation_ = turbulenceProperties.dissipationRate(inletVelocity_, charLength, kinematicViscosity);
        else
            dissipation_ = turbulenceProperties.dissipation(inletVelocity_, charLength, kinematicViscosity);
        std::cout << std::endl;

        problemName_  =  getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
        turbulenceModelName_ = turbulenceModelToString(ModelTraits::turbulenceModel());
        std::cout << "Using the "<< turbulenceModelName_ << " Turbulence Model. \n";
    }

   /*!
     * \name Problem parameters
     */
    // \{

    //! \brief The problem name.
    const std::string& name() const
    { return problemName_; }

    //! \brief The name of the RANS model.
    const std::string& ransName() const
    { return turbulenceModelName_; }

    //! \brief Return the temperature within the domain in [K].
    const Scalar temperature() const
    { return temperature_; }

    //! \brief Returns the inlet velocity.
    const Scalar inletVelocity() const
    { return inletVelocity_ ;}

    //! \brief Returns the inlet pressure.
    const Scalar inletPressure() const
    { return pressure_; }

    //! \brief Returns the inlet mass fraction.
    const Scalar inletMassFrac() const
    { return inletMassFrac_; }

    // \}

    /*!
     * \name Boundary Locations
     */
    // \{

   /*!
     * \brief Returns if the location is on the wall.
     *
     * \param scvf The sub control volume face
     */
    bool isOnWall(const SubControlVolumeFace &scvf) const
    {
        GlobalPosition globalPos = scvf.ipGlobal();
        return isOnWallAtPos(globalPos);
    }

   /*!
     * \brief Returns if the location is on the wall.
     *
     * \param globalPos The global position
     */
    bool isOnWallAtPos(const GlobalPosition &globalPos) const
    {
        return onLowerBoundary_(globalPos);
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
     * \param element The finite element
     * \param scvf The sub control volume face
     */
    BoundaryTypes boundaryTypes(const Element& element,
                                const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes bTypes;
        GlobalPosition globalPos = scvf.ipGlobal();

        // turbulence model-specific boundary types
        static constexpr auto numEq = numTurbulenceEq(ModelTraits::turbulenceModel());
        setBcTypes_(bTypes, globalPos, Dune::index_constant<numEq>{});

        // common boundary types for all turbulence models
        if(isInlet_(globalPos))
        {
            bTypes.setDirichlet(Indices::velocityXIdx);
            bTypes.setDirichlet(Indices::velocityYIdx);
            bTypes.setDirichlet(Indices::temperatureIdx);
            bTypes.setDirichlet(transportCompIdx);
        }
        else if (isOutlet_(globalPos))
        {
            bTypes.setDirichlet(Indices::pressureIdx);
            bTypes.setOutflow(transportEqIdx);
            bTypes.setOutflow(Indices::energyEqIdx);
        }
        else if (isOnWall(scvf))
        {
            bTypes.setDirichlet(Indices::velocityXIdx);
            bTypes.setDirichlet(Indices::velocityYIdx);
            bTypes.setNeumann(Indices::energyEqIdx);
            bTypes.setNeumann(Indices::conti0EqIdx + 1);
        }
        else
        {
            bTypes.setAllSymmetry();
        }

        if (isOnCouplingWall_(scvf))
        {
            bTypes.setCouplingNeumann(Indices::conti0EqIdx);
            bTypes.setCouplingNeumann(Indices::conti0EqIdx + 1);
            bTypes.setCouplingNeumann(Indices::energyEqIdx);
            bTypes.setCouplingNeumann(scvf.directionIndex());
            bTypes.setBJS(1 - scvf.directionIndex());
        }


        return bTypes;
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
      * \brief Evaluate the boundary conditions for a dirichlet values at the boundary.
      *
      * \param element The finite element
      * \param scvf the sub control volume face
      */
    PrimaryVariables dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    {
        const auto globalPos = scvf.ipGlobal();
        PrimaryVariables values(initialAtPos(globalPos));

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeomentry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     * \param scvf The boundary sub control volume face
     */
    template<class ElementVolumeVariables, class ElementFaceVariables>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFaceVariables& elemFaceVars,
                        const SubControlVolumeFace& scvf) const
    {
        PrimaryVariables values(0.0);
        if(couplingManager().isCoupledEntity(CouplingManager::stokesIdx, scvf))
        {
            values[Indices::momentumYBalanceIdx] = couplingManager().couplingData().momentumCouplingCondition(element, fvGeometry, elemVolVars, elemFaceVars, scvf);

            const auto massFlux = couplingManager().couplingData().massCouplingCondition(element, fvGeometry, elemVolVars, elemFaceVars, scvf, diffCoeffAvgType_);
            values[Indices::conti0EqIdx] = massFlux[0];
            values[Indices::conti0EqIdx + 1] = massFlux[1];

            values[Indices::energyEqIdx] = couplingManager().couplingData().energyCouplingCondition(element, fvGeometry, elemVolVars, elemFaceVars, scvf, diffCoeffAvgType_);

        }
        return values;
    }


    // \}


   /*!
     * \name Volume terms
     */
    // \{

   /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);

        values[Indices::velocityXIdx] = inletVelocity();

        if (isOnWallAtPos(globalPos))
            values[Indices::velocityXIdx] = 0.0;

        FluidState fluidState;
        fluidState.setPressure(0, pressure_);
        fluidState.setTemperature(temperature());
        fluidState.setMoleFraction(0, 0, 1);
        Scalar density = FluidSystem::density(fluidState, 0);
        values[Indices::pressureIdx] = pressure_ - density * this->gravity()[dimWorld-1]* (this->fvGridGeometry().bBoxMax()[1] - globalPos[1]);

        values[transportCompIdx] = inletMassFrac();
        values[Indices::temperatureIdx] = temperature();

        // turbulence model-specific initial conditions
        static constexpr auto numEq = numTurbulenceEq(ModelTraits::turbulenceModel());
        setInitialAtPos_(values, globalPos, Dune::index_constant<numEq>{});

        return values;
    }


    /*!
     * \brief Returns the intrinsic permeability of required as input parameter
              for the Beavers-Joseph-Saffman boundary condition
     */
    Scalar permeability(const Element& element, const SubControlVolumeFace& scvf) const
    {
        return couplingManager().couplingData().darcyPermeability(element, scvf);
    }

    /*!
     * \brief Returns the alpha value required as input parameter for the
              Beavers-Joseph-Saffman boundary condition.
     */
    Scalar alphaBJ(const SubControlVolumeFace& scvf) const
    {
        return couplingManager().problem(CouplingManager::darcyIdx).spatialParams().beaversJosephCoeffAtPos(scvf.center());
    }

    void setTimeLoop(TimeLoopPtr timeLoop)
    {
        timeLoop_ = timeLoop;
    }

    Scalar time() const
    {
        return timeLoop_->time();
    }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    // \}

private:
    bool isInlet_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + eps_; }

    bool isOutlet_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->fvGridGeometry().bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->fvGridGeometry().bBoxMax()[1] - eps_; }

    bool isOnCouplingWall_(const SubControlVolumeFace& scvf) const
    { return couplingManager().isCoupledEntity(CouplingManager::stokesIdx, scvf);}


   /*!
     * \name Turbulent Boundary and Initial Conditions
     */
    // \{

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


    // \}

    Scalar inletVelocity_;
    Scalar pressure_;
    Scalar inletMassFrac_;
    Scalar temperature_;

    Scalar viscosityTilde_;
    Scalar turbulentKineticEnergy_;
    Scalar dissipation_;

    DiffusionCoefficientAveragingType diffCoeffAvgType_;
    Scalar eps_;
    std::string problemName_;
    std::string turbulenceModelName_;
    TimeLoopPtr timeLoop_;
    std::shared_ptr<CouplingManager> couplingManager_;
};
} // end namespace Dumux

#endif // DUMUX_RANS_SUBPROBLEM_HH
