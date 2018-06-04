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
 * \ingroup PorousmediumNonEquilibriumModel
 * \brief This class contains the volume variables required for the
 *        modules which require the specific interfacial area between
 *        fluid phases.
 *
 * This files contains all specializations which use 'real'
 * interfacial areas.
 */
#ifndef DUMUX_NONEQUILIBRIUM_VOLUME_VARIABLES_HH
#define DUMUX_NONEQUILIBRIUM_VOLUME_VARIABLES_HH

#include <cassert>

#include <dumux/common/dimensionlessnumbers.hh>

namespace Dumux {

/*!
 * \file
 * \ingroup PorousmediumNonEquilibriumModel
 * \brief This class contains the volume variables required for the
 *        modules which require the specific interfacial area between
 *        fluid phases.
 */
template<class Traits, class EquilibriumVolumeVariables, bool enableChemicalNonEquilibrium ,bool enableThermalNonEquilibrium>
class NonEquilibriumVolumeVariablesImplementation;

template<class Traits, class EquilibriumVolumeVariables>
using NonEquilibriumVolumeVariables =
      NonEquilibriumVolumeVariablesImplementation< Traits,
                                                   EquilibriumVolumeVariables,
                                                   Traits::ModelTraits::enableChemicalNonEquilibrium(),
                                                   Traits::ModelTraits::enableThermalNonEquilibrium() >;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// specialization for the case of NO kinetic mass but kinetic energy transfer of a fluid mixture and solid
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class Traits, class EquilibriumVolumeVariables>
class NonEquilibriumVolumeVariablesImplementation< Traits,
                                                   EquilibriumVolumeVariables,
                                                   false/*chemicalNonEquilibrium?*/,
                                                   true/*thermalNonEquilibrium?*/ >
    :public EquilibriumVolumeVariables
{
    using ParentType = EquilibriumVolumeVariables;
    using ParameterCache = typename Traits::FluidSystem::ParameterCache;
    using Scalar = typename Traits::PrimaryVariables::value_type;

    using ModelTraits = typename Traits::ModelTraits;
    using Indices = typename ModelTraits::Indices;
    static constexpr auto numEnergyEqFluid = ModelTraits::numEnergyEqFluid();
    static constexpr auto numEnergyEqSolid = ModelTraits::numEnergyEqSolid();

    static_assert((numEnergyEqFluid < 2), "This model is a specialization for a energy transfer of a fluid mixture and a solid");

    using DimLessNum = DimensionlessNumbers<Scalar>;

public:
     using FluidState = typename Traits::FluidState;
    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol &elemSol,
                const Problem &problem,
                const Element &element,
                const Scv& scv)
    {
        // Update parent type (also completes the fluid state)
        ParentType::update(elemSol, problem, element, scv);

        ParameterCache paramCache;
        paramCache.updateAll(this->fluidState());
        updateInterfacialArea(elemSol, this->fluidState(), paramCache, problem, element, scv);
    }

    /*!
     * \brief Updates the volume specific interfacial area [m^2 / m^3] between the phases.
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param fluidState Container for all the secondary variables concerning the fluids
     * \param paramCache The parameter cache corresponding to the fluid state
     * \param problem The problem to be solved
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElemSol, class Problem, class Element, class Scv>
    void updateInterfacialArea(const ElemSol& elemSol,
                               const FluidState& fluidState,
                               const ParameterCache& paramCache,
                               const Problem& problem,
                               const Element& element,
                               const Scv& scv)
    {
        factorMassTransfer_  = problem.spatialParams().factorMassTransfer(element, scv, elemSol);
        factorEnergyTransfer_ = problem.spatialParams().factorEnergyTransfer(element, scv, elemSol);
        characteristicLength_ = problem.spatialParams().characteristicLength(element, scv, elemSol);

        // set the dimensionless numbers and obtain respective quantities
        const unsigned int vIdxGlobal = scv.dofIndex();
        for (int phaseIdx = 0; phaseIdx < ModelTraits::numPhases(); ++phaseIdx)
        {
            const auto darcyMagVelocity     = problem.gridVariables().volumeDarcyMagVelocity(phaseIdx, vIdxGlobal);
            const auto dynamicViscosity     = fluidState.viscosity(phaseIdx);
            const auto density              = fluidState.density(phaseIdx);
            const auto kinematicViscosity   = dynamicViscosity/density;

            using FluidSystem = typename Traits::FluidSystem;
            const auto heatCapacity        = FluidSystem::heatCapacity(fluidState, paramCache, phaseIdx);
            const auto thermalConductivity = FluidSystem::thermalConductivity(fluidState, paramCache, phaseIdx);
            const auto porosity            = this->porosity();

            reynoldsNumber_[phaseIdx] = DimLessNum::reynoldsNumber(darcyMagVelocity, characteristicLength_, kinematicViscosity);
            prandtlNumber_[phaseIdx]  = DimLessNum::prandtlNumber(dynamicViscosity, heatCapacity, thermalConductivity);
            nusseltNumber_[phaseIdx]  = DimLessNum::nusseltNumberForced(reynoldsNumber_[phaseIdx],
                                                                        prandtlNumber_[phaseIdx],
                                                                        porosity,
                                                                        ModelTraits::nusseltFormulation());
        }
    }

    //! access function Reynolds Number
    const Scalar reynoldsNumber(const unsigned int phaseIdx) const { return reynoldsNumber_[phaseIdx]; }
    //! access function Prandtl Number
    const Scalar prandtlNumber(const unsigned int phaseIdx) const { return prandtlNumber_[phaseIdx]; }
    //! access function Nusselt Number
    const Scalar nusseltNumber(const unsigned int phaseIdx) const { return nusseltNumber_[phaseIdx]; }
    //! access function characteristic length
    const Scalar characteristicLength() const { return characteristicLength_; }
    //! access function pre factor energy transfer
    const Scalar factorEnergyTransfer() const { return factorEnergyTransfer_; }
    //! access function pre factor mass transfer
    const Scalar factorMassTransfer() const { return factorMassTransfer_; }

private:
    //! dimensionless numbers
    Scalar reynoldsNumber_[ModelTraits::numPhases()];
    Scalar prandtlNumber_[ModelTraits::numPhases()];
    Scalar nusseltNumber_[ModelTraits::numPhases()];

    Scalar characteristicLength_;
    Scalar factorEnergyTransfer_;
    Scalar factorMassTransfer_;
    Scalar solidSurface_ ;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
// specialization for the case of (only) kinetic mass transfer. Be careful, we do not test this!
////////////////////////////////////////////////////////////////////////////////////////////////////
template<class Traits, class EquilibriumVolumeVariables>
class NonEquilibriumVolumeVariablesImplementation< Traits,
                                                   EquilibriumVolumeVariables,
                                                   true/*chemicalNonEquilibrium?*/,
                                                   false/*thermalNonEquilibrium?*/>
    :public EquilibriumVolumeVariables
{
    using ParentType = EquilibriumVolumeVariables;
    using FluidState = typename Traits::FluidState;
    using FS = typename Traits::FluidSystem;
    using ParameterCache = typename Traits::FluidSystem::ParameterCache;
    using Scalar = typename Traits::PrimaryVariables::value_type;

    using ModelTraits = typename Traits::ModelTraits;
    using Indices = typename ModelTraits::Indices;

    static constexpr auto phase0Idx = FS::phase0Idx;
    static constexpr auto phase1Idx = FS::phase1Idx;
    static constexpr auto wCompIdx = FS::comp0Idx;
    static constexpr auto nCompIdx = FS::comp1Idx;

    using DimLessNum = DimensionlessNumbers<Scalar>;
public:
    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol &elemSol,
                const Problem &problem,
                const Element &element,
                const Scv& scv)
    {
        // Update parent type (also completes the fluid state)
        ParentType::update(elemSol, problem, element, scv);

        ParameterCache paramCache;
        paramCache.updateAll(this->fluidState());
        updateInterfacialArea(elemSol, this->fluidState(), paramCache, problem, element, scv);
    }

    /*!
     * \brief Updates the volume specific interfacial area [m^2 / m^3] between the phases.
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param fluidState Container for all the secondary variables concerning the fluids
     * \param paramCache The parameter cache corresponding to the fluid state
     * \param problem The problem to be solved
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElemSol, class Problem, class Element, class Scv>
    void updateInterfacialArea(const ElemSol& elemSol,
                               const FluidState& fluidState,
                               const ParameterCache& paramCache,
                               const Problem& problem,
                               const Element& element,
                               const Scv& scv)
    {
        // obtain parameters for awnsurface and material law
        const auto& awnSurfaceParams = problem.spatialParams().aWettingNonWettingSurfaceParams(element, scv, elemSol) ;
        const auto& materialParams  = problem.spatialParams().materialLawParams(element, scv, elemSol) ;

        const auto Sw = fluidState.saturation(phase0Idx) ;
        const auto pc = fluidState.pressure(phase1Idx) - fluidState.pressure(phase0Idx);

        // so far there is only a model for kinetic mass transfer between fluid phases
        using AwnSurface = typename Problem::SpatialParams::AwnSurface;
        interfacialArea_ = AwnSurface::interfacialArea(awnSurfaceParams, materialParams, Sw, pc);

        factorMassTransfer_ = problem.spatialParams().factorMassTransfer(element, scv, elemSol);
        characteristicLength_ = problem.spatialParams().characteristicLength(element, scv, elemSol);

        // set the dimensionless numbers and obtain respective quantities.
        const auto globalVertexIdx = problem.fvGridGeometry().vertexMapper().subIndex(element, scv, Element::Geometry::myDimension);
        for (int phaseIdx = 0; phaseIdx < ModelTraits::numPhases(); ++phaseIdx)
        {
            const auto darcyMagVelocity     = problem.gridVariables().volumeDarcyMagVelocity(phaseIdx, globalVertexIdx);
            const auto dynamicViscosity     = fluidState.viscosity(phaseIdx);
            const auto density              = fluidState.density(phaseIdx);
            const auto kinematicViscosity   = dynamicViscosity/density;

            // diffusion coefficient of non-wetting component in wetting phase
            using FluidSystem = typename Traits::FluidSystem;
            const auto diffCoeff = FluidSystem::binaryDiffusionCoefficient(fluidState,
                                                                           paramCache,
                                                                           phaseIdx,
                                                                           wCompIdx,
                                                                           nCompIdx);

            reynoldsNumber_[phaseIdx] = DimLessNum::reynoldsNumber(darcyMagVelocity, characteristicLength_, kinematicViscosity);
            schmidtNumber_[phaseIdx]  = DimLessNum::schmidtNumber(dynamicViscosity, density, diffCoeff);
            sherwoodNumber_[phaseIdx] = DimLessNum::sherwoodNumber(reynoldsNumber_[phaseIdx],
                                                                   schmidtNumber_[phaseIdx],
                                                                   ModelTraits::sherwoodFormulation());
        }
    }

    /*!
     * \brief The specific interfacial area between two fluid phases [m^2 / m^3]
     */
    const Scalar interfacialArea(const unsigned int phaseIIdx, const unsigned int phaseJIdx) const
    {
        // so far there is only a model for kinetic mass transfer between fluid phases
        assert( (phaseIIdx == phase1Idx && phaseJIdx == phase0Idx)
                || (phaseIIdx == phase0Idx && phaseJIdx == phase1Idx) );
        return interfacialArea_;
    }

    //! access function Reynolds Number
    const Scalar reynoldsNumber(const unsigned int phaseIdx) const { return reynoldsNumber_[phaseIdx]; }
    //! access function Schmidt Number
    const Scalar schmidtNumber(const unsigned int phaseIdx) const { return schmidtNumber_[phaseIdx]; }
    //! access function Sherwood Number
    const Scalar sherwoodNumber(const unsigned int phaseIdx) const { return sherwoodNumber_[phaseIdx]; }
    //! access function characteristic length
    const Scalar characteristicLength() const { return characteristicLength_; }
    //! access function pre factor mass transfer
    const Scalar factorMassTransfer() const { return factorMassTransfer_; }

private:
    Scalar characteristicLength_;
    Scalar factorMassTransfer_;
    Scalar solidSurface_ ;
    Scalar interfacialArea_ ;
    Scalar sherwoodNumber_[ModelTraits::numPhases()] ;
    Scalar schmidtNumber_[ModelTraits::numPhases()] ;
    Scalar reynoldsNumber_[ModelTraits::numPhases()] ;
};

// Specialization where everything is in non-equilibrium.
template<class Traits, class EquilibriumVolumeVariables>
class NonEquilibriumVolumeVariablesImplementation< Traits,
                                                   EquilibriumVolumeVariables,
                                                   true/*chemicalNonEquilibrium?*/,
                                                   true/*thermalNonEquilibrium?*/>
   :public EquilibriumVolumeVariables
{
    using ParentType = EquilibriumVolumeVariables;
    using FluidState = typename Traits::FluidState;
    using FS = typename Traits::FluidSystem;
    using ParameterCache = typename Traits::FluidSystem::ParameterCache;
    using Scalar = typename Traits::PrimaryVariables::value_type;

    using ModelTraits = typename Traits::ModelTraits;
    using Indices = typename ModelTraits::Indices;
    static constexpr auto numEnergyEqFluid = ModelTraits::numEnergyEqFluid();
    static constexpr auto numEnergyEqSolid = ModelTraits::numEnergyEqSolid();

    static constexpr auto phase0Idx = FS::phase0Idx;
    static constexpr auto phase1Idx = FS::phase1Idx;
    static constexpr auto sPhaseIdx = FS::numPhases;
    static constexpr auto wCompIdx = FS::comp0Idx;
    static constexpr auto nCompIdx = FS::comp1Idx;

    using DimLessNum = DimensionlessNumbers<Scalar>;
    static_assert((numEnergyEqFluid > 1), "This model only deals with energy transfer between two fluids and one solid phase");

public:
    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol &elemSol,
                const Problem &problem,
                const Element &element,
                const Scv& scv)
    {
        // Update parent type (also completes the fluid state)
        ParentType::update(elemSol, problem, element, scv);

        ParameterCache paramCache;
        paramCache.updateAll(this->fluidState());
        updateInterfacialArea(elemSol, this->fluidState(), paramCache, problem, element, scv);
    }

    /*!
     * \brief Updates the volume specific interfacial area [m^2 / m^3] between the phases.
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param fluidState Container for all the secondary variables concerning the fluids
     * \param paramCache The parameter cache corresponding to the fluid state
     * \param problem The problem to be solved
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElemSol, class Problem, class Element, class Scv>
    void updateInterfacialArea(const ElemSol& elemSol,
                               const FluidState& fluidState,
                               const ParameterCache& paramCache,
                               const Problem& problem,
                               const Element& element,
                               const Scv& scv)
    {
        // obtain (standard) material parameters (needed for the residual saturations)
       const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);

        //obtain parameters for interfacial area constitutive relations
        const auto& aWettingNonWettingSurfaceParams = problem.spatialParams().aWettingNonWettingSurfaceParams(element, scv, elemSol);
        const auto& aNonWettingSolidSurfaceParams = problem.spatialParams().aNonWettingSolidSurfaceParams(element, scv, elemSol);

        const Scalar pc = fluidState.pressure(phase1Idx) - fluidState.pressure(phase0Idx);
        const Scalar Sw = fluidState.saturation(phase0Idx);

        Scalar awn;

        // TODO can we delete this??? awn is overwritten anyway!
#define AwnRegul 0
        // This regularizes the interfacial area between the fluid phases.
        // This makes sure, that
        // a) some saturation cannot be lost: Never leave two phase region.
        // b) We cannot leave the fit region: no crazy (e.g. negative) values possible

        // const Scalar Swr =  aWettingNonWettingSurfaceParams.Swr() ;
        // const Scalar Snr =  aWettingNonWettingSurfaceParams.Snr() ;

        // this just leads to a stalling newton error as soon as this kicks in.
        // May be a spline or sth like this would help, but I do not which derivatives
        // to specify.
#if AwnRegul
        if(Sw < 5e-3 ) // or Sw > (1.-1e-5 )
            awn = 0. ; // 10.; //
#endif

        using AwnSurface = typename Problem::SpatialParams::AwnSurface;
        awn = AwnSurface::interfacialArea(aWettingNonWettingSurfaceParams, materialParams, Sw, pc );
        interfacialArea_[phase0Idx][phase1Idx] = awn;
        interfacialArea_[phase1Idx][phase0Idx] = interfacialArea_[phase0Idx][phase1Idx];
        interfacialArea_[phase0Idx][phase0Idx] = 0.;

        using AnsSurface = typename Problem::SpatialParams::AnsSurface;
        Scalar ans = AnsSurface::interfacialArea(aNonWettingSolidSurfaceParams, materialParams,Sw, pc);

        // Switch for using a a_{wn} relations that has some "maximum capillary pressure" as parameter.
        // That value is obtained by regularization of the pc(Sw) function.
#if USE_PCMAX
        const Scalar pcMax = problem.spatialParams().pcMax(element, scv, elemSol);
        // I know the solid surface from the pore network. But it is more consistent to use the fit value.
        using AnsSurface = typename Problem::SpatialParams::AnsSurface;
        solidSurface_ = AnsSurface::interfacialArea(aNonWettingSolidSurfaceParams, materialParams, /*Sw=*/0., pcMax);

        const Scalar aws = solidSurface_ - ans;
        interfacialArea_[phase0Idx][sPhaseIdx] = aws;
        interfacialArea_[sPhaseIdx][phase0Idx] = interfacialArea_[phase0Idx][sPhaseIdx];
        interfacialArea_[sPhaseIdx][sPhaseIdx] = 0.;
#else
        using AwsSurface = typename Problem::SpatialParams::AwsSurface;
        const auto& aWettingSolidSurfaceParams = problem.spatialParams().aWettingSolidSurfaceParams();
        const auto aws = AwsSurface::interfacialArea(aWettingSolidSurfaceParams,materialParams, Sw, pc );
        interfacialArea_[phase0Idx][sPhaseIdx] = aws ;
        interfacialArea_[sPhaseIdx][phase0Idx] = interfacialArea_[phase0Idx][sPhaseIdx];
        interfacialArea_[sPhaseIdx][sPhaseIdx] = 0.;
#endif

        interfacialArea_[phase1Idx][sPhaseIdx] = ans;
        interfacialArea_[sPhaseIdx][phase1Idx] = interfacialArea_[phase1Idx][sPhaseIdx];
        interfacialArea_[phase1Idx][phase1Idx] = 0.;

        factorMassTransfer_ = problem.spatialParams().factorMassTransfer(element, scv, elemSol);
        factorEnergyTransfer_ = problem.spatialParams().factorEnergyTransfer(element, scv, elemSol);
        characteristicLength_ = problem.spatialParams().characteristicLength(element, scv, elemSol);

        const auto vIdxGlobal = scv.dofIndex();
        using FluidSystem = typename Traits::FluidSystem;
        for (int phaseIdx = 0; phaseIdx < ModelTraits::numPhases(); ++phaseIdx)
        {
            const auto darcyMagVelocity    = problem.gridVariables().volumeDarcyMagVelocity(phaseIdx, vIdxGlobal);
            const auto dynamicViscosity    = fluidState.viscosity(phaseIdx);
            const auto density             = fluidState.density(phaseIdx);
            const auto kinematicViscosity  = dynamicViscosity/density;
            const auto heatCapacity        = FluidSystem::heatCapacity(fluidState, paramCache, phaseIdx);
            const auto thermalConductivity = FluidSystem::thermalConductivity(fluidState, paramCache, phaseIdx);

            // diffusion coefficient of non-wetting component in wetting phase
            const auto porosity = this->porosity();
            const auto diffCoeff = FluidSystem::binaryDiffusionCoefficient(fluidState,
                                                                           paramCache,
                                                                           phaseIdx,
                                                                           wCompIdx,
                                                                           nCompIdx);

            reynoldsNumber_[phaseIdx] = DimLessNum::reynoldsNumber(darcyMagVelocity, characteristicLength_, kinematicViscosity);
            prandtlNumber_[phaseIdx]  = DimLessNum::prandtlNumber(dynamicViscosity, heatCapacity, thermalConductivity);
            schmidtNumber_[phaseIdx]  = DimLessNum::schmidtNumber(dynamicViscosity, density, diffCoeff);
            nusseltNumber_[phaseIdx]  = DimLessNum::nusseltNumberForced(reynoldsNumber_[phaseIdx],
                                                                        prandtlNumber_[phaseIdx],
                                                                        porosity,
                                                                        ModelTraits::nusseltFormulation());
            // If Diffusion is not enabled, Sherwood is divided by zero
            sherwoodNumber_[phaseIdx] = DimLessNum::sherwoodNumber(reynoldsNumber_[phaseIdx],
                                                                   schmidtNumber_[phaseIdx],
                                                                   ModelTraits::sherwoodFormulation());
        }
    }

    /*!
     * \brief The specific interfacial area between two fluid phases [m^2 / m^3]
     * \note This is _only_ required by the kinetic mass/energy modules
     */
    const Scalar interfacialArea(const unsigned int phaseIIdx, const unsigned int phaseJIdx) const
    {
        // there is no interfacial area between a phase and itself
        assert(phaseIIdx not_eq phaseJIdx);
        return interfacialArea_[phaseIIdx][phaseJIdx];
    }

    //! access function Reynolds Number
    const Scalar reynoldsNumber(const unsigned int phaseIdx) const { return reynoldsNumber_[phaseIdx]; }
    //! access function Prandtl Number
    const Scalar prandtlNumber(const unsigned int phaseIdx) const { return prandtlNumber_[phaseIdx]; }
    //! access function Nusselt Number
    const Scalar nusseltNumber(const unsigned int phaseIdx) const { return nusseltNumber_[phaseIdx]; }
    //! access function Schmidt Number
    const Scalar schmidtNumber(const unsigned int phaseIdx) const { return schmidtNumber_[phaseIdx]; }
    //! access function Sherwood Number
    const Scalar sherwoodNumber(const unsigned int phaseIdx) const { return sherwoodNumber_[phaseIdx]; }
    //! access function characteristic length
    const Scalar characteristicLength() const { return characteristicLength_; }
    //! access function pre factor energy transfer
    const Scalar factorEnergyTransfer() const { return factorEnergyTransfer_; }
    //! access function pre factor mass transfer
    const Scalar factorMassTransfer() const { return factorMassTransfer_; }

private:
    //! dimensionless numbers
    Scalar reynoldsNumber_[ModelTraits::numPhases()];
    Scalar prandtlNumber_[ModelTraits::numPhases()];
    Scalar nusseltNumber_[ModelTraits::numPhases()];
    Scalar schmidtNumber_[ModelTraits::numPhases()];
    Scalar sherwoodNumber_[ModelTraits::numPhases()];
    Scalar characteristicLength_;
    Scalar factorEnergyTransfer_;
    Scalar factorMassTransfer_;
    Scalar solidSurface_ ;
    Scalar interfacialArea_[ModelTraits::numPhases()+numEnergyEqSolid][ModelTraits::numPhases()+numEnergyEqSolid];
};

} // namespace Dumux

#endif
