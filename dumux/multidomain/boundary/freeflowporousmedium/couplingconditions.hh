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
 * \ingroup FreeFlowPorousMediumCoupling
 * \copydoc Dumux::FreeFlowPorousMediumCouplingData
 */

#ifndef DUMUX_MD_FREEFLOW_POROUSMEDIUM_COUPLINGCONDITIONS_HH
#define DUMUX_MD_FREEFLOW_POROUSMEDIUM_COUPLINGCONDITIONS_HH

#include <numeric>

#include <dune/common/exceptions.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>

#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>

#include <dumux/flux/darcyslaw_fwd.hh>
#include <dumux/flux/fickslaw_fwd.hh>
#include <dumux/flux/forchheimerslaw_fwd.hh>

namespace Dumux {

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief This structs holds a set of options which allow to modify the Stokes-Darcy
 *        coupling mechanism during runtime.
 */
struct FreeFlowPorousMediumCouplingOptions
{
    /*!
     * \brief Defines which kind of averanging of diffusion coefficients
     *        (moleculat diffusion or thermal conductance) at the interface
     *        between free flow and porous medium shall be used.
     */
    enum class DiffusionCoefficientAveragingType
    {
        harmonic, arithmethic, ffOnly, pmOnly
    };

    /*!
     * \brief Convenience function to convert user input given as std::string to the corresponding enum class used for chosing the type
     *        of averaging of the diffusion/conduction parameter at the interface between the two domains.
     */
    static DiffusionCoefficientAveragingType stringToEnum(DiffusionCoefficientAveragingType, const std::string& diffusionCoefficientAveragingType)
    {
        if (diffusionCoefficientAveragingType == "Harmonic")
            return DiffusionCoefficientAveragingType::harmonic;
        else if (diffusionCoefficientAveragingType == "Arithmethic")
            return DiffusionCoefficientAveragingType::arithmethic;
        else if (diffusionCoefficientAveragingType == "FreeFlowOnly")
            return DiffusionCoefficientAveragingType::ffOnly;
        else if (diffusionCoefficientAveragingType == "PorousMediumOnly")
            return DiffusionCoefficientAveragingType::pmOnly;
        else
            DUNE_THROW(Dune::IOError, "Unknown DiffusionCoefficientAveragingType");
    }

};

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief This structs helps to check if the two sub models use the same fluidsystem.
 *        Specialization for the case of using an adapter only for the free-flow model.
 * \tparam FFFS The free-flow fluidsystem
 * \tparam PMFS The porous-medium flow fluidsystem
 */
template<class FFFS, class PMFS>
struct IsSameFluidSystem
{
    static_assert(FFFS::numPhases == 1, "Only single-phase fluidsystems may be used for free flow.");
    static constexpr bool value = std::is_same<typename FFFS::MultiPhaseFluidSystem, PMFS>::value;
};

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief This structs helps to check if the two sub models use the same fluidsystem.
 * \tparam FS The fluidsystem
 */
template<class FS>
struct IsSameFluidSystem<FS, FS>
{
    static_assert(FS::numPhases == 1, "Only single-phase fluidsystems may be used for free flow.");
    static constexpr bool value = std::is_same<FS, FS>::value; // always true
};

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief This structs indicates that Fick's law is not used for diffusion.
 * \tparam DiffLaw The diffusion law.
 */
template<class DiffLaw>
struct IsFicksLaw : public std::false_type {};

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief This structs indicates that Fick's law is used for diffusion.
 * \tparam DiffLaw The diffusion law.
 */
template<class T>
struct IsFicksLaw<Dumux::FicksLaw<T>> : public std::true_type {};

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief Helper struct to choose the correct index for phases and components. This is need if the porous-medium-flow model
          features more fluid phases than the free-flow model.
 * \tparam stokesIdx The domain index of the free-flow model.
 * \tparam porousMediumIndex The domain index of the porous-medium-flow model.
 * \tparam FFFS The free-flow fluidsystem.
 * \tparam hasAdapter Specifies whether an adapter class for the fluidsystem is used.
 */
template<std::size_t stokesIdx, std::size_t porousMediumIndex, class FFFS, bool hasAdapter>
struct IndexHelper;

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief Helper struct to choose the correct index for phases and components. This is need if the porous-medium-flow model
          features more fluid phases than the free-flow model. Specialization for the case that no adapter is used.
 * \tparam stokesIdx The domain index of the free-flow model.
 * \tparam porousMediumIndex The domain index of the porous-medium-flow model.
 * \tparam FFFS The free-flow fluidsystem.
 */
template<std::size_t stokesIdx, std::size_t porousMediumIndex, class FFFS>
struct IndexHelper<stokesIdx, porousMediumIndex, FFFS, false>
{
    /*!
     * \brief No adapter is used, just return the input index.
     */
    template<std::size_t i>
    static constexpr auto couplingPhaseIdx(Dune::index_constant<i>, int coupledPhaseIdx = 0)
    { return coupledPhaseIdx; }

    /*!
     * \brief No adapter is used, just return the input index.
     */
    template<std::size_t i>
    static constexpr auto couplingCompIdx(Dune::index_constant<i>, int coupledCompdIdx)
    { return coupledCompdIdx; }
};

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief Helper struct to choose the correct index for phases and components. This is need if the porous-medium-flow model
          features more fluid phases than the free-flow model. Specialization for the case that a adapter is used.
 * \tparam stokesIdx The domain index of the free-flow model.
 * \tparam porousMediumIndex The domain index of the porous-medium-flow model.
 * \tparam FFFS The free-flow fluidsystem.
 */
template<std::size_t stokesIdx, std::size_t porousMediumIndex, class FFFS>
struct IndexHelper<stokesIdx, porousMediumIndex, FFFS, true>
{
    /*!
     * \brief The free-flow model always uses phase index 0.
     */
    static constexpr int couplingPhaseIdx(Dune::index_constant<stokesIdx>, int coupledPhaseIdx = 0)
    { return 0; }

    /*!
     * \brief The phase index of the porous-medium-flow model is given by the adapter fluidsytem (i.e., user input).
     */
    static constexpr auto couplingPhaseIdx(Dune::index_constant<porousMediumIndex>, int coupledPhaseIdx = 0)
    { return FFFS::multiphaseFluidsystemPhaseIdx; }

    /*!
     * \brief The free-flow model does not need any change of the component index.
     */
    static constexpr auto couplingCompIdx(Dune::index_constant<stokesIdx>, int coupledCompdIdx)
    { return coupledCompdIdx; }

    /*!
     * \brief The component index of the porous-medium-flow model is mapped by the adapter fluidsytem.
     */
    static constexpr auto couplingCompIdx(Dune::index_constant<porousMediumIndex>, int coupledCompdIdx)
    { return FFFS::compIdx(coupledCompdIdx); }
};

template<class MDTraits, class CouplingManager, bool enableEnergyBalance, bool isCompositional>
class FreeFlowPorousMediumCouplingConditionsImplementation;

/*!
* \ingroup BoundaryCoupling
* \brief Data for the coupling of a Darcy model (cell-centered finite volume)
*        with a (Navier-)Stokes model (staggerd grid).
*/
template<class MDTraits, class CouplingManager>
using FreeFlowPorousMediumCouplingConditions
    = FreeFlowPorousMediumCouplingConditionsImplementation<
        MDTraits, CouplingManager,
        GetPropType<typename MDTraits::template SubDomain<0>::TypeTag, Properties::ModelTraits>::enableEnergyBalance(),
        (GetPropType<typename MDTraits::template SubDomain<0>::TypeTag, Properties::ModelTraits>::numFluidComponents() > 1)
    >;

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief A base class which provides some common methods used for Stokes-Darcy coupling.
 */
template<class MDTraits, class CouplingManager>
class FreeFlowPorousMediumCouplingConditionsImplementationBase
{
    using Scalar = typename MDTraits::Scalar;

    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using Element = typename GridGeometry<id>::GridView::template Codim<0>::Entity;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolumeFace = typename GridGeometry<id>::LocalView::SubControlVolumeFace;
    template<std::size_t id> using SubControlVolume = typename GridGeometry<id>::LocalView::SubControlVolume;
    template<std::size_t id> using Indices = typename GetPropType<SubDomainTypeTag<id>, Properties::ModelTraits>::Indices;
    template<std::size_t id> using ElementVolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::LocalView;
    template<std::size_t id> using VolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::VolumeVariables;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using FluidSystem = GetPropType<SubDomainTypeTag<id>, Properties::FluidSystem>;
    template<std::size_t id> using ModelTraits = GetPropType<SubDomainTypeTag<id>, Properties::ModelTraits>;
    template<std::size_t id> using GlobalPosition = typename Element<id>::Geometry::GlobalCoordinate;
    template<std::size_t id> using NumEqVector = typename Problem<id>::Traits::NumEqVector;

public:
    static constexpr auto freeFlowMomentumIndex = CouplingManager::freeFlowMomentumIndex;
    static constexpr auto freeFlowMassIndex = CouplingManager::freeFlowMassIndex;
    static constexpr auto porousMediumIndex = CouplingManager::porousMediumIndex;
private:

    using AdvectionType = GetPropType<SubDomainTypeTag<porousMediumIndex>, Properties::AdvectionType>;
    using DarcysLaw = Dumux::DarcysLaw<SubDomainTypeTag<porousMediumIndex>>;
    using ForchheimersLaw = Dumux::ForchheimersLaw<SubDomainTypeTag<porousMediumIndex>>;

    static constexpr bool adapterUsed = ModelTraits<porousMediumIndex>::numFluidPhases() > 1;
    using IndexHelper = Dumux::IndexHelper<freeFlowMassIndex, porousMediumIndex, FluidSystem<freeFlowMassIndex>, adapterUsed>;

    static constexpr int enableEnergyBalance = GetPropType<SubDomainTypeTag<freeFlowMassIndex>, Properties::ModelTraits>::enableEnergyBalance();
    static_assert(GetPropType<SubDomainTypeTag<porousMediumIndex>, Properties::ModelTraits>::enableEnergyBalance() == enableEnergyBalance,
                  "All submodels must both be either isothermal or non-isothermal");

    static_assert(IsSameFluidSystem<FluidSystem<freeFlowMassIndex>,
                                    FluidSystem<porousMediumIndex>>::value,
                  "All submodels must use the same fluid system");

    using DiffusionCoefficientAveragingType = typename FreeFlowPorousMediumCouplingOptions::DiffusionCoefficientAveragingType;

public:
    /*!
     * \brief Returns the corresponding phase index needed for coupling.
     */
    template<std::size_t i>
    static constexpr auto couplingPhaseIdx(Dune::index_constant<i> id, int coupledPhaseIdx = 0)
    { return IndexHelper::couplingPhaseIdx(id, coupledPhaseIdx); }

    /*!
     * \brief Returns the corresponding component index needed for coupling.
     */
    template<std::size_t i>
    static constexpr auto couplingCompIdx(Dune::index_constant<i> id, int coupledCompdIdx)
    { return IndexHelper::couplingCompIdx(id, coupledCompdIdx); }

    /*!
     * \brief Returns the intrinsic permeability of the coupled Darcy element.
     */
    template<class Context>
    static auto darcyPermeability(const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                                  const SubControlVolumeFace<freeFlowMomentumIndex>& scvf,
                                  const Context& context)
    { return context.volVars.permeability(); }

    /*!
     * \brief Returns the momentum flux across the coupling boundary.
     *
     * For the normal momentum coupling, the porous medium side of the coupling condition
     * is evaluated, i.e. -[p n]^pm.
     *
     */
    template<class Context>
    static NumEqVector<freeFlowMomentumIndex> momentumCouplingCondition(const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                                                                        const SubControlVolumeFace<freeFlowMomentumIndex>& scvf,
                                                                        const ElementVolumeVariables<freeFlowMomentumIndex>& elemVolVars,
                                                                        const Context& context)
    {
        static constexpr auto numPhasesDarcy = GetPropType<SubDomainTypeTag<porousMediumIndex>, Properties::ModelTraits>::numFluidPhases();
        NumEqVector<freeFlowMomentumIndex> momentumFlux(0.0);

        if (!scvf.isFrontal())
            return momentumFlux;

        const auto pmPhaseIdx = couplingPhaseIdx(porousMediumIndex);

        if(numPhasesDarcy > 1)
            momentumFlux[scvf.normalAxis()] = context.volVars.pressure(pmPhaseIdx);
        else // use pressure reconstruction for single phase models
            momentumFlux[scvf.normalAxis()] = pressureAtInterface_(fvGeometry, scvf, elemVolVars, context);

        // TODO: generalize for permeability tensors

        // normalize pressure
        momentumFlux[scvf.normalAxis()] -= elemVolVars.gridVolVars().problem().referencePressure(fvGeometry.element(), fvGeometry, scvf);

        momentumFlux[scvf.normalAxis()] *= scvf.directionSign();

        return momentumFlux;
    }

    /*!
     * \brief Evaluate an advective flux across the interface and consider upwinding.
     */
    static Scalar advectiveFlux(const Scalar insideQuantity, const Scalar outsideQuantity, const Scalar volumeFlow, bool insideIsUpstream)
    {
        const Scalar upwindWeight = 1.0; //TODO use Flux.UpwindWeight or something like Coupling.UpwindWeight?

        if(insideIsUpstream)
            return (upwindWeight * insideQuantity + (1.0 - upwindWeight) * outsideQuantity) * volumeFlow;
        else
            return (upwindWeight * outsideQuantity + (1.0 - upwindWeight) * insideQuantity) * volumeFlow;
    }

protected:

    /*!
     * \brief Returns the transmissibility used for either molecular diffusion or thermal conductivity.
     */
    template<std::size_t i, std::size_t j>
    Scalar transmissibility_(Dune::index_constant<i> domainI,
                             Dune::index_constant<j> domainJ,
                             const Scalar insideDistance,
                             const Scalar outsideDistance,
                             const Scalar avgQuantityI,
                             const Scalar avgQuantityJ,
                             const DiffusionCoefficientAveragingType diffCoeffAvgType) const
    {
        const Scalar totalDistance = insideDistance + outsideDistance;
        if(diffCoeffAvgType == DiffusionCoefficientAveragingType::harmonic)
        {
            return harmonicMean(avgQuantityI, avgQuantityJ, insideDistance, outsideDistance)
                   / totalDistance;
        }
        else if(diffCoeffAvgType == DiffusionCoefficientAveragingType::arithmethic)
        {
            return arithmeticMean(avgQuantityI, avgQuantityJ, insideDistance, outsideDistance)
                   / totalDistance;
        }
        else if(diffCoeffAvgType == DiffusionCoefficientAveragingType::ffOnly)
            return domainI == freeFlowMassIndex
                            ? avgQuantityI / totalDistance
                            : avgQuantityJ / totalDistance;

        else // diffCoeffAvgType == DiffusionCoefficientAveragingType::pmOnly)
            return domainI == porousMediumIndex
                            ? avgQuantityI / totalDistance
                            : avgQuantityJ / totalDistance;
    }

    /*!
     * \brief Returns the distance between an scvf and the corresponding scv center.
     */
    template<class Scv, class Scvf>
    Scalar getDistance_(const Scv& scv, const Scvf& scvf) const
    {
        return (scv.dofPosition() - scvf.ipGlobal()).two_norm();
    }

    /*!
     * \brief Returns the conductive energy flux acorss the interface.
     */
    template<std::size_t i, std::size_t j, bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar conductiveEnergyFlux_(Dune::index_constant<i> domainI,
                                 Dune::index_constant<j> domainJ,
                                 const FVElementGeometry<i>& fvGeometryI,
                                 const FVElementGeometry<j>& fvGeometryJ,
                                 const SubControlVolumeFace<i>& scvfI,
                                 const SubControlVolume<i>& scvI,
                                 const SubControlVolume<j>& scvJ,
                                 const VolumeVariables<i>& volVarsI,
                                 const VolumeVariables<j>& volVarsJ,
                                 const DiffusionCoefficientAveragingType diffCoeffAvgType) const
    {
        const Scalar insideDistance = getDistance_(scvI, scvfI);
        const Scalar outsideDistance = getDistance_(scvJ, scvfI);

        const Scalar deltaT = volVarsJ.temperature() - volVarsI.temperature();
        const Scalar tij = transmissibility_(
            domainI, domainJ,
            insideDistance, outsideDistance,
            volVarsI.effectiveThermalConductivity(), volVarsJ.effectiveThermalConductivity(),
            diffCoeffAvgType
        );

        return -tij * deltaT;
    }

    /*!
     * \brief Returns the pressure at the interface
     */
    template<class CouplingContext>
    static Scalar pressureAtInterface_(const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                                       const SubControlVolumeFace<freeFlowMomentumIndex>& scvf,
                                       const ElementVolumeVariables<freeFlowMomentumIndex>& elemVolVars,
                                       const CouplingContext& context)
    {
        GlobalPosition<freeFlowMomentumIndex> velocity(0.0);
        const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
        velocity[scv.dofAxis()] = elemVolVars[scv].velocity();
        const auto& pmScvf = context.fvGeometry.scvf(context.porousMediumScvfIdx);
        return computeCouplingPhasePressureAtInterface_(context.fvGeometry, pmScvf, context.volVars, context.gravity, velocity, AdvectionType());
    }

    /*!
     * \brief Returns the pressure at the interface using Forchheimers's law for reconstruction
     */
    static Scalar computeCouplingPhasePressureAtInterface_(const FVElementGeometry<porousMediumIndex>& fvGeometry,
                                                           const SubControlVolumeFace<porousMediumIndex>& scvf,
                                                           const VolumeVariables<porousMediumIndex>& volVars,
                                                           const typename Element<freeFlowMomentumIndex>::Geometry::GlobalCoordinate& couplingPhaseVelocity,
                                                           ForchheimersLaw)
    {
        DUNE_THROW(Dune::NotImplemented, "Forchheimer's law pressure reconstruction");
    }

    /*!
     * \brief Returns the pressure at the interface using Darcy's law for reconstruction
     */
    static Scalar computeCouplingPhasePressureAtInterface_(const FVElementGeometry<porousMediumIndex>& fvGeometry,
                                                           const SubControlVolumeFace<porousMediumIndex>& scvf,
                                                           const VolumeVariables<porousMediumIndex>& volVars,
                                                           const typename Element<freeFlowMomentumIndex>::Geometry::GlobalCoordinate& gravity,
                                                           const typename Element<freeFlowMomentumIndex>::Geometry::GlobalCoordinate& couplingPhaseVelocity,
                                                           DarcysLaw)
    {
        const auto pmPhaseIdx = couplingPhaseIdx(porousMediumIndex);
        const Scalar couplingPhaseCellCenterPressure = volVars.pressure(pmPhaseIdx);
        const Scalar couplingPhaseMobility = volVars.mobility(pmPhaseIdx);
        const Scalar couplingPhaseDensity = volVars.density(pmPhaseIdx);
        const auto K = volVars.permeability();

        // A tpfa approximation yields (works if mobility != 0)
        // v*n = -kr/mu*K * (gradP - rho*g)*n = mobility*(ti*(p_center - p_interface) + rho*n^TKg)
        // -> p_interface = (1/mobility * (-v*n) + rho*n^TKg)/ti + p_center
        // where v is the free-flow velocity (couplingPhaseVelocity)
        const auto alpha = vtmv(scvf.unitOuterNormal(), K, gravity);

        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto ti = computeTpfaTransmissibility(scvf, insideScv, K, 1.0);

        return (-1/couplingPhaseMobility * (scvf.unitOuterNormal() * couplingPhaseVelocity) + couplingPhaseDensity * alpha)/ti
               + couplingPhaseCellCenterPressure;
    }
};

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief Coupling data specialization for non-compositional models.
 */
template<class MDTraits, class CouplingManager, bool enableEnergyBalance>
class FreeFlowPorousMediumCouplingConditionsImplementation<MDTraits, CouplingManager, enableEnergyBalance, false>
: public FreeFlowPorousMediumCouplingConditionsImplementationBase<MDTraits, CouplingManager>
{
    using ParentType = FreeFlowPorousMediumCouplingConditionsImplementationBase<MDTraits, CouplingManager>;
    using Scalar = typename MDTraits::Scalar;

    // the sub domain type tags
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using Element = typename GridGeometry<id>::GridView::template Codim<0>::Entity;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolumeFace = typename GridGeometry<id>::LocalView::SubControlVolumeFace;
    template<std::size_t id> using SubControlVolume = typename GridGeometry<id>::LocalView::SubControlVolume;
    template<std::size_t id> using Indices = typename GetPropType<SubDomainTypeTag<id>, Properties::ModelTraits>::Indices;
    template<std::size_t id> using ElementVolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::LocalView;
    template<std::size_t id> using VolumeVariables  = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::VolumeVariables;

    static_assert(GetPropType<SubDomainTypeTag<ParentType::porousMediumIndex>, Properties::ModelTraits>::numFluidComponents() == GetPropType<SubDomainTypeTag<ParentType::porousMediumIndex>, Properties::ModelTraits>::numFluidPhases(),
                  "Darcy Model must not be compositional");

    using DiffusionCoefficientAveragingType = typename FreeFlowPorousMediumCouplingOptions::DiffusionCoefficientAveragingType;

public:
    using ParentType::ParentType;
    using ParentType::couplingPhaseIdx;

    /*!
     * \brief Returns the mass flux across the coupling boundary as seen from the Darcy domain.
     */
    template<std::size_t i, std::size_t j, class CouplingContext>
    static Scalar massCouplingCondition(Dune::index_constant<i> domainI,
                                        Dune::index_constant<j> domainJ,
                                        const FVElementGeometry<i>& fvGeometry,
                                        const SubControlVolumeFace<i>& scvf,
                                        const ElementVolumeVariables<i>& insideVolVars,
                                        const CouplingContext& context)
    {
        const Scalar normalVelocity = context.velocity * scvf.unitOuterNormal();
        const Scalar darcyDensity = insideVolVars[scvf.insideScvIdx()].density(couplingPhaseIdx(ParentType::porousMediumIndex));
        const Scalar stokesDensity = context.volVars.density();
        const bool insideIsUpstream = normalVelocity > 0.0;
        return ParentType::advectiveFlux(darcyDensity, stokesDensity, normalVelocity, insideIsUpstream);
    }

    /*!
     * \brief Returns the energy flux across the coupling boundary as seen from the Darcy domain.
     */
    template<bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar energyCouplingCondition(const Element<ParentType::porousMediumIndex>& element,
                                   const FVElementGeometry<ParentType::porousMediumIndex>& fvGeometry,
                                   const ElementVolumeVariables<ParentType::porousMediumIndex>& darcyElemVolVars,
                                   const SubControlVolumeFace<ParentType::porousMediumIndex>& scvf,
                                   const DiffusionCoefficientAveragingType diffCoeffAvgType = DiffusionCoefficientAveragingType::ffOnly) const
    {
        DUNE_THROW(Dune::NotImplemented, "Energy coupling condition");
    }

    /*!
     * \brief Returns the energy flux across the coupling boundary as seen from the free-flow domain.
     */
    template<bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar energyCouplingCondition(const Element<ParentType::freeFlowMassIndex>& element,
                                   const FVElementGeometry<ParentType::freeFlowMassIndex>& fvGeometry,
                                   const ElementVolumeVariables<ParentType::freeFlowMassIndex>& stokesElemVolVars,
                                   const SubControlVolumeFace<ParentType::freeFlowMassIndex>& scvf,
                                   const DiffusionCoefficientAveragingType diffCoeffAvgType = DiffusionCoefficientAveragingType::ffOnly) const
    {
        DUNE_THROW(Dune::NotImplemented, "Energy coupling condition");
    }

private:



    /*!
     * \brief Evaluate the energy flux across the interface.
     */
    template<std::size_t i, std::size_t j, bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar energyFlux_(Dune::index_constant<i> domainI,
                       Dune::index_constant<j> domainJ,
                       const FVElementGeometry<i>& insideFvGeometry,
                       const FVElementGeometry<j>& outsideFvGeometry,
                       const SubControlVolumeFace<i>& scvf,
                       const VolumeVariables<i>& insideVolVars,
                       const VolumeVariables<j>& outsideVolVars,
                       const Scalar velocity,
                       const bool insideIsUpstream,
                       const DiffusionCoefficientAveragingType diffCoeffAvgType) const
    {
        Scalar flux(0.0);

        const auto& insideScv = (*scvs(insideFvGeometry).begin());
        const auto& outsideScv = (*scvs(outsideFvGeometry).begin());

        // convective fluxes
        const Scalar insideTerm = insideVolVars.density(couplingPhaseIdx(domainI)) * insideVolVars.enthalpy(couplingPhaseIdx(domainI));
        const Scalar outsideTerm = outsideVolVars.density(couplingPhaseIdx(domainJ)) * outsideVolVars.enthalpy(couplingPhaseIdx(domainJ));

        flux += this->advectiveFlux(insideTerm, outsideTerm, velocity, insideIsUpstream);

        flux += this->conductiveEnergyFlux_(domainI, domainJ, insideFvGeometry, outsideFvGeometry, scvf, insideScv, outsideScv, insideVolVars, outsideVolVars, diffCoeffAvgType);

        return flux;
    }

};

} // end namespace Dumux

#endif
