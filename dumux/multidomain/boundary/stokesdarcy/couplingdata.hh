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
 * \ingroup StokesDarcyCoupling
 * \copydoc Dumux::StokesDarcyCouplingData
 */

#ifndef DUMUX_STOKES_DARCY_COUPLINGDATA_HH
#define DUMUX_STOKES_DARCY_COUPLINGDATA_HH

#include <numeric>

#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>
#include <dumux/flux/referencesystemformulation.hh>
#include <dumux/multidomain/couplingmanager.hh>

namespace Dumux {

/*!
 * \ingroup StokesDarcyCoupling
 * \brief This structs holds a set of options which allow to modify the Stokes-Darcy
 *        coupling mechanism during runtime.
 */
struct StokesDarcyCouplingOptions
{
    /*!
     * \brief Defines which kind of averanging of diffusion coefficiencients
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
 * \ingroup StokesDarcyCoupling
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
 * \ingroup StokesDarcyCoupling
 * \brief This structs helps to check if the two sub models use the same fluidsystem.
 * \tparam FS The fluidsystem
 */
template<class FS>
struct IsSameFluidSystem<FS, FS>
{
    static_assert(FS::numPhases == 1, "Only single-phase fluidsystems may be used for free flow.");
    static constexpr bool value = std::is_same<FS, FS>::value; // always true
};

// forward declaration
template <class TypeTag, DiscretizationMethod discMethod, ReferenceSystemFormulation referenceSystem>
class FicksLawImplementation;

/*!
 * \ingroup StokesDarcyCoupling
 * \brief This structs indicates that Fick's law is not used for diffusion.
 * \tparam DiffLaw The diffusion law.
 */
template<class DiffLaw>
struct IsFicksLaw : public std::false_type {};

/*!
 * \ingroup StokesDarcyCoupling
 * \brief This structs indicates that Fick's law is used for diffusion.
 * \tparam DiffLaw The diffusion law.
 */
template<class T, DiscretizationMethod discMethod, ReferenceSystemFormulation referenceSystem>
struct IsFicksLaw<FicksLawImplementation<T, discMethod, referenceSystem>> : public std::true_type {};

/*!
 * \ingroup StokesDarcyCoupling
 * \brief Helper struct to choose the correct index for phases and components. This is need if the porous-medium-flow model
          features more fluid phases than the free-flow model.
 * \tparam stokesIdx The domain index of the free-flow model.
 * \tparam darcyIdx The domain index of the porous-medium-flow model.
 * \tparam FFFS The free-flow fluidsystem.
 * \tparam hasAdapter Specifies whether an adapter class for the fluidsystem is used.
 */
template<std::size_t stokesIdx, std::size_t darcyIdx, class FFFS, bool hasAdapter>
struct IndexHelper;

/*!
 * \ingroup StokesDarcyCoupling
 * \brief Helper struct to choose the correct index for phases and components. This is need if the porous-medium-flow model
          features more fluid phases than the free-flow model. Specialization for the case that no adapter is used.
 * \tparam stokesIdx The domain index of the free-flow model.
 * \tparam darcyIdx The domain index of the porous-medium-flow model.
 * \tparam FFFS The free-flow fluidsystem.
 */
template<std::size_t stokesIdx, std::size_t darcyIdx, class FFFS>
struct IndexHelper<stokesIdx, darcyIdx, FFFS, false>
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
 * \ingroup StokesDarcyCoupling
 * \brief Helper struct to choose the correct index for phases and components. This is need if the porous-medium-flow model
          features more fluid phases than the free-flow model. Specialization for the case that a adapter is used.
 * \tparam stokesIdx The domain index of the free-flow model.
 * \tparam darcyIdx The domain index of the porous-medium-flow model.
 * \tparam FFFS The free-flow fluidsystem.
 */
template<std::size_t stokesIdx, std::size_t darcyIdx, class FFFS>
struct IndexHelper<stokesIdx, darcyIdx, FFFS, true>
{
    /*!
     * \brief The free-flow model always uses phase index 0.
     */
    static constexpr auto couplingPhaseIdx(Dune::index_constant<stokesIdx>, int coupledPhaseIdx = 0)
    { return 0; }

    /*!
     * \brief The phase index of the porous-medium-flow model is given by the adapter fluidsytem (i.e., user input).
     */
    static constexpr auto couplingPhaseIdx(Dune::index_constant<darcyIdx>, int coupledPhaseIdx = 0)
    { return FFFS::multiphaseFluidsystemPhaseIdx; }

    /*!
     * \brief The free-flow model does not need any change of the component index.
     */
    static constexpr auto couplingCompIdx(Dune::index_constant<stokesIdx>, int coupledCompdIdx)
    { return coupledCompdIdx; }

    /*!
     * \brief The component index of the porous-medium-flow model is mapped by the adapter fluidsytem.
     */
    static constexpr auto couplingCompIdx(Dune::index_constant<darcyIdx>, int coupledCompdIdx)
    { return FFFS::compIdx(coupledCompdIdx); }
};

//! forward declare
template <class TypeTag, DiscretizationMethod discMethod>
class DarcysLawImplementation;

//! forward declare
template <class TypeTag, DiscretizationMethod discMethod>
class ForchheimersLawImplementation;


template<class MDTraits, class CouplingManager, bool enableEnergyBalance, bool isCompositional>
class StokesDarcyCouplingDataImplementation;

/*!
* \ingroup BoundaryCoupling
* \brief Data for the coupling of a Darcy model (cell-centered finite volume)
*        with a (Navier-)Stokes model (staggerd grid).
*/
template<class MDTraits, class CouplingManager>
using StokesDarcyCouplingData = StokesDarcyCouplingDataImplementation<MDTraits, CouplingManager,
                                                                      GetPropType<typename MDTraits::template SubDomain<0>::TypeTag, Properties::ModelTraits>::enableEnergyBalance(),
                                                                      (GetPropType<typename MDTraits::template SubDomain<0>::TypeTag, Properties::ModelTraits>::numFluidComponents() > 1)>;

/*!
 * \ingroup StokesDarcyCoupling
 * \brief A base class which provides some common methods used for Stokes-Darcy coupling.
 */
template<class MDTraits, class CouplingManager>
class StokesDarcyCouplingDataImplementationBase
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
    static constexpr auto stokesIdx = CouplingManager::stokesIdx;
    static constexpr auto darcyIdx = CouplingManager::darcyIdx;

    using AdvectionType = GetPropType<SubDomainTypeTag<darcyIdx>, Properties::AdvectionType>;
    using DarcysLaw = DarcysLawImplementation<SubDomainTypeTag<darcyIdx>, GridGeometry<darcyIdx>::discMethod>;
    using ForchheimersLaw = ForchheimersLawImplementation<SubDomainTypeTag<darcyIdx>, GridGeometry<darcyIdx>::discMethod>;

    static constexpr bool adapterUsed = ModelTraits<darcyIdx>::numFluidPhases() > 1;
    using IndexHelper = Dumux::IndexHelper<stokesIdx, darcyIdx, FluidSystem<stokesIdx>, adapterUsed>;

    static constexpr int enableEnergyBalance = GetPropType<SubDomainTypeTag<stokesIdx>, Properties::ModelTraits>::enableEnergyBalance();
    static_assert(GetPropType<SubDomainTypeTag<darcyIdx>, Properties::ModelTraits>::enableEnergyBalance() == enableEnergyBalance,
                  "All submodels must both be either isothermal or non-isothermal");

    static_assert(IsSameFluidSystem<FluidSystem<stokesIdx>,
                                    FluidSystem<darcyIdx>>::value,
                  "All submodels must use the same fluid system");

    using DiffusionCoefficientAveragingType = typename StokesDarcyCouplingOptions::DiffusionCoefficientAveragingType;

public:
    StokesDarcyCouplingDataImplementationBase(const CouplingManager& couplingmanager): couplingManager_(couplingmanager) {}

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
     * \brief Returns a reference to the coupling manager.
     */
    const CouplingManager& couplingManager() const
    { return couplingManager_; }

    /*!
     * \brief Returns the intrinsic permeability of the coupled Darcy element.
     */
    auto darcyPermeability(const Element<stokesIdx>& element, const SubControlVolumeFace<stokesIdx>& scvf) const
    {
        const auto& stokesContext = couplingManager().stokesCouplingContext(element, scvf);
        return stokesContext.volVars.permeability();
    }

    /*!
     * \brief Returns the momentum flux across the coupling boundary.
     *
     * For the normal momentum coupling, the porous medium side of the coupling condition
     * is evaluated, i.e. -[p n]^pm.
     *
     */
    template<class ElementFaceVariables>
    Scalar momentumCouplingCondition(const Element<stokesIdx>& element,
                                     const FVElementGeometry<stokesIdx>& fvGeometry,
                                     const ElementVolumeVariables<stokesIdx>& stokesElemVolVars,
                                     const ElementFaceVariables& stokesElemFaceVars,
                                     const SubControlVolumeFace<stokesIdx>& scvf) const
    {
        static constexpr auto numPhasesDarcy = GetPropType<SubDomainTypeTag<darcyIdx>, Properties::ModelTraits>::numFluidPhases();

        Scalar momentumFlux(0.0);
        const auto& stokesContext = couplingManager_.stokesCouplingContext(element, scvf);
        const auto darcyPhaseIdx = couplingPhaseIdx(darcyIdx);

        // - p_pm * n_pm = p_pm * n_ff
        const Scalar darcyPressure = stokesContext.volVars.pressure(darcyPhaseIdx);

        if(numPhasesDarcy > 1)
            momentumFlux = darcyPressure;
        else // use pressure reconstruction for single phase models
            momentumFlux = pressureAtInterface_(element, scvf, stokesElemFaceVars, stokesContext);
        // TODO: generalize for permeability tensors

        // normalize pressure
        if(getPropValue<SubDomainTypeTag<stokesIdx>, Properties::NormalizePressure>())
            momentumFlux -= couplingManager_.problem(stokesIdx).initial(scvf)[Indices<stokesIdx>::pressureIdx];

        momentumFlux *= scvf.directionSign();

        return momentumFlux;
    }

    /*!
     * \brief Evaluate an advective flux across the interface and consider upwinding.
     */
    Scalar advectiveFlux(const Scalar insideQuantity, const Scalar outsideQuantity, const Scalar volumeFlow, bool insideIsUpstream) const
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
            return domainI == stokesIdx
                            ? avgQuantityI / totalDistance
                            : avgQuantityJ / totalDistance;

        else // diffCoeffAvgType == DiffusionCoefficientAveragingType::pmOnly)
            return domainI == darcyIdx
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
        const Scalar tij = transmissibility_(domainI,
                                             domainJ,
                                             insideDistance,
                                             outsideDistance,
                                             volVarsI.effectiveThermalConductivity(),
                                             volVarsJ.effectiveThermalConductivity(),
                                             diffCoeffAvgType);

        return -tij * deltaT;
    }

    /*!
     * \brief Returns the pressure at the interface
     */
    template<class ElementFaceVariables, class CouplingContext>
    Scalar pressureAtInterface_(const Element<stokesIdx>& element,
                                const SubControlVolumeFace<stokesIdx>& scvf,
                                const ElementFaceVariables& elemFaceVars,
                                const CouplingContext& context) const
    {
        GlobalPosition<stokesIdx> velocity(0.0);
        velocity[scvf.directionIndex()] = elemFaceVars[scvf].velocitySelf();
        const auto& darcyScvf = context.fvGeometry.scvf(context.darcyScvfIdx);
        return computeCouplingPhasePressureAtInterface_(context.element, context.fvGeometry, darcyScvf, context.volVars, velocity, AdvectionType());
    }

    /*!
     * \brief Returns the pressure at the interface using Forchheimers's law for reconstruction
     */
     Scalar computeCouplingPhasePressureAtInterface_(const Element<darcyIdx>& element,
                                                     const FVElementGeometry<darcyIdx>& fvGeometry,
                                                     const SubControlVolumeFace<darcyIdx>& scvf,
                                                     const VolumeVariables<darcyIdx>& volVars,
                                                     const typename Element<stokesIdx>::Geometry::GlobalCoordinate& couplingPhaseVelocity,
                                                     ForchheimersLaw) const
    {
        const auto darcyPhaseIdx = couplingPhaseIdx(darcyIdx);
        const Scalar cellCenterPressure = volVars.pressure(darcyPhaseIdx);
        using std::sqrt;

        // v + (cF*sqrt(K)*rho/mu*|v|) * v  = - K/mu grad(p - rho g)
        // multiplying with n and using a tpfa for the right-hand side yields
        // v*n + (cF*sqrt(K)*rho/mu*|v|) * (v*n) =  1/mu * (ti*(p_center - p_interface) + rho*n^TKg)
        // --> p_interface = (-mu*v*n + (cF*sqrt(K)*rho*|v|) * (-v*n) + rho*n^TKg)/ti + p_center
        const auto velocity = couplingPhaseVelocity;
        const Scalar mu = volVars.viscosity(darcyPhaseIdx);
        const Scalar rho = volVars.density(darcyPhaseIdx);
        const auto K = volVars.permeability();
        const auto alpha = vtmv(scvf.unitOuterNormal(), K, couplingManager_.problem(darcyIdx).spatialParams().gravity(scvf.center()));

        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto ti = computeTpfaTransmissibility(scvf, insideScv, K, 1.0);

        // get the Forchheimer coefficient
        Scalar cF = couplingManager_.problem(darcyIdx).spatialParams().forchCoeff(scvf);

        const Scalar interfacePressure = ((-mu*(scvf.unitOuterNormal() * velocity))
                                        + (-(scvf.unitOuterNormal() * velocity) * velocity.two_norm() * rho * sqrt(darcyPermeability(element, scvf)) * cF)
                                        +  rho * alpha)/ti + cellCenterPressure;
        return interfacePressure;
    }

    /*!
     * \brief Returns the pressure at the interface using Darcy's law for reconstruction
     */
    Scalar computeCouplingPhasePressureAtInterface_(const Element<darcyIdx>& element,
                                                    const FVElementGeometry<darcyIdx>& fvGeometry,
                                                    const SubControlVolumeFace<darcyIdx>& scvf,
                                                    const VolumeVariables<darcyIdx>& volVars,
                                                    const typename Element<stokesIdx>::Geometry::GlobalCoordinate& couplingPhaseVelocity,
                                                    DarcysLaw) const
    {
        const auto darcyPhaseIdx = couplingPhaseIdx(darcyIdx);
        const Scalar couplingPhaseCellCenterPressure = volVars.pressure(darcyPhaseIdx);
        const Scalar couplingPhaseMobility = volVars.mobility(darcyPhaseIdx);
        const Scalar couplingPhaseDensity = volVars.density(darcyPhaseIdx);
        const auto K = volVars.permeability();

        // A tpfa approximation yields (works if mobility != 0)
        // v*n = -kr/mu*K * (gradP - rho*g)*n = mobility*(ti*(p_center - p_interface) + rho*n^TKg)
        // -> p_interface = (1/mobility * (-v*n) + rho*n^TKg)/ti + p_center
        // where v is the free-flow velocity (couplingPhaseVelocity)
        const auto alpha = vtmv(scvf.unitOuterNormal(), K, couplingManager_.problem(darcyIdx).spatialParams().gravity(scvf.center()));

        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto ti = computeTpfaTransmissibility(scvf, insideScv, K, 1.0);

        return (-1/couplingPhaseMobility * (scvf.unitOuterNormal() * couplingPhaseVelocity) + couplingPhaseDensity * alpha)/ti
               + couplingPhaseCellCenterPressure;
    }


private:
    const CouplingManager& couplingManager_;

};

/*!
 * \ingroup StokesDarcyCoupling
 * \brief Coupling data specialization for non-compositional models.
 */
template<class MDTraits, class CouplingManager, bool enableEnergyBalance>
class StokesDarcyCouplingDataImplementation<MDTraits, CouplingManager, enableEnergyBalance, false>
: public StokesDarcyCouplingDataImplementationBase<MDTraits, CouplingManager>
{
    using ParentType = StokesDarcyCouplingDataImplementationBase<MDTraits, CouplingManager>;
    using Scalar = typename MDTraits::Scalar;
    static constexpr auto stokesIdx = CouplingManager::stokesIdx;
    static constexpr auto darcyIdx = CouplingManager::darcyIdx;
    static constexpr auto stokesFaceIdx = CouplingManager::stokesFaceIdx;
    static constexpr auto stokesCellCenterIdx = CouplingManager::stokesCellCenterIdx;

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
    template<std::size_t id> using ElementFaceVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridFaceVariables>::LocalView;
    template<std::size_t id> using VolumeVariables  = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::VolumeVariables;

    static_assert(GetPropType<SubDomainTypeTag<darcyIdx>, Properties::ModelTraits>::numFluidComponents() == GetPropType<SubDomainTypeTag<darcyIdx>, Properties::ModelTraits>::numFluidPhases(),
                  "Darcy Model must not be compositional");

    using DiffusionCoefficientAveragingType = typename StokesDarcyCouplingOptions::DiffusionCoefficientAveragingType;

public:
    using ParentType::ParentType;
    using ParentType::couplingPhaseIdx;

    /*!
     * \brief Returns the mass flux across the coupling boundary as seen from the Darcy domain.
     */
    Scalar massCouplingCondition(const Element<darcyIdx>& element,
                                 const FVElementGeometry<darcyIdx>& fvGeometry,
                                 const ElementVolumeVariables<darcyIdx>& darcyElemVolVars,
                                 const SubControlVolumeFace<darcyIdx>& scvf) const
    {
        const auto& darcyContext = this->couplingManager().darcyCouplingContext(element, scvf);
        const Scalar velocity = darcyContext.velocity * scvf.unitOuterNormal();
        const Scalar darcyDensity = darcyElemVolVars[scvf.insideScvIdx()].density(couplingPhaseIdx(darcyIdx));
        const Scalar stokesDensity = darcyContext.volVars.density();
        const bool insideIsUpstream = velocity > 0.0;

        return massFlux_(velocity, darcyDensity, stokesDensity, insideIsUpstream);
    }

    /*!
     * \brief Returns the mass flux across the coupling boundary as seen from the free-flow domain.
     */
    Scalar massCouplingCondition(const Element<stokesIdx>& element,
                                 const FVElementGeometry<stokesIdx>& fvGeometry,
                                 const ElementVolumeVariables<stokesIdx>& stokesElemVolVars,
                                 const ElementFaceVariables<stokesIdx>& stokesElemFaceVars,
                                 const SubControlVolumeFace<stokesIdx>& scvf) const
    {
        const auto& stokesContext = this->couplingManager().stokesCouplingContext(element, scvf);
        const Scalar velocity = stokesElemFaceVars[scvf].velocitySelf();
        const Scalar stokesDensity = stokesElemVolVars[scvf.insideScvIdx()].density();
        const Scalar darcyDensity = stokesContext.volVars.density(couplingPhaseIdx(darcyIdx));
        const bool insideIsUpstream = sign(velocity) == scvf.directionSign();

        return massFlux_(velocity * scvf.directionSign(), stokesDensity, darcyDensity, insideIsUpstream);
    }

    /*!
     * \brief Returns the energy flux across the coupling boundary as seen from the Darcy domain.
     */
    template<bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar energyCouplingCondition(const Element<darcyIdx>& element,
                                   const FVElementGeometry<darcyIdx>& fvGeometry,
                                   const ElementVolumeVariables<darcyIdx>& darcyElemVolVars,
                                   const SubControlVolumeFace<darcyIdx>& scvf,
                                   const DiffusionCoefficientAveragingType diffCoeffAvgType = DiffusionCoefficientAveragingType::ffOnly) const
    {
        const auto& darcyContext = this->couplingManager().darcyCouplingContext(element, scvf);
        const auto& darcyVolVars = darcyElemVolVars[scvf.insideScvIdx()];
        const auto& stokesVolVars = darcyContext.volVars;

        const Scalar velocity = darcyContext.velocity * scvf.unitOuterNormal();
        const bool insideIsUpstream = velocity > 0.0;

        return energyFlux_(darcyIdx, stokesIdx, fvGeometry, darcyContext.fvGeometry, scvf,
                           darcyVolVars, stokesVolVars, velocity, insideIsUpstream, diffCoeffAvgType);
    }

    /*!
     * \brief Returns the energy flux across the coupling boundary as seen from the free-flow domain.
     */
    template<bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar energyCouplingCondition(const Element<stokesIdx>& element,
                                   const FVElementGeometry<stokesIdx>& fvGeometry,
                                   const ElementVolumeVariables<stokesIdx>& stokesElemVolVars,
                                   const ElementFaceVariables<stokesIdx>& stokesElemFaceVars,
                                   const SubControlVolumeFace<stokesIdx>& scvf,
                                   const DiffusionCoefficientAveragingType diffCoeffAvgType = DiffusionCoefficientAveragingType::ffOnly) const
    {
        const auto& stokesContext = this->couplingManager().stokesCouplingContext(element, scvf);
        const auto& stokesVolVars = stokesElemVolVars[scvf.insideScvIdx()];
        const auto& darcyVolVars = stokesContext.volVars;

        const Scalar velocity = stokesElemFaceVars[scvf].velocitySelf();
        const bool insideIsUpstream = sign(velocity) == scvf.directionSign();

        return energyFlux_(stokesIdx, darcyIdx, fvGeometry, stokesContext.fvGeometry, scvf,
                           stokesVolVars, darcyVolVars, velocity * scvf.directionSign(), insideIsUpstream, diffCoeffAvgType);
    }

private:

    /*!
     * \brief Evaluate the mole/mass flux across the interface.
     */
    Scalar massFlux_(const Scalar velocity,
                     const Scalar insideDensity,
                     const Scalar outSideDensity,
                     bool insideIsUpstream) const
    {
        return this->advectiveFlux(insideDensity, outSideDensity, velocity, insideIsUpstream);
    }

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

/*!
 * \ingroup StokesDarcyCoupling
 * \brief Coupling data specialization for compositional models.
 */
template<class MDTraits, class CouplingManager, bool enableEnergyBalance>
class StokesDarcyCouplingDataImplementation<MDTraits, CouplingManager, enableEnergyBalance, true>
: public StokesDarcyCouplingDataImplementationBase<MDTraits, CouplingManager>
{
    using ParentType = StokesDarcyCouplingDataImplementationBase<MDTraits, CouplingManager>;
    using Scalar = typename MDTraits::Scalar;
    static constexpr auto stokesIdx = CouplingManager::stokesIdx;
    static constexpr auto darcyIdx = CouplingManager::darcyIdx;
    static constexpr auto stokesFaceIdx = CouplingManager::stokesFaceIdx;
    static constexpr auto stokesCellCenterIdx = CouplingManager::stokesCellCenterIdx;

    // the sub domain type tags
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using Element = typename GridGeometry<id>::GridView::template Codim<0>::Entity;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolumeFace = typename FVElementGeometry<id>::SubControlVolumeFace;
    template<std::size_t id> using SubControlVolume = typename GridGeometry<id>::LocalView::SubControlVolume;
    template<std::size_t id> using Indices = typename GetPropType<SubDomainTypeTag<id>, Properties::ModelTraits>::Indices;
    template<std::size_t id> using ElementVolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::LocalView;
    template<std::size_t id> using ElementFaceVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridFaceVariables>::LocalView;
    template<std::size_t id> using VolumeVariables  = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::VolumeVariables;
    template<std::size_t id> using FluidSystem  = GetPropType<SubDomainTypeTag<id>, Properties::FluidSystem>;

    static constexpr auto numComponents = GetPropType<SubDomainTypeTag<stokesIdx>, Properties::ModelTraits>::numFluidComponents();
    static constexpr auto replaceCompEqIdx = GetPropType<SubDomainTypeTag<stokesIdx>, Properties::ModelTraits>::replaceCompEqIdx();
    static constexpr bool useMoles = GetPropType<SubDomainTypeTag<stokesIdx>, Properties::ModelTraits>::useMoles();
    static constexpr auto referenceSystemFormulation = GetPropType<SubDomainTypeTag<stokesIdx>, Properties::MolecularDiffusionType>::referenceSystemFormulation();

    static_assert(GetPropType<SubDomainTypeTag<darcyIdx>, Properties::ModelTraits>::numFluidComponents() == numComponents, "Both submodels must use the same number of components");
    static_assert(getPropValue<SubDomainTypeTag<darcyIdx>, Properties::UseMoles>() == useMoles, "Both submodels must either use moles or not");
    static_assert(getPropValue<SubDomainTypeTag<darcyIdx>, Properties::ReplaceCompEqIdx>() == replaceCompEqIdx, "Both submodels must use the same replaceCompEqIdx");
    static_assert(GetPropType<SubDomainTypeTag<darcyIdx>, Properties::MolecularDiffusionType>::referenceSystemFormulation() == referenceSystemFormulation,
                  "Both submodels must use the same reference system formulation for diffusion");

    using NumEqVector = Dune::FieldVector<Scalar, numComponents>;

    using DiffusionCoefficientAveragingType = typename StokesDarcyCouplingOptions::DiffusionCoefficientAveragingType;

    static constexpr bool isFicksLaw = IsFicksLaw<GetPropType<SubDomainTypeTag<stokesIdx>, Properties::MolecularDiffusionType>>();
    static_assert(isFicksLaw == IsFicksLaw<GetPropType<SubDomainTypeTag<darcyIdx>, Properties::MolecularDiffusionType>>(),
                  "Both submodels must use the same diffusion law.");

    using ReducedComponentVector = Dune::FieldVector<Scalar, numComponents-1>;
    using ReducedComponentMatrix = Dune::FieldMatrix<Scalar, numComponents-1, numComponents-1>;

    using MolecularDiffusionType = GetPropType<SubDomainTypeTag<stokesIdx>, Properties::MolecularDiffusionType>;
public:
    using ParentType::ParentType;
    using ParentType::couplingPhaseIdx;
    using ParentType::couplingCompIdx;

    /*!
     * \brief Returns the mass flux across the coupling boundary as seen from the Darcy domain.
     */
    NumEqVector massCouplingCondition(const Element<darcyIdx>& element,
                                      const FVElementGeometry<darcyIdx>& fvGeometry,
                                      const ElementVolumeVariables<darcyIdx>& darcyElemVolVars,
                                      const SubControlVolumeFace<darcyIdx>& scvf,
                                      const DiffusionCoefficientAveragingType diffCoeffAvgType = DiffusionCoefficientAveragingType::ffOnly) const
    {
        NumEqVector flux(0.0);
        const auto& darcyContext = this->couplingManager().darcyCouplingContext(element, scvf);
        const auto& darcyVolVars = darcyElemVolVars[scvf.insideScvIdx()];
        const auto& stokesVolVars = darcyContext.volVars;
        const auto& outsideScv = (*scvs(darcyContext.fvGeometry).begin());

        const Scalar velocity = darcyContext.velocity * scvf.unitOuterNormal();
        const bool insideIsUpstream = velocity > 0.0;

        return massFlux_(darcyIdx, stokesIdx, fvGeometry,
                         scvf, darcyVolVars, stokesVolVars,
                         outsideScv, velocity, insideIsUpstream,
                         diffCoeffAvgType);
    }

    /*!
     * \brief Returns the mass flux across the coupling boundary as seen from the free-flow domain.
     */
    NumEqVector massCouplingCondition(const Element<stokesIdx>& element,
                                      const FVElementGeometry<stokesIdx>& fvGeometry,
                                      const ElementVolumeVariables<stokesIdx>& stokesElemVolVars,
                                      const ElementFaceVariables<stokesIdx>& stokesElemFaceVars,
                                      const SubControlVolumeFace<stokesIdx>& scvf,
                                      const DiffusionCoefficientAveragingType diffCoeffAvgType = DiffusionCoefficientAveragingType::ffOnly) const
    {
        NumEqVector flux(0.0);
        const auto& stokesContext = this->couplingManager().stokesCouplingContext(element, scvf);
        const auto& stokesVolVars = stokesElemVolVars[scvf.insideScvIdx()];
        const auto& darcyVolVars = stokesContext.volVars;
        const auto& outsideScv = (*scvs(stokesContext.fvGeometry).begin());

        const Scalar velocity = stokesElemFaceVars[scvf].velocitySelf();
        const bool insideIsUpstream = sign(velocity) == scvf.directionSign();

        return massFlux_(stokesIdx, darcyIdx, fvGeometry,
                         scvf, stokesVolVars, darcyVolVars,
                         outsideScv, velocity * scvf.directionSign(),
                         insideIsUpstream, diffCoeffAvgType);
    }

    /*!
     * \brief Returns the energy flux across the coupling boundary as seen from the Darcy domain.
     */
    template<bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar energyCouplingCondition(const Element<darcyIdx>& element,
                                   const FVElementGeometry<darcyIdx>& fvGeometry,
                                   const ElementVolumeVariables<darcyIdx>& darcyElemVolVars,
                                   const SubControlVolumeFace<darcyIdx>& scvf,
                                   const DiffusionCoefficientAveragingType diffCoeffAvgType = DiffusionCoefficientAveragingType::ffOnly) const
    {
        const auto& darcyContext = this->couplingManager().darcyCouplingContext(element, scvf);
        const auto& darcyVolVars = darcyElemVolVars[scvf.insideScvIdx()];
        const auto& stokesVolVars = darcyContext.volVars;

        const Scalar velocity = darcyContext.velocity * scvf.unitOuterNormal();
        const bool insideIsUpstream = velocity > 0.0;

        return energyFlux_(darcyIdx, stokesIdx, fvGeometry, darcyContext.fvGeometry, scvf,
                           darcyVolVars, stokesVolVars, velocity, insideIsUpstream, diffCoeffAvgType);
    }

    /*!
     * \brief Returns the energy flux across the coupling boundary as seen from the free-flow domain.
     */
    template<bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar energyCouplingCondition(const Element<stokesIdx>& element,
                                   const FVElementGeometry<stokesIdx>& fvGeometry,
                                   const ElementVolumeVariables<stokesIdx>& stokesElemVolVars,
                                   const ElementFaceVariables<stokesIdx>& stokesElemFaceVars,
                                   const SubControlVolumeFace<stokesIdx>& scvf,
                                   const DiffusionCoefficientAveragingType diffCoeffAvgType = DiffusionCoefficientAveragingType::ffOnly) const
    {
        const auto& stokesContext = this->couplingManager().stokesCouplingContext(element, scvf);
        const auto& stokesVolVars = stokesElemVolVars[scvf.insideScvIdx()];
        const auto& darcyVolVars = stokesContext.volVars;

        const Scalar velocity = stokesElemFaceVars[scvf].velocitySelf();
        const bool insideIsUpstream = sign(velocity) == scvf.directionSign();

        return energyFlux_(stokesIdx, darcyIdx, fvGeometry, stokesContext.fvGeometry, scvf,
                           stokesVolVars, darcyVolVars, velocity * scvf.directionSign(), insideIsUpstream, diffCoeffAvgType);
    }

protected:

    /*!
     * \brief Evaluate the compositional mole/mass flux across the interface.
     */
    template<std::size_t i, std::size_t j>
    NumEqVector massFlux_(Dune::index_constant<i> domainI,
                          Dune::index_constant<j> domainJ,
                          const FVElementGeometry<i>& insideFvGeometry,
                          const SubControlVolumeFace<i>& scvf,
                          const VolumeVariables<i>& insideVolVars,
                          const VolumeVariables<j>& outsideVolVars,
                          const SubControlVolume<j>& outsideScv,
                          const Scalar velocity,
                          const bool insideIsUpstream,
                          const DiffusionCoefficientAveragingType diffCoeffAvgType) const
    {
        NumEqVector flux(0.0);
        NumEqVector diffusiveFlux(0.0);

        auto moleOrMassFraction = [&](const auto& volVars, int phaseIdx, int compIdx)
        { return useMoles ? volVars.moleFraction(phaseIdx, compIdx) : volVars.massFraction(phaseIdx, compIdx); };

        auto moleOrMassDensity = [&](const auto& volVars, int phaseIdx)
        { return useMoles ? volVars.molarDensity(phaseIdx) : volVars.density(phaseIdx); };

        // treat the advective fluxes
        auto insideTerm = [&](int compIdx)
        { return moleOrMassFraction(insideVolVars, couplingPhaseIdx(domainI), compIdx) * moleOrMassDensity(insideVolVars, couplingPhaseIdx(domainI)); };

        auto outsideTerm = [&](int compIdx)
        { return moleOrMassFraction(outsideVolVars, couplingPhaseIdx(domainJ), compIdx) * moleOrMassDensity(outsideVolVars, couplingPhaseIdx(domainJ)); };


        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            const int domainICompIdx = couplingCompIdx(domainI, compIdx);
            const int domainJCompIdx = couplingCompIdx(domainJ, compIdx);
            flux[domainICompIdx] += this->advectiveFlux(insideTerm(domainICompIdx), outsideTerm(domainJCompIdx), velocity, insideIsUpstream);
        }

        // treat the diffusive fluxes
        const auto& insideScv = insideFvGeometry.scv(scvf.insideScvIdx());

        if (isFicksLaw)
            diffusiveFlux += diffusiveMolecularFluxFicksLaw_(domainI, domainJ, scvf, insideScv, outsideScv, insideVolVars, outsideVolVars, diffCoeffAvgType);
        else //maxwell stefan
            diffusiveFlux += diffusiveMolecularFluxMaxwellStefan_(domainI, domainJ, scvf, insideScv, outsideScv, insideVolVars, outsideVolVars);

        //convert to correct units if necessary
        if (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged && useMoles)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                const int domainICompIdx = couplingCompIdx(domainI, compIdx);
                diffusiveFlux[domainICompIdx] *= 1/FluidSystem<i>::molarMass(domainICompIdx);
            }
        }
        if (referenceSystemFormulation == ReferenceSystemFormulation::molarAveraged && !useMoles)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                const int domainICompIdx = couplingCompIdx(domainI, compIdx);
                diffusiveFlux[domainICompIdx] *= FluidSystem<i>::molarMass(domainICompIdx);
            }
        }

        flux += diffusiveFlux;
        // convert to total mass/mole balance, if set be user
        if (replaceCompEqIdx < numComponents)
            flux[replaceCompEqIdx] = std::accumulate(flux.begin(), flux.end(), 0.0);

        return flux;
    }

    Scalar getComponentEnthalpy(const VolumeVariables<stokesIdx>& volVars, int phaseIdx, int compIdx) const
    {
        return FluidSystem<stokesIdx>::componentEnthalpy(volVars.fluidState(), 0, compIdx);
    }

    Scalar getComponentEnthalpy(const VolumeVariables<darcyIdx>& volVars, int phaseIdx, int compIdx) const
    {
        return FluidSystem<darcyIdx>::componentEnthalpy(volVars.fluidState(), phaseIdx, compIdx);
    }

    /*!
     * \brief Evaluate the diffusive mole/mass flux across the interface.
     */
    template<std::size_t i, std::size_t j>
    NumEqVector diffusiveMolecularFluxMaxwellStefan_(Dune::index_constant<i> domainI,
                                                     Dune::index_constant<j> domainJ,
                                                     const SubControlVolumeFace<i>& scvfI,
                                                     const SubControlVolume<i>& scvI,
                                                     const SubControlVolume<j>& scvJ,
                                                     const VolumeVariables<i>& volVarsI,
                                                     const VolumeVariables<j>& volVarsJ) const
    {
        NumEqVector diffusiveFlux(0.0);

        const Scalar insideDistance = this->getDistance_(scvI, scvfI);
        const Scalar outsideDistance = this->getDistance_(scvJ, scvfI);

        ReducedComponentVector moleFracInside(0.0);
        ReducedComponentVector moleFracOutside(0.0);
        ReducedComponentVector reducedFlux(0.0);
        ReducedComponentMatrix reducedDiffusionMatrixInside(0.0);
        ReducedComponentMatrix reducedDiffusionMatrixOutside(0.0);

        //calculate the mole fraction vectors
        for (int compIdx = 0; compIdx < numComponents-1; compIdx++)
        {
            const int domainICompIdx = couplingCompIdx(domainI, compIdx);
            const int domainJCompIdx = couplingCompIdx(domainJ, compIdx);

            assert(FluidSystem<i>::componentName(domainICompIdx) == FluidSystem<j>::componentName(domainJCompIdx));

            //calculate x_inside
            const Scalar xInside = volVarsI.moleFraction(couplingPhaseIdx(domainI), domainICompIdx);
            //calculate outside molefraction with the respective transmissibility
            const Scalar xOutside = volVarsJ.moleFraction(couplingPhaseIdx(domainJ), domainJCompIdx);
            moleFracInside[domainICompIdx] = xInside;
            moleFracOutside[domainICompIdx] = xOutside;
        }

        //now we have to do the tpfa: J_i = -J_j which leads to: J_i = -rho_i Bi^-1 omegai(x*-xi) with x* = (omegai rho_i Bi^-1 + omegaj rho_j Bj^-1)^-1 (xi omegai rho_i Bi^-1 + xj omegaj rho_j Bj^-1) with i inside and j outside.

        //first set up the matrices containing the binary diffusion coefficients and mole fractions

        //inside matrix. KIdx and LIdx are the indices for the k and l-th component, N for the n-th component
        for (int compKIdx = 0; compKIdx < numComponents-1; compKIdx++)
        {
            const int domainICompKIdx = couplingCompIdx(domainI, compKIdx);
            const Scalar xk = volVarsI.moleFraction(couplingPhaseIdx(domainI), domainICompKIdx);
            const Scalar avgMolarMass = volVarsI.averageMolarMass(couplingPhaseIdx(domainI));
            const Scalar Mn = FluidSystem<i>::molarMass(numComponents-1);
            const Scalar tkn = volVarsI.effectiveDiffusionCoefficient(couplingPhaseIdx(domainI), domainICompKIdx, couplingCompIdx(domainI, numComponents-1));

            // set the entries of the diffusion matrix of the diagonal
            reducedDiffusionMatrixInside[domainICompKIdx][domainICompKIdx] += xk*avgMolarMass/(tkn*Mn);

            for (int compLIdx = 0; compLIdx < numComponents; compLIdx++)
            {
                const int domainICompLIdx = couplingCompIdx(domainI, compLIdx);

                // we don't want to calculate e.g. water in water diffusion
                if (domainICompKIdx == domainICompLIdx)
                    continue;

                const Scalar xl = volVarsI.moleFraction(couplingPhaseIdx(domainI), domainICompLIdx);
                const Scalar Mk = FluidSystem<i>::molarMass(domainICompKIdx);
                const Scalar Ml = FluidSystem<i>::molarMass(domainICompLIdx);
                const Scalar tkl = volVarsI.effectiveDiffusionCoefficient(couplingPhaseIdx(domainI), domainICompKIdx, domainICompLIdx);
                reducedDiffusionMatrixInside[domainICompKIdx][domainICompKIdx] += xl*avgMolarMass/(tkl*Mk);
                reducedDiffusionMatrixInside[domainICompKIdx][domainICompLIdx] += xk*(avgMolarMass/(tkn*Mn) - avgMolarMass/(tkl*Ml));
            }
        }
        //outside matrix
        for (int compKIdx = 0; compKIdx < numComponents-1; compKIdx++)
        {
            const int domainJCompKIdx = couplingCompIdx(domainJ, compKIdx);
            const int domainICompKIdx = couplingCompIdx(domainI, compKIdx);

            const Scalar xk = volVarsJ.moleFraction(couplingPhaseIdx(domainJ), domainJCompKIdx);
            const Scalar avgMolarMass = volVarsJ.averageMolarMass(couplingPhaseIdx(domainJ));
            const Scalar Mn = FluidSystem<j>::molarMass(numComponents-1);
            const Scalar tkn = volVarsJ.effectiveDiffusionCoefficient(couplingPhaseIdx(domainJ), domainJCompKIdx, couplingCompIdx(domainJ, numComponents-1));

            // set the entries of the diffusion matrix of the diagonal
            reducedDiffusionMatrixOutside[domainICompKIdx][domainICompKIdx] +=  xk*avgMolarMass/(tkn*Mn);

            for (int compLIdx = 0; compLIdx < numComponents; compLIdx++)
            {
                const int domainJCompLIdx = couplingCompIdx(domainJ, compLIdx);
                const int domainICompLIdx = couplingCompIdx(domainI, compLIdx);

                // we don't want to calculate e.g. water in water diffusion
                if (domainJCompLIdx == domainJCompKIdx)
                    continue;

                const Scalar xl = volVarsJ.moleFraction(couplingPhaseIdx(domainJ), domainJCompLIdx);
                const Scalar Mk = FluidSystem<j>::molarMass(domainJCompKIdx);
                const Scalar Ml = FluidSystem<j>::molarMass(domainJCompLIdx);
                const Scalar tkl = volVarsJ.effectiveDiffusionCoefficient(couplingPhaseIdx(domainJ), domainJCompKIdx, domainJCompLIdx);
                reducedDiffusionMatrixOutside[domainICompKIdx][domainICompKIdx] += xl*avgMolarMass/(tkl*Mk);
                reducedDiffusionMatrixOutside[domainICompKIdx][domainICompLIdx] += xk*(avgMolarMass/(tkn*Mn) - avgMolarMass/(tkl*Ml));
            }
        }

        const Scalar omegai = 1/insideDistance;
        const Scalar omegaj = 1/outsideDistance;

        reducedDiffusionMatrixInside.invert();
        reducedDiffusionMatrixInside *= omegai*volVarsI.density(couplingPhaseIdx(domainI));
        reducedDiffusionMatrixOutside.invert();
        reducedDiffusionMatrixOutside *= omegaj*volVarsJ.density(couplingPhaseIdx(domainJ));

        //in the helpervector we store the values for x*
        ReducedComponentVector helperVector(0.0);
        ReducedComponentVector gradientVectori(0.0);
        ReducedComponentVector gradientVectorj(0.0);

        reducedDiffusionMatrixInside.mv(moleFracInside, gradientVectori);
        reducedDiffusionMatrixOutside.mv(moleFracOutside, gradientVectorj);

        auto gradientVectorij = (gradientVectori + gradientVectorj);

        //add the two matrixes to each other
        reducedDiffusionMatrixOutside += reducedDiffusionMatrixInside;

        reducedDiffusionMatrixOutside.solve(helperVector, gradientVectorij);

        //Bi^-1 omegai rho_i (x*-xi). As we previously multiplied rho_i and omega_i wit the insidematrix, this does not need to be done again
        helperVector -=moleFracInside;
        reducedDiffusionMatrixInside.mv(helperVector, reducedFlux);

        reducedFlux *= -1;

        for (int compIdx = 0; compIdx < numComponents-1; compIdx++)
        {
            const int domainICompIdx = couplingCompIdx(domainI, compIdx);
            diffusiveFlux[domainICompIdx] = reducedFlux[domainICompIdx];
            diffusiveFlux[couplingCompIdx(domainI, numComponents-1)] -= reducedFlux[domainICompIdx];
        }
        return diffusiveFlux;
    }

    template<std::size_t i, std::size_t j>
    NumEqVector diffusiveMolecularFluxFicksLaw_(Dune::index_constant<i> domainI,
                                                Dune::index_constant<j> domainJ,
                                                const SubControlVolumeFace<i>& scvfI,
                                                const SubControlVolume<i>& scvI,
                                                const SubControlVolume<j>& scvJ,
                                                const VolumeVariables<i>& volVarsI,
                                                const VolumeVariables<j>& volVarsJ,
                                                const DiffusionCoefficientAveragingType diffCoeffAvgType) const
    {
        NumEqVector diffusiveFlux(0.0);

        const Scalar rhoInside = massOrMolarDensity(volVarsI, referenceSystemFormulation, couplingPhaseIdx(domainI));
        const Scalar rhoOutside = massOrMolarDensity(volVarsJ, referenceSystemFormulation, couplingPhaseIdx(domainJ));
        const Scalar avgDensity = 0.5 * rhoInside + 0.5 * rhoOutside;

        const Scalar insideDistance = this->getDistance_(scvI, scvfI);
        const Scalar outsideDistance = this->getDistance_(scvJ, scvfI);

        for (int compIdx = 1; compIdx < numComponents; ++compIdx)
        {
            const int domainIMainCompIdx = couplingPhaseIdx(domainI);
            const int domainJMainCompIdx = couplingPhaseIdx(domainJ);
            const int domainICompIdx = couplingCompIdx(domainI, compIdx);
            const int domainJCompIdx = couplingCompIdx(domainJ, compIdx);

            assert(FluidSystem<i>::componentName(domainICompIdx) == FluidSystem<j>::componentName(domainJCompIdx));

            const Scalar massOrMoleFractionInside = massOrMoleFraction(volVarsI, referenceSystemFormulation, couplingPhaseIdx(domainI), domainICompIdx);
            const Scalar massOrMoleFractionOutside = massOrMoleFraction(volVarsJ, referenceSystemFormulation, couplingPhaseIdx(domainJ), domainJCompIdx);

            const Scalar deltaMassOrMoleFrac = massOrMoleFractionOutside - massOrMoleFractionInside;
            const Scalar tij = this->transmissibility_(domainI,
                                                       domainJ,
                                                       insideDistance,
                                                       outsideDistance,
                                                       volVarsI.effectiveDiffusionCoefficient(couplingPhaseIdx(domainI), domainIMainCompIdx, domainICompIdx),
                                                       volVarsJ.effectiveDiffusionCoefficient(couplingPhaseIdx(domainJ), domainJMainCompIdx, domainJCompIdx),
                                                       diffCoeffAvgType);
            diffusiveFlux[domainICompIdx] += -avgDensity * tij * deltaMassOrMoleFrac;
        }

        const Scalar cumulativeFlux = std::accumulate(diffusiveFlux.begin(), diffusiveFlux.end(), 0.0);
        diffusiveFlux[couplingCompIdx(domainI, 0)] = -cumulativeFlux;

        return diffusiveFlux;
    }

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

        auto diffusiveFlux = isFicksLaw ? diffusiveMolecularFluxFicksLaw_(domainI, domainJ, scvf, insideScv, outsideScv, insideVolVars, outsideVolVars, diffCoeffAvgType)
                                        : diffusiveMolecularFluxMaxwellStefan_(domainI, domainJ, scvf, insideScv, outsideScv, insideVolVars, outsideVolVars);


        for (int compIdx = 0; compIdx < diffusiveFlux.size(); ++compIdx)
        {
            const int domainICompIdx = couplingCompIdx(domainI, compIdx);
            const int domainJCompIdx = couplingCompIdx(domainJ, compIdx);

            const bool insideDiffFluxIsUpstream = diffusiveFlux[domainICompIdx] > 0;
            const Scalar componentEnthalpy = insideDiffFluxIsUpstream ?
                                             getComponentEnthalpy(insideVolVars, couplingPhaseIdx(domainI), domainICompIdx)
                                           : getComponentEnthalpy(outsideVolVars, couplingPhaseIdx(domainJ), domainJCompIdx);

            if (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged)
                flux += diffusiveFlux[domainICompIdx] * componentEnthalpy;
            else
                flux += diffusiveFlux[domainICompIdx] * FluidSystem<i>::molarMass(domainICompIdx) * componentEnthalpy;
        }

        return flux;
    }
};

} // end namespace Dumux

#endif
