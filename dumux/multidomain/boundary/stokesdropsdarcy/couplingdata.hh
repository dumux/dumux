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
 * \ingroup StokesDropsDarcyCoupling
 * \copydoc Dumux::StokesDropsDarcyCouplingData
 */

#ifndef DUMUX_STOKES_DROPS_DARCY_COUPLINGDATA_HH
#define DUMUX_STOKES_DROPS_DARCY_COUPLINGDATA_HH

#include <numeric>

#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>
#include <dumux/discretization/method.hh>
#include <dumux/multidomain/couplingmanager.hh>

namespace Dumux {

/*!
 * \ingroup StokesDarcyCoupling
 * \brief This structs holds a set of options which allow to modify the Stokes-Darcy
 *        coupling mechanism during runtime.
 */
//struct StokesDarcyCouplingOptions
//{
//    /*!
//     * \brief Defines which kind of averanging of diffusion coefficiencients
//     *        (moleculat diffusion or thermal conductance) at the interface
//     *        between free flow and porous medium shall be used.
//     */
//    enum class DiffusionCoefficientAveragingType
//    {
//        harmonic, arithmethic, ffOnly, pmOnly
//    };
//
//    /*!
//     * \brief Convenience function to convert user input given as std::string to the corresponding enum class used for chosing the type
//     *        of averaging of the diffusion/conduction parameter at the interface between the two domains.
//     */
//    static DiffusionCoefficientAveragingType stringToEnum(DiffusionCoefficientAveragingType, const std::string& diffusionCoefficientAveragingType)
//    {
//        if (diffusionCoefficientAveragingType == "Harmonic")
//            return DiffusionCoefficientAveragingType::harmonic;
//        else if (diffusionCoefficientAveragingType == "Arithmethic")
//            return DiffusionCoefficientAveragingType::arithmethic;
//        else if (diffusionCoefficientAveragingType == "FreeFlowOnly")
//            return DiffusionCoefficientAveragingType::ffOnly;
//        else if (diffusionCoefficientAveragingType == "PorousMediumOnly")
//            return DiffusionCoefficientAveragingType::pmOnly;
//        else
//            DUNE_THROW(Dune::IOError, "Unknown DiffusionCoefficientAveragingType");
//    }
//
//};

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

/*!
 * \ingroup StokesDarcyCoupling
 * \brief Helper struct to choose the correct index for phases and components. This is need if the porous-medium-flow model
          features more fluid phases than the free-flow model.
 * \tparam stokesIdx The domain index of the free-flow model.
 * \tparam interfaceIdx The domain index of the interface model.
 * \tparam darcyIdx The domain index of the porous-medium-flow model.
 * \tparam FFFS The free-flow fluidsystem.
 * \tparam hasAdapter Specifies whether an adapter class for the fluidsystem is used.
 */
template<std::size_t stokesIdx, std::size_t interfaceIdx, std::size_t darcyIdx, class FFFS, bool hasAdapter>
struct IndexHelper;

/*!
 * \ingroup StokesDarcyCoupling
 * \brief Helper struct to choose the correct index for phases and components. This is need if the porous-medium-flow model
          features more fluid phases than the free-flow model. Specialization for the case that no adapter is used.
 * \tparam stokesIdx The domain index of the free-flow model.
 * \tparam darcyIdx The domain index of the porous-medium-flow model.
 * \tparam FFFS The free-flow fluidsystem.
 */
template<std::size_t stokesIdx, std::size_t interfaceIdx, std::size_t darcyIdx, class FFFS>
struct IndexHelper<stokesIdx, interfaceIdx, darcyIdx, FFFS, false>
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
template<std::size_t stokesIdx, std::size_t interfaceIdx, std::size_t darcyIdx, class FFFS>
struct IndexHelper<stokesIdx, interfaceIdx, darcyIdx, FFFS, true>
{
    /*!
     * \brief The free-flow model always uses phase index 0.
     */
    static constexpr auto couplingPhaseIdx(Dune::index_constant<stokesIdx>, int coupledPhaseIdx = 0)
    { return 0; }

    /*!
     * \brief The phase index of the porous-medium-flow model is given by the adapter fluidsytem (i.e., user input).
     */
    static constexpr auto couplingPhaseIdx(Dune::index_constant<interfaceIdx>, int coupledPhaseIdx = 0)
    { return FFFS::multiphaseFluidsystemPhaseIdx; }

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

template<class MDTraits, class CouplingManager, bool enableEnergyBalance, bool isCompositional>
class StokesDropsDarcyCouplingDataImplementation;

/*!
* \ingroup BoundaryCoupling
* \brief Data for the coupling of a Darcy model (cell-centered finite volume)
*        with a (Navier-)Stokes model (staggerd grid).
*/
template<class MDTraits, class CouplingManager>
using StokesDropsDarcyCouplingData = StokesDropsDarcyCouplingDataImplementation<MDTraits, CouplingManager,
                                                                      GetPropType<typename MDTraits::template SubDomain<0>::TypeTag, Properties::ModelTraits>::enableEnergyBalance(),
                                                                      (GetPropType<typename MDTraits::template SubDomain<0>::TypeTag, Properties::ModelTraits>::numFluidComponents() > 1)>;

// TODO StokesDarcyCouplingOptions --> inherits from StokesDarcyCouplingData ?
// TODO Indexhelper --> inherits from StokesDarcyCouplingData ?
// TODO isSameFluidSystem --> inherits from StokesDarcyCouplingData ?
/*!
 * \ingroup StokesDarcyCoupling
 * \brief A base class which provides some common methods used for Stokes-Darcy coupling.
 */
template<class MDTraits, class CouplingManager>
class StokesDropsDarcyCouplingDataImplementationBase
{
    using Scalar = typename MDTraits::Scalar;

    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using FVGridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::FVGridGeometry>;
    template<std::size_t id> using Element = typename FVGridGeometry<id>::GridView::template Codim<0>::Entity;
    template<std::size_t id> using FVElementGeometry = typename FVGridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolumeFace = typename FVGridGeometry<id>::LocalView::SubControlVolumeFace;
//    template<std::size_t id> using SubControlVolume = typename FVGridGeometry<id>::LocalView::SubControlVolume;
    template<std::size_t id> using Indices = typename GetPropType<SubDomainTypeTag<id>, Properties::ModelTraits>::Indices;
    template<std::size_t id> using ElementVolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::LocalView;
//    template<std::size_t id> using VolumeVariables  = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::VolumeVariables;
//    template<std::size_t id> using Problem  = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using FluidSystem  = GetPropType<SubDomainTypeTag<id>, Properties::FluidSystem>;
    template<std::size_t id> using ModelTraits  = GetPropType<SubDomainTypeTag<id>, Properties::ModelTraits>;

    static constexpr auto stokesIdx = CouplingManager::stokesIdx;
    static constexpr auto interfaceIdx = CouplingManager::interfaceIdx;
    static constexpr auto darcyIdx = CouplingManager::darcyIdx;

    static constexpr bool adapterUsed = ModelTraits<darcyIdx>::numFluidPhases() > 1;
    using IndexHelper = Dumux::IndexHelper<stokesIdx, interfaceIdx, darcyIdx, FluidSystem<stokesIdx>, adapterUsed>;

    static constexpr int enableEnergyBalance = GetPropType<SubDomainTypeTag<stokesIdx>, Properties::ModelTraits>::enableEnergyBalance();
    static_assert(GetPropType<SubDomainTypeTag<darcyIdx>, Properties::ModelTraits>::enableEnergyBalance() == enableEnergyBalance,
                  "All submodels must both be either isothermal or non-isothermal");
    // TODO assert interface -> (non-)isothermal

    static_assert(IsSameFluidSystem<FluidSystem<stokesIdx>,
                                    FluidSystem<darcyIdx>>::value,
                  "All submodels must use the same fluid system");
    // TODO assert interface -> is same fluid system

//    using DiffusionCoefficientAveragingType = typename StokesDarcyCouplingOptions::DiffusionCoefficientAveragingType;

public:
    StokesDropsDarcyCouplingDataImplementationBase(const CouplingManager& couplingmanager): couplingManager_(couplingmanager) {}

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

    // TODO copied from stokesdarcycouplingmanager, inherit ?! --> interface permeability ???
    /*!
     * \brief Returns the intrinsic permeability of the coupled interface element.
     */
    Scalar interfacePermeability(const Element<stokesIdx>& element, const SubControlVolumeFace<stokesIdx>& scvf) const
    {
        const auto& stokesContext = couplingManager().stokesCouplingContext(element, scvf);
        return stokesContext.volVars.permeability();
    }

    // TODO copied from stokesdarcycouplingmanager, inherit ?!
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
        static constexpr auto numPhasesInterface = GetPropType<SubDomainTypeTag<interfaceIdx>, Properties::ModelTraits>::numFluidPhases();

        Scalar momentumFlux(0.0);
        const auto& stokesContext = couplingManager_.stokesCouplingContext(element, scvf);
        const auto interfacePhaseIdx = couplingPhaseIdx(interfaceIdx);

        // - p_pm * n_pm = p_pm * n_ff
        const Scalar interfacePressure = stokesContext.volVars.pressure(interfacePhaseIdx);

        if(numPhasesInterface > 1)
        {
            momentumFlux = interfacePressure;
//            std::cout << "** couplingdata: momentumCouplingCondition: p_if = " << momentumFlux << std::endl;
        }
        else // use pressure reconstruction for single phase models
        {
            // v = -K/mu * (gradP + rho*g)
            const Scalar velocity = stokesElemFaceVars[scvf].velocitySelf();
            const Scalar mu = stokesContext.volVars.viscosity(interfacePhaseIdx);
            const Scalar rho = stokesContext.volVars.density(interfacePhaseIdx);
            const Scalar distance = (stokesContext.element.geometry().center() - scvf.center()).two_norm();
            const Scalar g = -scvf.directionSign() * couplingManager_.problem(darcyIdx).gravity()[scvf.directionIndex()];
            const Scalar reconstructedInterfacePressure = ((scvf.directionSign() * velocity * (mu/interfacePermeability(element, scvf))) + rho * g) * distance + interfacePressure;
            momentumFlux = reconstructedInterfacePressure;
        }

        // normalize pressure
        if(getPropValue<SubDomainTypeTag<stokesIdx>, Properties::NormalizePressure>())
            momentumFlux -= couplingManager_.problem(stokesIdx).initial(scvf)[Indices<stokesIdx>::pressureIdx];

        momentumFlux *= scvf.directionSign();

        return momentumFlux;
    }

    // TODO copied from stokesdarcycouplingmanager, inherit ?!
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
//    /*!
//     * \brief Returns the distance between an scvf and the corresponding scv center.
//     */
//    template<class Scv, class Scvf>
//    Scalar getDistance_(const Scv& scv, const Scvf& scvf) const
//    {
//        return (scv.dofPosition() - scvf.ipGlobal()).two_norm();
//    }

private:
    const CouplingManager& couplingManager_;

};

/*!
 * \ingroup StokesDarcyCoupling
 * \brief Coupling data specialization for non-compositional models.
 */
template<class MDTraits, class CouplingManager, bool enableEnergyBalance>
class StokesDropsDarcyCouplingDataImplementation<MDTraits, CouplingManager, enableEnergyBalance, false>
: public StokesDropsDarcyCouplingDataImplementationBase<MDTraits, CouplingManager>
{
    using ParentType = StokesDropsDarcyCouplingDataImplementationBase<MDTraits, CouplingManager>;
    using Scalar = typename MDTraits::Scalar;
    static constexpr auto stokesIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto stokesCellCenterIdx = stokesIdx;
    static constexpr auto stokesFaceIdx = typename MDTraits::template SubDomain<1>::Index();
    static constexpr auto interfaceIdx = typename MDTraits::template SubDomain<2>::Index();
    static constexpr auto darcyIdx = typename MDTraits::template SubDomain<3>::Index();

    static constexpr auto interfaceNormal = 1.0; // interface parallel to x-direction, normal to y-direction

    // the sub domain type tags
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    template<std::size_t id> using FVGridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::FVGridGeometry>;
    template<std::size_t id> using Element = typename FVGridGeometry<id>::GridView::template Codim<0>::Entity;
    template<std::size_t id> using FVElementGeometry = typename FVGridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolumeFace = typename FVGridGeometry<id>::LocalView::SubControlVolumeFace;
    template<std::size_t id> using SubControlVolume = typename FVGridGeometry<id>::LocalView::SubControlVolume;
//    template<std::size_t id> using Indices = typename GetPropType<SubDomainTypeTag<id>, Properties::ModelTraits>::Indices;
    template<std::size_t id> using ElementVolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::LocalView;
    template<std::size_t id> using ElementFaceVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridFaceVariables>::LocalView;
//    template<std::size_t id> using VolumeVariables  = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::VolumeVariables;
    template<std::size_t id> using NumEqVector = GetPropType<SubDomainTypeTag<id>, Properties::NumEqVector>;
    template<std::size_t id> using LocalResidual = GetPropType<SubDomainTypeTag<id>, Properties::LocalResidual>;

    static_assert(GetPropType<SubDomainTypeTag<darcyIdx>, Properties::ModelTraits>::numFluidComponents() == GetPropType<SubDomainTypeTag<darcyIdx>, Properties::ModelTraits>::numFluidPhases(),
                  "Darcy Model must not be compositional");
    // TODO assert interface model is not compositional

//    using DiffusionCoefficientAveragingType = typename StokesDarcyCouplingOptions::DiffusionCoefficientAveragingType;

public:
    using ParentType::ParentType;
    using ParentType::couplingPhaseIdx;

    /*!
      * \name Coupling conditions between free flow / interface / porous medium
      */
     // \{

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
        const Scalar interfaceDensity = stokesContext.volVars.density(couplingPhaseIdx(interfaceIdx));
        const bool insideIsUpstream = sign(velocity) == scvf.directionSign();

//        return massFlux_(velocity * scvf.directionSign(), stokesDensity, interfaceDensity, insideIsUpstream);
        const auto massFlux = massFlux_(velocity * scvf.directionSign(), stokesDensity, interfaceDensity, insideIsUpstream);
        const Scalar usedVelocity = velocity * scvf.directionSign();
//        std::cout << "** couplingdata: massCouplingCondition for Stokes: massFlux = " << massFlux
//                  << ", v = " << usedVelocity << ", rho_ff = " << stokesDensity << ", rho_if = " << interfaceDensity
//                  << std::endl;
        return massFlux;
    }

    /*!
     * \brief Returns the mass flux across the coupling boundary as seen from the interface domain.
     */
    Scalar massCouplingCondition(const Element<interfaceIdx>& element,
                                 const FVElementGeometry<interfaceIdx>& fvGeometry,
                                 const ElementVolumeVariables<interfaceIdx>& interfaceElemVolVars,
                                 const SubControlVolume<interfaceIdx>& scv) const
    {
        const auto& interfaceContext = this->couplingManager().interfaceCouplingContext(element, scv);

        // flux from free flow
        const Scalar topVelocity = -1.0 * interfaceContext.stokesVelocity[interfaceNormal];
        const Scalar interfaceDensity = interfaceElemVolVars[scv.dofIndex()].density(couplingPhaseIdx(interfaceIdx));
        const Scalar stokesDensity = interfaceContext.stokesVolVars.density();
        const bool topInsideIsUpstream = topVelocity > 0.0;

        auto massFlux = -1.0 * massFlux_(topVelocity, interfaceDensity, stokesDensity, topInsideIsUpstream);
//        std::cout << "** couplingdata: massCouplingCondition for Interface: massFluxFF = " << massFlux
//                  << ", v = " << topVelocity << std::endl;

        // compute velocity at the interface (pressure gradient between interface and porous medium)
        const Scalar diffP = interfaceElemVolVars[scv.dofIndex()].pressure(couplingPhaseIdx(interfaceIdx)) - interfaceContext.darcyVolVars.pressure(couplingPhaseIdx(darcyIdx));
        // TODO sign of diffP determines volvars = if-volvars / pm-volvars !!
        const Scalar diffY = scv.center()[interfaceNormal] - interfaceContext.darcyElementCenter[interfaceNormal];
        const Scalar darcyDensity = interfaceContext.darcyVolVars.density(couplingPhaseIdx(darcyIdx));
        const Scalar gravity = -1.0 * this->couplingManager().problem(interfaceIdx).gravity()[interfaceNormal];
        const Scalar permeability = interfaceContext.darcyVolVars.permeability();
        const Scalar viscosity = interfaceContext.darcyVolVars.viscosity(couplingPhaseIdx(darcyIdx));
        const Scalar relPermeability = 1.0; // TODO
        const Scalar bottomVelocity = -1.0 * permeability * relPermeability / viscosity * (diffP/diffY - gravity * darcyDensity);

        // TODO flux from porous medium
        const bool bottomInsideIsUpstream = bottomVelocity > 0.0;
        auto massFluxPM = -1.0 * massFlux_(bottomVelocity, interfaceDensity, darcyDensity, bottomInsideIsUpstream);
//        std::cout << "** couplingdata: massCouplingCondition for Interface: massFluxPM = " << massFluxPM
//                  << ", v = " << bottomVelocity << std::endl;
        massFlux += massFluxPM;

        return massFlux;
    }

    /*!
     * \brief Returns the mass flux across the coupling boundary as seen from the interface domain.
     */
    Scalar massCouplingCondition(const Element<darcyIdx>& element,
                                 const FVElementGeometry<darcyIdx>& fvGeometry,
                                 const ElementVolumeVariables<darcyIdx>& darcyElemVolVars,
                                 const SubControlVolumeFace<darcyIdx>& scvf) const
    {
        const auto& darcyContext = this->couplingManager().darcyCouplingContext(element, scvf);
        // compute velocity at the interface (pressure gradient between interface and porous medium)
        const Scalar diffP = darcyContext.volVars.pressure(couplingPhaseIdx(interfaceIdx)) - darcyElemVolVars[scvf.insideScvIdx()].pressure(couplingPhaseIdx(darcyIdx));
        // TODO sign of diffP determines volvars = if-volvars / pm-volvars !!
        const auto& scv = (*scvs(fvGeometry).begin());
        const Scalar diffY = darcyContext.elementCenter - scv.center()[interfaceNormal];
        const Scalar darcyDensity = darcyElemVolVars[scvf.insideScvIdx()].density(couplingPhaseIdx(darcyIdx));
        const Scalar gravity = -1.0 * this->couplingManager().problem(darcyIdx).gravity()[interfaceNormal];
        const Scalar permeability = darcyElemVolVars[scvf.insideScvIdx()].permeability();
        const Scalar viscosity = darcyElemVolVars[scvf.insideScvIdx()].viscosity(couplingPhaseIdx(darcyIdx));
        const Scalar relPermeability = 1.0; // TODO
        const Scalar velocity = - 1.0 * permeability * relPermeability / viscosity * (diffP/diffY - gravity * darcyDensity); // TODO use Darcy's law?
        const Scalar interfaceDensity = darcyContext.volVars.density(couplingPhaseIdx(interfaceIdx));
        const bool insideIsUpstream = velocity > 0.0;

        return massFlux_(velocity, darcyDensity, interfaceDensity, insideIsUpstream);
//        const auto massFlux = massFlux_(velocity, darcyDensity, interfaceDensity, insideIsUpstream);
//        std::cout << "** couplingdata: massCouplingCondition for Darcy: massFlux = " << massFlux
//                  << ", v = " << velocity << std::endl;
//        return massFlux;
    }
    // \}

private:
    /*!
     * \brief Evaluate the mole/mass flux between Stokes and interface domain.
     */
    // TODO ParentType::massFlux_?
    Scalar massFlux_(const Scalar velocity,
                     const Scalar insideDensity,
                     const Scalar outSideDensity,
                     bool insideIsUpstream) const
    {
        return this->advectiveFlux(insideDensity, outSideDensity, velocity, insideIsUpstream);
    }

};

// TODO add specialization for compositional flow (see StokesDarcyCouplingData)

} // end namespace Dumux

#endif // DUMUX_STOKES_DROPS_DARCY_COUPLINGDATA_HH
