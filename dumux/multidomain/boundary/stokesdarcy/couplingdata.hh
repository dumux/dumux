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
     * \brief Defines which kind of averaging of diffusion coefficients
     *        (molecular diffusion or thermal conductance) at the interface
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
 * \tparam freeFlowIdx The domain index of the free-flow model.
 * \tparam porousMediumIdx The domain index of the porous-medium-flow model.
 * \tparam FFFS The free-flow fluidsystem.
 * \tparam hasAdapter Specifies whether an adapter class for the fluidsystem is used.
 */
template<std::size_t freeFlowIdx, std::size_t porousMediumIdx, class FFFS, bool hasAdapter>
struct IndexHelper;

/*!
 * \ingroup StokesDarcyCoupling
 * \brief Helper struct to choose the correct index for phases and components. This is need if the porous-medium-flow model
          features more fluid phases than the free-flow model. Specialization for the case that no adapter is used.
 * \tparam freeFlowIdx The domain index of the free-flow model.
 * \tparam porousMediumIdx The domain index of the porous-medium-flow model.
 * \tparam FFFS The free-flow fluidsystem.
 */
template<std::size_t freeFlowIdx, std::size_t porousMediumIdx, class FFFS>
struct IndexHelper<freeFlowIdx, porousMediumIdx, FFFS, false>
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
 * \tparam freeFlowIdx The domain index of the free-flow model.
 * \tparam porousMediumIdx The domain index of the porous-medium-flow model.
 * \tparam FFFS The free-flow fluidsystem.
 */
template<std::size_t freeFlowIdx, std::size_t porousMediumIdx, class FFFS>
struct IndexHelper<freeFlowIdx, porousMediumIdx, FFFS, true>
{
    /*!
     * \brief The free-flow model always uses phase index 0.
     */
    static constexpr auto couplingPhaseIdx(Dune::index_constant<freeFlowIdx>, int coupledPhaseIdx = 0)
    { return 0; }

    /*!
     * \brief The phase index of the porous-medium-flow model is given by the adapter fluidsytem (i.e., user input).
     */
    static constexpr auto couplingPhaseIdx(Dune::index_constant<porousMediumIdx>, int coupledPhaseIdx = 0)
    { return FFFS::multiphaseFluidsystemPhaseIdx; }

    /*!
     * \brief The free-flow model does not need any change of the component index.
     */
    static constexpr auto couplingCompIdx(Dune::index_constant<freeFlowIdx>, int coupledCompdIdx)
    { return coupledCompdIdx; }

    /*!
     * \brief The component index of the porous-medium-flow model is mapped by the adapter fluidsytem.
     */
    static constexpr auto couplingCompIdx(Dune::index_constant<porousMediumIdx>, int coupledCompdIdx)
    { return FFFS::compIdx(coupledCompdIdx); }
};

//! forward declare
template <class TypeTag, DiscretizationMethod discMethod>
class DarcysLawImplementation;

//! forward declare
template <class TypeTag, DiscretizationMethod discMethod>
class ForchheimersLawImplementation;


template<class MDTraits, class CouplingManager, bool enableEnergyBalance, bool isCompositional, DiscretizationMethod darcyDM>
class StokesDarcyCouplingDataImplementation;

/*!
* \ingroup BoundaryCoupling
* \brief Data for the coupling of a Darcy model (cell-centered finite volume)
*        with a (Navier-)Stokes model (staggerd grid).
*/
template<class MDTraits, class CouplingManager>
using StokesDarcyCouplingData = StokesDarcyCouplingDataImplementation<MDTraits, CouplingManager,
                                                                      GetPropType<typename MDTraits::template SubDomain<CouplingManager::freeFlowIdx>::TypeTag, Properties::ModelTraits>::enableEnergyBalance(),
                                                                      (GetPropType<typename MDTraits::template SubDomain<CouplingManager::freeFlowIdx>::TypeTag, Properties::ModelTraits>::numFluidComponents() > 1),
                                                                      GetPropType<typename MDTraits::template SubDomain<CouplingManager::porousMediumIdx>::TypeTag, Properties::GridGeometry>::discMethod>;

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

    static constexpr auto freeFlowIdx = CouplingManager::freeFlowIdx;
    static constexpr auto porousMediumIdx = CouplingManager::porousMediumIdx;

    static constexpr bool adapterUsed = ModelTraits<porousMediumIdx>::numFluidPhases() > 1;
    using IndexHelper = Dumux::IndexHelper<freeFlowIdx, porousMediumIdx, FluidSystem<freeFlowIdx>, adapterUsed>;

    static constexpr int enableEnergyBalance = GetPropType<SubDomainTypeTag<freeFlowIdx>, Properties::ModelTraits>::enableEnergyBalance();
    static_assert(GetPropType<SubDomainTypeTag<porousMediumIdx>, Properties::ModelTraits>::enableEnergyBalance() == enableEnergyBalance,
                  "All submodels must both be either isothermal or non-isothermal");

    static_assert(IsSameFluidSystem<FluidSystem<freeFlowIdx>,
                                    FluidSystem<porousMediumIdx>>::value,
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
    template<bool scalarPerm = std::is_same<typename Problem<porousMediumIdx>::SpatialParams::PermeabilityType, Scalar>::value,
             std::enable_if_t<scalarPerm, int> = 0>
    Scalar darcyPermeability(const Element<freeFlowIdx>& element, const SubControlVolumeFace<freeFlowIdx>& scvf) const
    {
        const auto& stokesContext = couplingManager().stokesCouplingContext(element, scvf);
        const auto perm = stokesContext.volVars.permeability();

        return perm;
    }

    /*!
     * \brief Returns the intrinsic permeability of the coupled Darcy element.
     */
    template<bool scalarPerm = std::is_same<typename Problem<porousMediumIdx>::SpatialParams::PermeabilityType, Scalar>::value,
             std::enable_if_t<!scalarPerm, int> = 0>
    Scalar darcyPermeability(const Element<freeFlowIdx>& element, const SubControlVolumeFace<freeFlowIdx>& scvf) const
    {
        const auto& stokesContext = couplingManager().stokesCouplingContext(element, scvf);
        const auto perm = stokesContext.volVars.permeability();
        const auto dirIdx = 1 - scvf.directionIndex();

        return perm[dirIdx][dirIdx];
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
            return domainI == freeFlowIdx
                            ? avgQuantityI / totalDistance
                            : avgQuantityJ / totalDistance;

        else // diffCoeffAvgType == DiffusionCoefficientAveragingType::pmOnly)
            return domainI == porousMediumIdx
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

private:
    const CouplingManager& couplingManager_;

};

} // end namespace Dumux

#include  <dumux/multidomain/boundary/stokesdarcy/cellcentered/tpfa/couplingdata.hh>
#include  <dumux/multidomain/boundary/stokesdarcy/box/couplingdata.hh>

#endif // DUMUX_STOKES_DARCY_COUPLINGDATA_HH
