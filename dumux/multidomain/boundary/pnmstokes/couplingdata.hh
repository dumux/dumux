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
 * \ingroup BoubdaryCoupling
 * \copydoc Dumux::PNMStokesCouplingData
 */

#ifndef DUMUX_PNM_STOKES_COUPLINGDATA_HH
#define DUMUX_PNM_STOKES_COUPLINGDATA_HH

#include <type_traits>
#include <dumux/common/properties.hh>
#include <dumux/multidomain/couplingmanager.hh>
#include <dumux/freeflow/navierstokes/staggered/velocitygradients.hh>
#include <dumux/flux/referencesystemformulation.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/typetraits/problem.hh>
#include "geometry.hh"

namespace Dumux {

// forward declaration
namespace FluidSystems {
template <class MPFluidSystem, int phase>
class OnePAdapter;
}

template<class T>
struct UsesAdapter : std::false_type {};

template<class MPFluidSystem, int phase>
struct UsesAdapter<FluidSystems::OnePAdapter<MPFluidSystem, phase>> : std::true_type {};

/*!
 * \ingroup MultiDomain
 * \ingroup BoundaryCoupling
 * \ingroup StokesDarcyCoupling
 * \brief This structs helps to check if the two sub models use the same fluidsystem.
 *        Specialization for the case of using an adapter only for the free-flow model.
 * \tparam FFFS The free-flow fluidsystem
 * \tparam PMFS The porous-medium flow fluidsystem
 */
template<class FFFS, class PMFS>
struct IsSameFluidSystem
{
    static_assert(UsesAdapter<FFFS>(), "Free flow fluid system must use an adapter!");
    static_assert(FFFS::numPhases == 1, "Only single-phase fluidsystems may be used for free flow.");
    static constexpr bool value = std::is_same<typename FFFS::MultiPhaseFluidSystem, PMFS>::value;
};

/*!
 * \ingroup MultiDomain
 * \ingroup BoundaryCoupling
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
 * \ingroup MultiDomain
 * \ingroup BoundaryCoupling
 * \ingroup StokesDarcyCoupling
 * \brief Helper struct to choose the correct index for phases and components. This is need if the porous-medium-flow model
          features more fluid phases than the free-flow model.
 * \tparam stokesIdx The domain index of the free-flow model.
 * \tparam darcyIdx The domain index of the porous-medium-flow model.
 * \tparam FFFS The free-flow fluidsystem.
 * \tparam indexRequiresMapping Specifies whether the phase or component indices need to be mapped (e.g, in an adapt fluidsystem is used)
 */
template<std::size_t stokesIdx, std::size_t darcyIdx, class FFFS, bool indexRequiresMapping>
struct IndexHelper;

/*!
 * \ingroup MultiDomain
 * \ingroup BoundaryCoupling
 * \ingroup StokesDarcyCoupling
 * \brief Helper struct to choose the correct index for phases and components. This is needed if the porous-medium-flow model
          features more fluid phases than the free-flow model. Specialization for the case that no adapter is used or both models use the adapter.
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
 * \ingroup MultiDomain
 * \ingroup BoundaryCoupling
 * \ingroup StokesDarcyCoupling
 * \brief Helper struct to choose the correct index for phases and components. This is need if the porous-medium-flow model
          features more fluid phases than the free-flow model. Specialization for the case that only one model uses an adapter.
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

template<class MDTraits, class CouplingManager, bool enableEnergyBalance, bool isCompositional>
class PNMStokesCouplingDataImplementation;

/*!
* \ingroup MultiDomain
* \ingroup BoundaryCoupling
* \brief Data for the coupling of a Darcy model (cell-centered finite volume)
*        with a (Navier-)Stokes model (staggerd grid).
*/
template<class MDTraits, class CouplingManager>
using PNMStokesCouplingData = PNMStokesCouplingDataImplementation<MDTraits, CouplingManager,
                                                                  GetPropType<typename MDTraits::template SubDomain<0>::TypeTag, Properties::ModelTraits>::enableEnergyBalance(),
                                                                  (GetPropType<typename MDTraits::template SubDomain<0>::TypeTag, Properties::ModelTraits>::numFluidComponents() > 1)>;

/*!
 * \ingroup MixedDimension
 * \ingroup MixedDimensionBoundary
 * \brief Data for the coupling of a pore-network model with a staggerd grid Navier-Stokes model
 */
template<class MDTraits, class CouplingManager>
class PNMStokesCouplingDataImplementationBase
{
    using Scalar = typename MDTraits::Scalar;
    static constexpr auto bulkIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto lowDimIdx = typename MDTraits::template SubDomain<2>::Index();

    // the sub domain type tags
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using NumEqVector = Dumux::NumEqVector<GetPropType<SubDomainTypeTag<id>, Properties::PrimaryVariables>>;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolume = typename FVElementGeometry<id>::SubControlVolume;
    template<std::size_t id> using SubControlVolumeFace = typename FVElementGeometry<id>::SubControlVolumeFace;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using VolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVariables>::GridVolumeVariables::VolumeVariables;
    template<std::size_t id> using FluidSystem = typename VolumeVariables<id>::FluidSystem;
    template<std::size_t id> using ModelTraits = GetPropType<SubDomainTypeTag<id>, Properties::ModelTraits>;
    template<std::size_t id> using Indices = typename ModelTraits<id>::Indices;
    template<std::size_t id> using BoundaryTypes = typename ProblemTraits<Problem<id>>::BoundaryTypes;

    using VelocityVector = Dune::FieldVector<Scalar, GridView<bulkIdx>::dimension>;
    using VelocityGradients = StaggeredVelocityGradients<Scalar, GridGeometry<bulkIdx>, BoundaryTypes<bulkIdx>, Indices<bulkIdx>>;

    static constexpr bool indexRequiresMapping = UsesAdapter<FluidSystem<bulkIdx>>() != UsesAdapter<FluidSystem<lowDimIdx>>();
    using IndexHelper = Dumux::IndexHelper<bulkIdx, lowDimIdx, FluidSystem<bulkIdx>, indexRequiresMapping>;
    static constexpr bool enableEnergyBalance = ModelTraits<bulkIdx>::enableEnergyBalance();
    static_assert(ModelTraits<lowDimIdx>::enableEnergyBalance() == enableEnergyBalance,
                  "All submodels must both be either isothermal or non-isothermal");

    static_assert(IsSameFluidSystem<FluidSystem<bulkIdx>,
                                    FluidSystem<lowDimIdx>
                                    >::value,
                  "All submodels must use the same fluid system");

public:

    PNMStokesCouplingDataImplementationBase(const CouplingManager& couplingmanager)
    : couplingManager_(couplingmanager)
    , lowDimProblem_(&couplingmanager.problem(lowDimIdx)) {}

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
     * \brief Evaluate an advective flux across the interface and consider upwinding.
     */
    Scalar advectiveFlux(const Scalar insideQuantity, const Scalar outsideQuantity, const Scalar volumeFlow, bool insideIsUpstream) const
    {
        const Scalar upwindWeight = 1.0; //TODO use Implicit.UpwindWeight or something like Coupling.UpwindWeight?

        if (insideIsUpstream)
            return (upwindWeight * insideQuantity + (1.0 - upwindWeight) * outsideQuantity) * volumeFlow;
        else
            return (upwindWeight * outsideQuantity + (1.0 - upwindWeight) * insideQuantity) * volumeFlow;
    }


    //! Compute the velocity and orientation of the fluid leaving a pore throat at the coupling interface
    VelocityVector boundaryVelocity(const Element<bulkIdx>& element,
                                    const SubControlVolumeFace<bulkIdx>& scvf,
                                    const bool verbose = false) const
    {
        assert(couplingManager_.bulkCouplingContext(element, scvf).size() == 1);
        const auto& context = couplingManager_.bulkCouplingContext(element, scvf)[0];

        const auto& lowDimElement = context.element;
        const auto& lowDimFvGeometry = context.fvGeometry;
        const auto& lowDimScvf = lowDimFvGeometry.scvf(0);
        const auto& lowDimElemVolVars = context.elemVolVars;
        const auto lowDimPhaseIdx = couplingPhaseIdx(lowDimIdx);

        const Scalar area = context.elemFluxVarsCache[lowDimScvf].throatCrossSectionalArea(lowDimPhaseIdx);

        // only proceed if area > 0 in order to prevent division by zero (e.g., when the throat was not invaded yet)
        if (area > 0.0)
        {
            using LowDimFluxVariables = GetPropType<SubDomainTypeTag<lowDimIdx>, Properties::FluxVariables>;
            LowDimFluxVariables fluxVars;
            fluxVars.init(*lowDimProblem_, lowDimElement, lowDimFvGeometry, lowDimElemVolVars, lowDimScvf, context.elemFluxVarsCache);

            const Scalar flux = fluxVars.advectiveFlux(lowDimPhaseIdx, [lowDimPhaseIdx](const auto& volVars){ return volVars.mobility(lowDimPhaseIdx);});

            // account for the orientation of the bulk face.
            VelocityVector velocity = (lowDimElement.geometry().corner(1) - lowDimElement.geometry().corner(0));
            velocity /= velocity.two_norm();
            velocity *= flux / area;

            // TODO: Multiple throats connected to the same pore?
            return velocity;
        }
        else
            return VelocityVector(0.0);
    }

    //! Evaluate the pressure value in a free flow grid cell at the coupling interface TODO: reconstruct value directly at interface
    auto bulkPrivar(const Element<lowDimIdx>& element, const SubControlVolume<lowDimIdx>& scv) const
    {
        Scalar pressure = 0.0;
        const auto& context = couplingManager_.lowDimCouplingContext(element, scv);

        for (const auto& i : context)
            pressure += i.volVars.pressure();

        pressure /= context.size();

        return pressure;
    }

    Scalar coupledRadius(const Element<bulkIdx>& element,
                         const SubControlVolumeFace<bulkIdx>& scvf) const
    {
        static const bool coupleOverPoreRadius = getParamFromGroup<bool>(couplingManager_.problem(CouplingManager::bulkIdx).paramGroup(), "Grid.CoupleOverPoreRadius", false);
        using GlobalPosition = typename Element<lowDimIdx>::Geometry::GlobalCoordinate;
        static const auto couplingPlaneNormal = getParamFromGroup<GlobalPosition>(couplingManager_.problem(CouplingManager::bulkIdx).paramGroup(),
                                                                                  "Grid.CouplingPlaneNormal",
                                                                                  [](){ GlobalPosition tmp(0.0); tmp[tmp.size()-1] = 1.0; return tmp; }());

        const auto& stokesContext = couplingManager_.bulkCouplingContext(element, scvf)[0];
        for (auto&& scv : scvs(stokesContext.fvGeometry))
        {
            if (couplingManager_.isCoupledDof(lowDimIdx, scv.dofIndex()))
            {
                if (coupleOverPoreRadius)
                    return stokesContext.elemVolVars[scv].poreInscribedRadius();
                else
                {
                    return projectedThroatRadius(lowDimProblem_->spatialParams().throatInscribedRadius(stokesContext.element, stokesContext.elemVolVars),
                                                 stokesContext.element, couplingPlaneNormal);
                }
            }
        }
        DUNE_THROW(Dune::InvalidStateException, "No coupled pore found");
    }



    template<class ElementVolumeVariables, class ElementFaceVariables>
    Scalar momentumCouplingCondition(const Element<bulkIdx>& element,
                                     const FVElementGeometry<bulkIdx>& fvGeometry,
                                     const ElementVolumeVariables& stokesElemVolVars,
                                     const ElementFaceVariables& stokesElemFaceVars,
                                     const SubControlVolumeFace<bulkIdx>& scvf) const
    {
        // set p_freeflow = p_PNM
        Scalar momentumFlux(0.0);
        const auto& stokesContext = couplingManager_.bulkCouplingContext(element, scvf)[0];
        const auto lowDimPhaseIdx = couplingPhaseIdx(lowDimIdx);
        Scalar pnmPressure = 0;

        for (auto&& scv : scvs(stokesContext.fvGeometry))
        {
            if (couplingManager_.isCoupledDof(lowDimIdx, scv.dofIndex()))
                pnmPressure = stokesContext.elemVolVars[scv].pressure(lowDimPhaseIdx);
            // TODO: reconstruct pressure for inclined throats
        }

        momentumFlux = pnmPressure;

        // normalize pressure
        if (getPropValue<SubDomainTypeTag<bulkIdx>, Properties::NormalizePressure>())
            momentumFlux -= couplingManager_.problem(bulkIdx).initial(scvf)[Indices<bulkIdx>::pressureIdx];

        // Explicitly account for dv_i/dx_i, which is NOT part of the actual coupling condition. We do it here for convenience so
        // we do not forget to set it in the problem. We assume that the velocity gradient at the boundary towards the interface is the same
        // as the one in the center of the element.
        momentumFlux += VelocityGradients::velocityGradII(scvf, stokesElemFaceVars[scvf]) * stokesElemVolVars[scvf.insideScvIdx()].effectiveViscosity();

        // We do NOT consider the intertia term here. If included, Newton convergence decreases drastically and the solution even does not converge to a reference solution.
        // We furthermore assume creeping flow within the boundary layer thus neglecting this term is physically justified.

        momentumFlux *= scvf.directionSign();

        return momentumFlux;
    }

    /*!
     * \brief Returns a reference to the coupling manager.
     */
    const CouplingManager& couplingManager() const
    { return couplingManager_; }

    // }
protected:

    /*!
     * \brief Evaluate the diffusive mole/mass flux across the interface.
     */
    template<std::size_t i, std::size_t j, bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar conductiveEnergyFlux_(Dune::index_constant<i> domainI,
                                 Dune::index_constant<j> domainJ,
                                 const SubControlVolumeFace<bulkIdx>& scvf,
                                 const SubControlVolume<i>& scvI,
                                 const SubControlVolume<j>& scvJ,
                                 const VolumeVariables<i>& volVarsI,
                                 const VolumeVariables<j>& volVarsJ) const
    {
        const auto& bulkVolVars = this->getBulkVolVars_(volVarsI, volVarsJ);

        const Scalar distance = this->getDistance_(scvI, scvJ, scvf);

        const Scalar deltaT = volVarsJ.temperature() - volVarsI.temperature();
        const Scalar tij = bulkVolVars.thermalConductivity() / distance;

        return -deltaT * tij;
    }

    template<class ScvI,class ScvJ>
    Scalar getDistance_(const ScvI& scvI, const ScvJ& scvJ, const SubControlVolumeFace<bulkIdx>& scvf) const
    {
        if (std::is_same<ScvI, SubControlVolume<bulkIdx>>::value)
            return  (scvI.center() - scvf.center()).two_norm();
        else if (std::is_same<ScvJ, SubControlVolume<bulkIdx>>::value)
            return (scvJ.center() - scvf.center()).two_norm();
        else
            DUNE_THROW(Dune::InvalidStateException, "None of the scvs are bulk scvs");
    }

    const VolumeVariables<bulkIdx>& getBulkVolVars_(const VolumeVariables<bulkIdx>& bulkVolVars, const VolumeVariables<lowDimIdx>&) const
    {
        return bulkVolVars;
    }

    const VolumeVariables<bulkIdx>& getBulkVolVars_(const VolumeVariables<lowDimIdx>&, const VolumeVariables<bulkIdx>& bulkVolVars) const
    {
        return bulkVolVars;
    }

    const CouplingManager& couplingManager_;
    const Problem<lowDimIdx>* lowDimProblem_;
};

/*!
 * \ingroup MultiDomain
 * \ingroup BoundaryCoupling
 * \ingroup StokesDarcyCoupling
 * \brief Coupling data specialization for non-compositional models.
 */
template<class MDTraits, class CouplingManager, bool enableEnergyBalance>
class PNMStokesCouplingDataImplementation<MDTraits, CouplingManager, enableEnergyBalance, false>
: public PNMStokesCouplingDataImplementationBase<MDTraits, CouplingManager>
{
    static constexpr auto bulkIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto lowDimIdx = typename MDTraits::template SubDomain<2>::Index();

    // the sub domain type tags
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    using ParentType = PNMStokesCouplingDataImplementationBase<MDTraits, CouplingManager>;
    using Scalar = typename MDTraits::Scalar;

    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using Element = typename GridGeometry<id>::GridView::template Codim<0>::Entity;
    template<std::size_t id> using SubControlVolume = typename FVElementGeometry<id>::SubControlVolume;
    template<std::size_t id> using SubControlVolumeFace = typename FVElementGeometry<id>::SubControlVolumeFace;

    template<std::size_t id> using ElementVolumeVariables  = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::LocalView;
    template<std::size_t id> using VolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::VolumeVariables;

public:

    using ParentType::ParentType;
    using ParentType::couplingPhaseIdx;

    /*!
     * \brief Returns the mass flux across the coupling boundary as seen from the pore-network domain.
     */
    Scalar massCouplingCondition(const Element<lowDimIdx>& element,
                                 const FVElementGeometry<lowDimIdx>& fvGeometry,
                                 const ElementVolumeVariables<lowDimIdx>& elemVolVars,
                                 const SubControlVolume<lowDimIdx>& scv) const
    {
        const auto& lowDimContext = this->couplingManager().lowDimCouplingContext(element, scv);
        Scalar massFlux(0.0);

        for (const auto& i : lowDimContext)
        {
            const auto& velocity = i.velocity;
            const Scalar lowDimDensity = elemVolVars[scv].density(couplingPhaseIdx(lowDimIdx));
            const Scalar bulkDensity = i.volVars.density(couplingPhaseIdx(bulkIdx));
            const auto& bulkScvf = i.getBulkScvf();
            const Scalar area = bulkScvf.area() * i.volVars.extrusionFactor();

            const bool lowDimIsUpstream = sign(velocity) != bulkScvf.directionSign();

            // This value is used as a source term which implies a sign flip. When used as a Neumann flux, the correct sign would be -bulkScvf.directionSign().
            massFlux += this->advectiveFlux(lowDimDensity, bulkDensity, velocity, lowDimIsUpstream) * area * bulkScvf.directionSign();
        }

        return massFlux;
    }


    /*!
     * \brief Returns the mass flux across the coupling boundary as seen from the free-flow domain.
     */
    template<class ElementFaceVariables>
    Scalar massCouplingCondition(const Element<bulkIdx>& element,
                                 const FVElementGeometry<bulkIdx>& fvGeometry,
                                 const ElementVolumeVariables<bulkIdx>& elemVolVars,
                                 const ElementFaceVariables& elemFaceVars,
                                 const SubControlVolumeFace<bulkIdx>& scvf) const
    {
        const auto& stokesContext = this->couplingManager().bulkCouplingContext(element, scvf)[0];
        const Scalar velocity = elemFaceVars[scvf].velocitySelf();
        const Scalar stokesDensity = elemVolVars[scvf.insideScvIdx()].density();

        for (auto&& scv : scvs(stokesContext.fvGeometry))
        {
            if (this->couplingManager().isCoupledDof(lowDimIdx, scv.dofIndex()))
            {
                const Scalar pnmDensity = stokesContext.elemVolVars[scv].density(couplingPhaseIdx(lowDimIdx));
                const bool insideIsUpstream = sign(velocity) == scvf.directionSign();
                return this->advectiveFlux(stokesDensity, pnmDensity, velocity, insideIsUpstream) * scvf.directionSign();
            }
        }
        DUNE_THROW(Dune::InvalidStateException, "No coupled scvf found");
    }
};

/*!
 * \ingroup MultiDomain
 * \ingroup BoundaryCoupling
 * \ingroup StokesDarcyCoupling
 * \brief Coupling data specialization for compositional models.
 */
template<class MDTraits, class CouplingManager, bool enableEnergyBalance>
class PNMStokesCouplingDataImplementation<MDTraits, CouplingManager, enableEnergyBalance, true>
: public PNMStokesCouplingDataImplementationBase<MDTraits, CouplingManager>
{
    static constexpr auto bulkIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto lowDimIdx = typename MDTraits::template SubDomain<2>::Index();

    // the sub domain type tags
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    using ParentType = PNMStokesCouplingDataImplementationBase<MDTraits, CouplingManager>;
    using Scalar = typename MDTraits::Scalar;

    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using Element = typename GridGeometry<id>::GridView::template Codim<0>::Entity;
    template<std::size_t id> using SubControlVolume = typename FVElementGeometry<id>::SubControlVolume;
    template<std::size_t id> using SubControlVolumeFace = typename FVElementGeometry<id>::SubControlVolumeFace;

    template<std::size_t id> using ElementVolumeVariables  = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::LocalView;
    template<std::size_t id> using VolumeVariables  = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::VolumeVariables;
    template<std::size_t id> using FluidSystem = typename VolumeVariables<id>::FluidSystem;

    static constexpr auto replaceCompEqIdx = GetPropType<SubDomainTypeTag<bulkIdx>, Properties::ModelTraits>::replaceCompEqIdx();
    static constexpr bool useMoles = GetPropType<SubDomainTypeTag<bulkIdx>, Properties::ModelTraits>::useMoles();
    static constexpr auto numComponents = GetPropType<SubDomainTypeTag<bulkIdx>, Properties::ModelTraits>::numFluidComponents();
    static constexpr auto referenceSystemFormulation = GetPropType<SubDomainTypeTag<bulkIdx>, Properties::MolecularDiffusionType>::referenceSystemFormulation();

    static_assert(GetPropType<SubDomainTypeTag<lowDimIdx>, Properties::ModelTraits>::numFluidComponents() == numComponents, "Submodels must use same number of components");
    static_assert(GetPropType<SubDomainTypeTag<lowDimIdx>, Properties::ModelTraits>::useMoles() == useMoles, "Both models must either use moles or not");
    static_assert(GetPropType<SubDomainTypeTag<lowDimIdx>, Properties::ModelTraits>::replaceCompEqIdx() == replaceCompEqIdx, "Both models must use the same replaceCompEqIdx");
    static_assert(GetPropType<SubDomainTypeTag<lowDimIdx>, Properties::MolecularDiffusionType>::referenceSystemFormulation() == referenceSystemFormulation,
                  "Both models must use the same reference system formulation for diffusion");

    using NumEqVector = Dune::FieldVector<Scalar, numComponents>;

public:

    using ParentType::ParentType;
    using ParentType::couplingPhaseIdx;
    using ParentType::couplingCompIdx;

    /*!
     * \brief Returns the mass flux across the coupling boundary as seen from the pore-network domain.
     */
    NumEqVector massCouplingCondition(const Element<lowDimIdx>& element,
                                      const FVElementGeometry<lowDimIdx>& fvGeometry,
                                      const ElementVolumeVariables<lowDimIdx>& elemVolVars,
                                      const SubControlVolume<lowDimIdx>& scv) const
    {
        const auto& lowDimContext = this->couplingManager().lowDimCouplingContext(element, scv);
        NumEqVector massFlux(0.0);

        for (const auto& i : lowDimContext)
        {
            const Scalar velocity = i.velocity;
            const auto& bulkScvf = i.getBulkScvf();
            const Scalar area = bulkScvf.area();

            const bool insideIsUpstream = sign(velocity) != bulkScvf.directionSign();

            auto flux = massFlux_(lowDimIdx, bulkIdx, bulkScvf, scv, i.fvGeometry.scv(bulkScvf.insideScvIdx()), elemVolVars[scv], i.volVars, velocity, insideIsUpstream);
            flux *= area * i.volVars.extrusionFactor();
            flux *= -1.0; // flip the sign because the flux is used as a source term (different sign convention)

            massFlux += flux;
        }

        return massFlux;
    }

    /*!
     * \brief Returns the mass flux across the coupling boundary as seen from the free-flow domain.
     */
    template<class ElementFaceVariables>
    NumEqVector massCouplingCondition(const Element<bulkIdx>& element,
                                      const FVElementGeometry<bulkIdx>& fvGeometry,
                                      const ElementVolumeVariables<bulkIdx>& elemVolVars,
                                      const ElementFaceVariables& elemFaceVars,
                                      const SubControlVolumeFace<bulkIdx>& scvf) const
    {
        const auto& stokesContext = this->couplingManager().bulkCouplingContext(element, scvf)[0];
        const Scalar velocity = elemFaceVars[scvf].velocitySelf();
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());

        for (auto&& scv : scvs(stokesContext.fvGeometry))
        {
            if (this->couplingManager().isCoupledDof(lowDimIdx, scv.dofIndex()))
            {
                const bool insideIsUpstream = sign(velocity) == scvf.directionSign();
                return massFlux_(bulkIdx, lowDimIdx, scvf, insideScv, scv, elemVolVars[insideScv], stokesContext.elemVolVars[scv], velocity, insideIsUpstream);
            }
        }
        DUNE_THROW(Dune::InvalidStateException, "No coupled scvf found");
    }

    /*!
     * \brief Returns the mass flux across the coupling boundary as seen from the pore-network domain.
     */
    template<bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar energyCouplingCondition(const Element<lowDimIdx>& element,
                                   const FVElementGeometry<lowDimIdx>& fvGeometry,
                                   const ElementVolumeVariables<lowDimIdx>& elemVolVars,
                                   const SubControlVolume<lowDimIdx>& scv) const
    {
        const auto& lowDimContext = this->couplingManager().lowDimCouplingContext(element, scv);
        Scalar energyFlux(0.0);

        for (const auto& i : lowDimContext)
        {
            const auto& velocity = i.velocity;
            const auto& bulkScvf = i.getBulkScvf();
            const Scalar area = bulkScvf.area();

            const bool lowDimIsUpstream = sign(velocity) != bulkScvf.directionSign();

            auto flux = energyFlux_(lowDimIdx, bulkIdx, bulkScvf, scv, i.fvGeometry.scv(bulkScvf.insideScvIdx()), elemVolVars[scv], i.volVars, velocity, lowDimIsUpstream);
            flux *= area * i.volVars.extrusionFactor();
            flux *= -1.0; // flip the sign because the flux is used as a source term (different sign convention)

            energyFlux += flux;
        }

        return energyFlux;
    }

    /*!
     * \brief Returns the mass flux across the coupling boundary as seen from the free-flow domain.
     */
    template<class ElementVolumeVariables, class ElementFaceVariables, bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar energyCouplingCondition(const Element<bulkIdx>& element,
                                   const FVElementGeometry<bulkIdx>& fvGeometry,
                                   const ElementVolumeVariables& stokesElemVolVars,
                                   const ElementFaceVariables& stokesElemFaceVars,
                                   const SubControlVolumeFace<bulkIdx>& scvf) const
    {
        const auto& stokesContext = this->couplingManager().bulkCouplingContext(element, scvf)[0];
        const Scalar velocity = stokesElemFaceVars[scvf].velocitySelf();
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());

        for(auto&& scv : scvs(stokesContext.fvGeometry))
        {
            if(this->couplingManager().isCoupledDof(lowDimIdx, scv.dofIndex()))
            {
                const bool insideIsUpstream = sign(velocity) == scvf.directionSign();
                return energyFlux_(bulkIdx, lowDimIdx, scvf, insideScv, scv, stokesElemVolVars[insideScv], stokesContext.elemVolVars[scv], velocity, insideIsUpstream);
            }
        }
        DUNE_THROW(Dune::InvalidStateException, "No coupled scvf found");
    }

private:
    /*!
     * \brief Evaluate the compositional mole/mass flux across the interface.
     */
    template<std::size_t i, std::size_t j>
    NumEqVector massFlux_(Dune::index_constant<i> domainI,
                          Dune::index_constant<j> domainJ,
                          const SubControlVolumeFace<bulkIdx>& scvf,
                          const SubControlVolume<i>& insideScv,
                          const SubControlVolume<j>& outsideScv,
                          const VolumeVariables<i>& insideVolVars,
                          const VolumeVariables<j>& outsideVolVars,
                          const Scalar velocity,
                          const bool insideIsUpstream) const
    {
        NumEqVector flux(0.0);

        auto moleOrMassFraction = [&](const auto& volVars, int phaseIdx, int compIdx)
        { return useMoles ? volVars.moleFraction(phaseIdx, compIdx) : volVars.massFraction(phaseIdx, compIdx); };

        auto moleOrMassDensity = [&](const auto& volVars, int phaseIdx)
        { return useMoles ? volVars.molarDensity(phaseIdx) : volVars.density(phaseIdx); };

        // treat the advective fluxes
        auto insideTerm = [&](int compIdx)
        { return moleOrMassFraction(insideVolVars, couplingPhaseIdx(domainI), compIdx) * moleOrMassDensity(insideVolVars, couplingPhaseIdx(domainI)); };

        auto outsideTerm = [&](int compIdx)
        { return moleOrMassFraction(outsideVolVars, couplingPhaseIdx(domainJ), compIdx) * moleOrMassDensity(outsideVolVars, couplingPhaseIdx(domainJ)); };

        static const bool debugOutput = getParam<bool>("CouplingManager.DebugOutputMass", false);

        // TODO remove later
        if (debugOutput)
        {
            std::cout << "in domain " << i << std::endl;
            std::cout << "domainI phase " << FluidSystem<i>::phaseName(couplingPhaseIdx(domainI)) << ", domainJ phase " << FluidSystem<j>::phaseName(couplingPhaseIdx(domainJ)) << std::endl;
        }

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            const int domainICompIdx = couplingCompIdx(domainI, compIdx);
            const int domainJCompIdx = couplingCompIdx(domainJ, compIdx);

            assert(FluidSystem<i>::componentName(domainICompIdx) == FluidSystem<j>::componentName(domainJCompIdx));

            if (debugOutput)
                std::cout << "domainICompIdx " << domainICompIdx << ", " << FluidSystem<i>::componentName(domainICompIdx) << ", domainJCompIdx " << domainJCompIdx << ", " << FluidSystem<j>::componentName(domainJCompIdx) << std::endl;

            const Scalar sign = (domainI == lowDimIdx) ? -scvf.directionSign() : scvf.directionSign(); // compute the flux towards the interface
            flux[domainICompIdx] += sign * this->advectiveFlux(insideTerm(domainICompIdx), outsideTerm(domainJCompIdx), velocity, insideIsUpstream);
        }

        if (debugOutput)
        {
            std::cout << "advectiveFlux: " << flux << ", diffusive flux  " << diffusiveMolecularFlux_(domainI, domainJ, scvf, insideScv, outsideScv, insideVolVars, outsideVolVars) << std::endl;
            std::cout << "ratio adv/diff " << flux[0]/diffusiveMolecularFlux_(domainI, domainJ, scvf, insideScv, outsideScv, insideVolVars, outsideVolVars)[0] << std::endl;
            std::cout << "ratio adv/diff " << flux[1]/diffusiveMolecularFlux_(domainI, domainJ, scvf, insideScv, outsideScv, insideVolVars, outsideVolVars)[1] << std::endl;
        }

        NumEqVector diffusiveFlux = diffusiveMolecularFlux_(domainI, domainJ, scvf, insideScv, outsideScv, insideVolVars, outsideVolVars);

        //convert to correct units if necessary
        if (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged && useMoles)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                const int domainICompIdx = couplingCompIdx(domainI, compIdx);
                diffusiveFlux[domainICompIdx] /= FluidSystem<i>::molarMass(domainICompIdx);
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

        if (debugOutput)
            std::cout << "total flux " << flux << std::endl;

        // convert to total mass/mole balance, if set be user
        if (replaceCompEqIdx < numComponents)
            flux[replaceCompEqIdx] = std::accumulate(flux.begin(), flux.end(), 0.0);

        return flux;
    }

    /*!
     * \brief Evaluate the diffusive mole/mass flux across the interface.
     */
    template<std::size_t i, std::size_t j, bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar energyFlux_(Dune::index_constant<i> domainI,
                       Dune::index_constant<j> domainJ,
                       const SubControlVolumeFace<bulkIdx>& scvf,
                       const SubControlVolume<i>& scvI,
                       const SubControlVolume<j>& scvJ,
                       const VolumeVariables<i>& insideVolVars,
                       const VolumeVariables<j>& outsideVolVars,
                       const Scalar velocity,
                       const bool insideIsUpstream) const
    {
        Scalar flux(0.0);

        // convective fluxes
        const Scalar insideTerm = insideVolVars.density(couplingPhaseIdx(domainI)) * insideVolVars.enthalpy(couplingPhaseIdx(domainI));
        const Scalar outsideTerm = outsideVolVars.density(couplingPhaseIdx(domainJ)) * outsideVolVars.enthalpy(couplingPhaseIdx(domainJ));

        static const bool debugOutput = getParam<bool>("CouplingManager.DebugOutputEnergy", false);

        // TODO remove later
        if (debugOutput)
        {
            std::cout << "in domain " << i << std::endl;
            std::cout << "domainI phase " << FluidSystem<i>::phaseName(couplingPhaseIdx(domainI)) << ", domainJ phase " << FluidSystem<j>::phaseName(couplingPhaseIdx(domainJ)) << std::endl;
        }

        const Scalar sign = (domainI == lowDimIdx) ? -scvf.directionSign() : scvf.directionSign(); // compute the flux towards the interface
        flux += sign * this->advectiveFlux(insideTerm, outsideTerm, velocity, insideIsUpstream);

        if (debugOutput)
        {
            std::cout << "advective energy flux: " << flux << ", conductive flux  " << this->conductiveEnergyFlux_(domainI, domainJ, scvf, scvI, scvJ, insideVolVars, outsideVolVars) << std::endl;
            std::cout << "ratio adv/diff " << flux/this->conductiveEnergyFlux_(domainI, domainJ, scvf, scvI, scvJ, insideVolVars, outsideVolVars) << std::endl;
        }

        flux += this->conductiveEnergyFlux_(domainI, domainJ, scvf, scvI, scvJ, insideVolVars, outsideVolVars);

        auto diffusiveFlux = diffusiveMolecularFlux_(domainI, domainJ, scvf, scvI, scvJ, insideVolVars, outsideVolVars);
        for (int compIdx = 0; compIdx < diffusiveFlux.size(); ++compIdx)
        {
            const int domainICompIdx = couplingCompIdx(domainI, compIdx);
            const int domainJCompIdx = couplingCompIdx(domainJ, compIdx);

            const bool insideDiffFluxIsUpstream = diffusiveFlux[domainICompIdx] > 0.0;
            const Scalar componentEnthalpy = insideDiffFluxIsUpstream ?
                                             getComponentEnthalpy_(insideVolVars, couplingPhaseIdx(domainI), domainICompIdx)
                                           : getComponentEnthalpy_(outsideVolVars, couplingPhaseIdx(domainJ), domainJCompIdx);

           if (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged)
               flux += diffusiveFlux[domainICompIdx] * componentEnthalpy;
           else
               flux += diffusiveFlux[domainICompIdx] * FluidSystem<i>::molarMass(domainICompIdx) * componentEnthalpy;

            if (debugOutput)
                std::cout << "diffusive energy transfer for " << FluidSystem<i>::componentName(domainICompIdx)  << ": " << diffusiveFlux[domainICompIdx] * componentEnthalpy << std::endl;
        }

        if (debugOutput)
            std::cout << "total energy flux " << flux << std::endl;

        return flux;
    }

    /*!
     * \brief Evaluate the diffusive mole/mass flux across the interface.
     */
    template<std::size_t i, std::size_t j>
    NumEqVector diffusiveMolecularFlux_(Dune::index_constant<i> domainI,
                                        Dune::index_constant<j> domainJ,
                                        const SubControlVolumeFace<bulkIdx>& scvf,
                                        const SubControlVolume<i>& scvI,
                                        const SubControlVolume<j>& scvJ,
                                        const VolumeVariables<i>& volVarsI,
                                        const VolumeVariables<j>& volVarsJ) const
    {
        NumEqVector diffusiveFlux(0.0);
        const Scalar avgDensity = 0.5*(massOrMolarDensity(volVarsI, referenceSystemFormulation, couplingPhaseIdx(domainI))
                                     + massOrMolarDensity(volVarsJ, referenceSystemFormulation, couplingPhaseIdx(domainJ)));

        const auto& bulkVolVars = this->getBulkVolVars_(volVarsI, volVarsJ);
        const Scalar distance = this->getDistance_(scvI, scvJ, scvf);

        for (int compIdx = 1; compIdx < numComponents; ++compIdx)
        {
            const int bulkMainCompIdx = couplingPhaseIdx(bulkIdx);
            const int domainICompIdx = couplingCompIdx(domainI, compIdx);
            const int domainJCompIdx = couplingCompIdx(domainJ, compIdx);

            assert(FluidSystem<i>::componentName(domainICompIdx) == FluidSystem<j>::componentName(domainJCompIdx));

            const Scalar massOrMoleFractionI = massOrMoleFraction(volVarsI, referenceSystemFormulation, couplingPhaseIdx(domainI), domainICompIdx);
            const Scalar massOrMoleFractionJ = massOrMoleFraction(volVarsJ, referenceSystemFormulation, couplingPhaseIdx(domainJ), domainJCompIdx);
            const Scalar deltaMassOrMoleFrac = massOrMoleFractionJ - massOrMoleFractionI;

            const Scalar tij = bulkVolVars.effectiveDiffusionCoefficient(couplingPhaseIdx(bulkIdx), bulkMainCompIdx, couplingCompIdx(bulkIdx, compIdx)) / distance;
            diffusiveFlux[domainICompIdx] += -avgDensity * tij * deltaMassOrMoleFrac;
        }

        const Scalar cumulativeFlux = std::accumulate(diffusiveFlux.begin(), diffusiveFlux.end(), 0.0);
        diffusiveFlux[couplingCompIdx(domainI, 0)] = -cumulativeFlux;

        return diffusiveFlux;
    }

    Scalar getComponentEnthalpy_(const VolumeVariables<bulkIdx>& volVars, int phaseIdx, int compIdx) const
    {
        return FluidSystem<bulkIdx>::componentEnthalpy(volVars.fluidState(), 0, compIdx);
    }

    Scalar getComponentEnthalpy_(const VolumeVariables<lowDimIdx>& volVars, int phaseIdx, int compIdx) const
    {
        return FluidSystem<lowDimIdx>::componentEnthalpy(volVars.fluidState(), phaseIdx, compIdx);
    }
};

} // end namespace Dumux

#endif
