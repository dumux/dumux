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
 * \copydoc Dumux::FreeFlowPorousMediumCouplingManager
 */

#ifndef DUMUX_MULTIDOMAIN_BOUNDARY_FREEFLOW_POROUSMEDIUM_COUPLINGMANAGER_HH
#define DUMUX_MULTIDOMAIN_BOUNDARY_FREEFLOW_POROUSMEDIUM_COUPLINGMANAGER_HH

#include <utility>
#include <memory>

#include <dune/common/indices.hh>
#include <dumux/common/properties.hh>
#include <dumux/multidomain/boundary/freeflowporousmedium/ffmasspm/couplingmanager.hh>
#include <dumux/multidomain/boundary/freeflowporousmedium/ffmomentumpm/couplingmanager.hh>
#include <dumux/multidomain/staggeredfreeflow/couplingmanager.hh>
#include <dumux/multidomain/multibinarycouplingmanager.hh>

#include "couplingconditions.hh"

namespace Dumux {

namespace FreeFlowPorousMediumDetail {

// global subdomain indices
static constexpr auto freeFlowMomentumIndex = Dune::index_constant<0>();
static constexpr auto freeFlowMassIndex = Dune::index_constant<1>();
static constexpr auto porousMediumIndex = Dune::index_constant<2>();

// coupling indices
static constexpr auto freeFlowMassToFreeFlowMomentumIndex = Dune::index_constant<0>();
static constexpr auto freeFlowMomentumToPorousMediumIndex = Dune::index_constant<1>();
static constexpr auto freeFlowMassToPorousMediumIndex = Dune::index_constant<2>();
static constexpr auto noCouplingIdx = Dune::index_constant<99>();

constexpr auto makeCouplingManagerMap()
{
    auto map = std::array<std::array<std::size_t, 3>, 3>{};

    // free flow (momentum-mass)
    map[freeFlowMomentumIndex][freeFlowMassIndex] = freeFlowMassToFreeFlowMomentumIndex;
    map[freeFlowMassIndex][freeFlowMomentumIndex] = freeFlowMassToFreeFlowMomentumIndex;

    // free flow momentum - porous medium
    map[freeFlowMomentumIndex][porousMediumIndex] = freeFlowMomentumToPorousMediumIndex;
    map[porousMediumIndex][freeFlowMomentumIndex] = freeFlowMomentumToPorousMediumIndex;

    // free flow mass - porous medium
    map[freeFlowMassIndex][porousMediumIndex] = freeFlowMassToPorousMediumIndex;
    map[porousMediumIndex][freeFlowMassIndex] = freeFlowMassToPorousMediumIndex;

    return map;
}

template<std::size_t i>
constexpr auto coupledDomains(Dune::index_constant<i> domainI)
{
    if constexpr (i == freeFlowMomentumIndex)
        return std::make_tuple(freeFlowMassIndex, porousMediumIndex);
    else if constexpr (i == freeFlowMassIndex)
        return std::make_tuple(freeFlowMomentumIndex, porousMediumIndex);
    else // i == porousMediumIndex
        return std::make_tuple(freeFlowMomentumIndex, freeFlowMassIndex);
}

template<std::size_t i, std::size_t j>
constexpr auto globalToLocalDomainIndices(Dune::index_constant<i>, Dune::index_constant<j>)
{
    static_assert(i <= 2 && j <= 2);
    static_assert(i != j);

    if constexpr (i < j)
        return std::pair<Dune::index_constant<0>, Dune::index_constant<1>>{};
    else
        return std::pair<Dune::index_constant<1>, Dune::index_constant<0>>{};
}

struct CouplingMaps
{
    static constexpr auto managerMap()
    {
        return  FreeFlowPorousMediumDetail::makeCouplingManagerMap();
    }

    template<std::size_t i, std::size_t j>
    static constexpr auto globalToLocal(Dune::index_constant<i> domainI, Dune::index_constant<j> domainJ)
    {
        return FreeFlowPorousMediumDetail::globalToLocalDomainIndices(domainI, domainJ);
    }

    template<std::size_t i>
    static constexpr auto coupledDomains(Dune::index_constant<i> domainI)
    {
        return FreeFlowPorousMediumDetail::coupledDomains(domainI);
    }
};

template<class MDTraits>
struct CouplingManagers
{
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    using FreeFlowTraits = MultiDomainTraits<
        SubDomainTypeTag<freeFlowMomentumIndex>, SubDomainTypeTag<freeFlowMassIndex>
    >;

    using FreeFlowMomentumPorousMediumTraits = MultiDomainTraits<
        SubDomainTypeTag<freeFlowMomentumIndex>, SubDomainTypeTag<porousMediumIndex>
    >;

    using FreeFlowMassPorousMediumTraits = MultiDomainTraits<
        SubDomainTypeTag<freeFlowMassIndex>, SubDomainTypeTag<porousMediumIndex>
    >;

    using FreeFlowCouplingManager
        = Dumux::StaggeredFreeFlowCouplingManager<FreeFlowTraits>;
    using FreeFlowMomentumPorousMediumCouplingManager
        = Dumux::FreeFlowMomentumPorousMediumCouplingManager<FreeFlowMomentumPorousMediumTraits>;
    using FreeFlowMassPorousMediumCouplingManager
        = Dumux::FreeFlowMassPorousMediumCouplingManager<FreeFlowMassPorousMediumTraits>;
};

} // end namespace Detail


/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief Coupling manager for coupling freeflow and porous medium flow models
 */
template<class MDTraits>
class FreeFlowPorousMediumCouplingManager
: public MultiBinaryCouplingManager<
    MDTraits,
    FreeFlowPorousMediumDetail::CouplingMaps,
    typename FreeFlowPorousMediumDetail::CouplingManagers<MDTraits>::FreeFlowCouplingManager,
    typename FreeFlowPorousMediumDetail::CouplingManagers<MDTraits>::FreeFlowMomentumPorousMediumCouplingManager,
    typename FreeFlowPorousMediumDetail::CouplingManagers<MDTraits>::FreeFlowMassPorousMediumCouplingManager
>
{
    using ParentType = MultiBinaryCouplingManager<
        MDTraits,
        FreeFlowPorousMediumDetail::CouplingMaps,
        typename FreeFlowPorousMediumDetail::CouplingManagers<MDTraits>::FreeFlowCouplingManager,
        typename FreeFlowPorousMediumDetail::CouplingManagers<MDTraits>::FreeFlowMomentumPorousMediumCouplingManager,
        typename FreeFlowPorousMediumDetail::CouplingManagers<MDTraits>::FreeFlowMassPorousMediumCouplingManager
    >;

    using Scalar = typename MDTraits::Scalar;

    // the sub domain type tags
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolumeFace = typename FVElementGeometry<id>::SubControlVolumeFace;
    template<std::size_t id> using SubControlVolume = typename FVElementGeometry<id>::SubControlVolume;
    template<std::size_t id> using ElementVolumeVariables = typename GetPropType<SubDomainTypeTag<id>, Properties::GridVolumeVariables>::LocalView;
    template<std::size_t id> using NumEqVector = typename Problem<id>::Traits::NumEqVector;

    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    using SolutionVector = typename MDTraits::SolutionVector;

    using CouplingConditions = FreeFlowPorousMediumCouplingConditions<MDTraits, FreeFlowPorousMediumCouplingManager<MDTraits>>;

public:

    template<std::size_t i, std::size_t j>
    using SubCouplingManager = typename ParentType::template SubCouplingManager<i, j>;

    static constexpr auto freeFlowMomentumIndex = FreeFlowPorousMediumDetail::freeFlowMomentumIndex;
    static constexpr auto freeFlowMassIndex = FreeFlowPorousMediumDetail::freeFlowMassIndex;
    static constexpr auto porousMediumIndex = FreeFlowPorousMediumDetail::porousMediumIndex;

public:
    using ParentType::ParentType;

    template<class GridVarsTuple>
    void init(std::shared_ptr<Problem<freeFlowMomentumIndex>> freeFlowMomentumProblem,
              std::shared_ptr<Problem<freeFlowMassIndex>> freeFlowMassProblem,
              std::shared_ptr<Problem<porousMediumIndex>> porousMediumProblem,
              GridVarsTuple&& gridVarsTuple,
              const SolutionVector& curSol)
    {
        this->updateSolution(curSol); // generic coupling manager stores tuple of shared_ptr

        // initialize the binary sub coupling managers
        typename SubCouplingManager<freeFlowMomentumIndex, freeFlowMassIndex>::SolutionVectorStorage ffSolVecTuple;
        std::get<0>(ffSolVecTuple) = std::get<freeFlowMomentumIndex>(this->curSol());
        std::get<1>(ffSolVecTuple) = std::get<freeFlowMassIndex>(this->curSol());
        this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).init(
            freeFlowMomentumProblem, freeFlowMassProblem,
            std::make_tuple(std::get<freeFlowMomentumIndex>(gridVarsTuple), std::get<freeFlowMassIndex>(gridVarsTuple)),
            ffSolVecTuple
        );

        typename SubCouplingManager<freeFlowMassIndex, porousMediumIndex>::SolutionVectorStorage ffMassPmSolVecTuple;
        std::get<0>(ffMassPmSolVecTuple) = std::get<freeFlowMassIndex>(this->curSol());
        std::get<1>(ffMassPmSolVecTuple) = std::get<porousMediumIndex>(this->curSol());
        this->subCouplingManager(freeFlowMassIndex, porousMediumIndex).init(
            freeFlowMassProblem, porousMediumProblem, ffMassPmSolVecTuple
        );

        typename SubCouplingManager<freeFlowMomentumIndex, porousMediumIndex>::SolutionVectorStorage ffMomentumPmSolVecTuple;
        std::get<0>(ffMomentumPmSolVecTuple) = std::get<freeFlowMomentumIndex>(this->curSol());
        std::get<1>(ffMomentumPmSolVecTuple) = std::get<porousMediumIndex>(this->curSol());
        this->subCouplingManager(freeFlowMomentumIndex, porousMediumIndex).init(
            freeFlowMomentumProblem, porousMediumProblem, ffMomentumPmSolVecTuple
        );
    }

    /*!
     * \brief Returns the mass flux across the coupling boundary.
     */
    template<std::size_t i, std::size_t j>
    auto massCouplingCondition(Dune::index_constant<i> domainI, Dune::index_constant<j> domainJ,
                               const FVElementGeometry<i>& fvGeometry,
                               const typename FVElementGeometry<i>::SubControlVolumeFace& scvf,
                               const ElementVolumeVariables<i>& elemVolVars) const
    {
        static_assert(domainI != freeFlowMomentumIndex && domainJ != freeFlowMomentumIndex);

        const auto& couplingContext = this->subApply(domainI, domainJ, [&](const auto& cm, auto&& ii, auto&& jj) -> const auto& {
            return cm.couplingContext(ii, fvGeometry, scvf);
        });

        const auto& freeFlowElement = [&]
        {
            if constexpr (domainI == freeFlowMassIndex)
                return fvGeometry.element();
            else
                return couplingContext.fvGeometry.element();
        }();

        const auto& freeFlowScvf = [&]
        {
            if constexpr (domainI == freeFlowMassIndex)
                return scvf;
            else
                return couplingContext.fvGeometry.scvf(couplingContext.freeFlowMassScvfIdx);

        }();

        // todo revise velocity (see ff mom pm mgr)

        couplingContext.velocity = this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).faceVelocity(freeFlowElement, freeFlowScvf);
        return CouplingConditions::massCouplingCondition(domainI, domainJ, fvGeometry, scvf, elemVolVars, couplingContext);
    }


    //////////////////////// Conditions for FreeFlowMomentum - PorousMedium coupling //////////
    ///////////////////////////////////////////////////////////////////////////////////////////

    NumEqVector<freeFlowMomentumIndex> momentumCouplingCondition(Dune::index_constant<freeFlowMomentumIndex> domainI,
                                                                 Dune::index_constant<porousMediumIndex> domainJ,
                                                                 const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                                                                 const typename FVElementGeometry<freeFlowMomentumIndex>::SubControlVolumeFace& scvf,
                                                                 const ElementVolumeVariables<freeFlowMomentumIndex>& elemVolVars) const
    {
        if (scvf.isLateral())
            return NumEqVector<freeFlowMomentumIndex>(0.0);

        const auto& context = this->subCouplingManager(freeFlowMomentumIndex, porousMediumIndex).couplingContext(
            domainI, fvGeometry, scvf
        );

        return CouplingConditions::momentumCouplingCondition(fvGeometry, scvf, elemVolVars, context);
    }

    /*!
     * \brief Returns the intrinsic permeability of the coupled Darcy element.
     */
    auto darcyPermeability(const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                           const SubControlVolumeFace<freeFlowMomentumIndex>& scvf) const
    {
        if (scvf.isFrontal())
        {
            const auto& context = this->subCouplingManager(freeFlowMomentumIndex, porousMediumIndex).couplingContext(
                Dune::index_constant<freeFlowMomentumIndex>(), fvGeometry, scvf
            );

            return CouplingConditions::darcyPermeability(fvGeometry, scvf, context);
        }
        else
        {
            const auto& orthogonalScvf = fvGeometry.lateralOrthogonalScvf(scvf);
            const auto& orthogonalScv = fvGeometry.scv(orthogonalScvf.insideScvIdx());
            const auto& frontalScvfOnBoundary = fvGeometry.frontalScvfOnBoundary(orthogonalScv);
            const auto& context = this->subCouplingManager(freeFlowMomentumIndex, porousMediumIndex).couplingContext(
                Dune::index_constant<freeFlowMomentumIndex>(), fvGeometry, frontalScvfOnBoundary
            );

            return CouplingConditions::darcyPermeability(fvGeometry, frontalScvfOnBoundary, context);
        }
    }

    //////////////////////// Conditions for FreeFlowMomentum - FreeFlowMass coupling //////////
    ///////////////////////////////////////////////////////////////////////////////////////////

    /*!
     * \brief Returns the pressure at a given sub control volume face.
     */
    Scalar pressure(const Element<freeFlowMomentumIndex>& element,
                    const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                    const SubControlVolumeFace<freeFlowMomentumIndex>& scvf) const
    {
        return this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).pressure(
            element, fvGeometry, scvf
        );
    }

    /*!
     * \brief Returns the pressure at the center of a sub control volume corresponding to a given sub control volume face.
     *        This is used for setting a Dirichlet pressure for the mass model when a fixed pressure for the momentum balance is set at another
     *        boundary. Since the the pressure at the given scvf is solution-dependent and thus unknown a priori, we just use the value
     *        of the interior cell here.
     */
    Scalar cellPressure(const Element<freeFlowMomentumIndex>& element,
                        const SubControlVolumeFace<freeFlowMomentumIndex>& scvf) const
    {
        return this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).cellPressure(
            element, scvf
        );
    }

    /*!
     * \brief Returns the density at a given sub control volume face.
     */
    Scalar density(const Element<freeFlowMomentumIndex>& element,
                   const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                   const SubControlVolumeFace<freeFlowMomentumIndex>& scvf,
                   const bool considerPreviousTimeStep = false) const
    {
        return this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).density(
            element, fvGeometry, scvf, considerPreviousTimeStep
        );
    }

    auto insideAndOutsideDensity(const Element<freeFlowMomentumIndex>& element,
                                 const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                                 const SubControlVolumeFace<freeFlowMomentumIndex>& scvf,
                                 const bool considerPreviousTimeStep = false) const
    {
        return this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).insideAndOutsideDensity(
            element, fvGeometry, scvf, considerPreviousTimeStep
        );
    }

    /*!
     * \brief Returns the density at a given sub control volume.
     */
    Scalar density(const Element<freeFlowMomentumIndex>& element,
                   const SubControlVolume<freeFlowMomentumIndex>& scv,
                   const bool considerPreviousTimeStep = false) const
    {
        return this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).density(
            element, scv, considerPreviousTimeStep
        );
    }

    /*!
     * \brief Returns the pressure at a given sub control volume face.
     */
    Scalar effectiveViscosity(const Element<freeFlowMomentumIndex>& element,
                              const FVElementGeometry<freeFlowMomentumIndex>& fvGeometry,
                              const SubControlVolumeFace<freeFlowMomentumIndex>& scvf) const
    {
        return this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).effectiveViscosity(
            element, fvGeometry, scvf
        );
    }

    /*!
     * \brief Returns the velocity at a given sub control volume face.
     */
    auto faceVelocity(const Element<freeFlowMassIndex>& element,
                      const SubControlVolumeFace<freeFlowMassIndex>& scvf) const
    {
        return this->subCouplingManager(freeFlowMomentumIndex, freeFlowMassIndex).faceVelocity(
            element, scvf
        );
    }

    template<std::size_t i>
    const Problem<i>& problem(Dune::index_constant<i> domainI) const
    {
        return this->subApply(domainI, [&](const auto& cm, auto&& ii) -> const auto& {
            return cm.problem(ii);
        });
    }

    template<std::size_t i, std::size_t j>
    bool isCoupled(Dune::index_constant<i> domainI,
                   Dune::index_constant<j> domainJ,
                   const SubControlVolumeFace<i>& scvf) const
    {
        return this->subApply(domainI, domainJ, [&](const auto& cm, auto&& ii, auto&& jj){
            return cm.isCoupled(ii, scvf);
        });
    }

    /*!
     * \brief If the boundary entity is on a coupling boundary
     * \param domainI the domain index of domain i for which to compute the flux
     * \param domainJ the domain index of domain j for which to compute the flux
     * \param scv the sub control volume
     */
    template<std::size_t i, std::size_t j>
    bool isCoupled(Dune::index_constant<i> domainI,
                   Dune::index_constant<j> domainJ,
                   const SubControlVolume<i>& scv) const
    {
        return this->subApply(domainI, domainJ, [&](const auto& cm, auto&& ii, auto&& jj){
            return cm.isCoupled(ii, scv);
        });
    }

    /*!
     * \brief Returns whether a given scvf is coupled to the other domain
     */
    bool isCoupledLateralScvf(Dune::index_constant<freeFlowMomentumIndex> domainI,
                              Dune::index_constant<porousMediumIndex> domainJ,
                              const SubControlVolumeFace<freeFlowMomentumIndex>& scvf) const
    {
        return this->subApply(domainI, domainJ, [&](const auto& cm, auto&& ii, auto&& jj){
            return cm.isCoupledLateralScvf(ii, scvf);
        });
    }


    using ParentType::couplingStencil;
    /*!
     * \brief returns an iterable container of all indices of degrees of freedom of domain j
     *        that couple with / influence the residual of the given sub-control volume of domain i
     *
     * \param domainI the domain index of domain i
     * \param elementI the coupled element of domain í
     * \param scvI the sub-control volume of domain i
     * \param domainJ the domain index of domain j
     */
    template<std::size_t j>
    const auto& couplingStencil(Dune::index_constant<freeFlowMomentumIndex> domainI,
                                const Element<freeFlowMomentumIndex>& elementI,
                                const SubControlVolume<freeFlowMomentumIndex>& scvI,
                                Dune::index_constant<j> domainJ) const
    {
        static_assert(freeFlowMomentumIndex != j);
        return this->subApply(domainI, domainJ, [&](const auto& cm, auto&& ii, auto&& jj) -> const auto& {
            return cm.couplingStencil(ii, elementI, scvI, jj);
        });
    }
};

} // end namespace Dumux

#endif
