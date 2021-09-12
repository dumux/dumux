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
 * \ingroup MultiDomain
 * \copydoc Dumux::MultiBinaryCouplingManager
 */

#ifndef DUMUX_MULTIDOMAIN_MULTIBINARY_COUPLINGMANAGER_HH
#define DUMUX_MULTIDOMAIN_MULTIBINARY_COUPLINGMANAGER_HH

#include <utility>
#include <memory>
#include <dune/common/hybridutilities.hh>
#include <dumux/common/properties.hh>
#include <dumux/multidomain/traits.hh>

namespace Dumux {

namespace Detail {

template <std::size_t, typename Tuple>
struct HasIndex;

template <std::size_t i, typename... Indices>
struct HasIndex<i, std::tuple<Indices...>>
: std::disjunction<std::is_same<Dune::index_constant<i>, Indices>...>
{};

template<class Map, std::size_t i, std::size_t j>
constexpr bool isCoupled(Dune::index_constant<i> I, Dune::index_constant<j>)
{ return HasIndex<j, std::decay_t<decltype(Map::coupledDomains(I))>>::value; }

} // end namespace Detail

/*!
 * \ingroup MultiDomain
 * \brief Coupling manager that combines an arbitrary number of binary coupling manager (coupling two domains each)
 */
template<class MDTraits, class CouplingManagerMaps, class ...CouplingMgrs>
class MultiBinaryCouplingManager
{
    // the sub domain type tags
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolumeFace = typename FVElementGeometry<id>::SubControlVolumeFace;
    template<std::size_t id> using SubControlVolume = typename FVElementGeometry<id>::SubControlVolume;
    using CouplingStencil = std::vector<std::size_t>;

    template<std::size_t id>
    using SubCouplingManagerT = typename std::tuple_element_t<id, std::tuple<CouplingMgrs...>>;

    using CMIndices = std::make_index_sequence<sizeof...(CouplingMgrs)>;
    using CouplingManagers = typename Detail::MultiDomainTupleSharedPtr<SubCouplingManagerT, CMIndices>::type;

    template<std::size_t id>
    using SubSolutionVector = std::decay_t<decltype(std::declval<typename MDTraits::SolutionVector>()[Dune::index_constant<id>()])>;
    using SolutionVectors = typename MDTraits::template TupleOfSharedPtr<SubSolutionVector>;

public:

    static constexpr auto couplingManagerMap = CouplingManagerMaps::makeCouplingManagerMap();

    template<std::size_t i, std::size_t j>
    using SubCouplingManager = SubCouplingManagerT<couplingManagerMap[i][j]>;

    MultiBinaryCouplingManager()
    {
        using namespace Dune::Hybrid;
        forEach(couplingManagers_, [&](auto&& couplingManager)
        {
            couplingManager = std::make_shared<typename std::decay_t<decltype(couplingManager)>::element_type>();
        });

        forEach(solutionVectors_, [&](auto&& solutionVector)
        {
            solutionVector = std::make_shared<typename std::decay_t<decltype(solutionVector)>::element_type>();
        });
    }

    //! Returns the local domain indices used within the binary sub coupling managers.
    template<std::size_t i, std::size_t j>
    static constexpr auto globalToLocalDomainIndices(Dune::index_constant<i> domainI, Dune::index_constant<j> domainJ)
    {
        static_assert(i <= MDTraits::numSubDomains && j <= MDTraits::numSubDomains);
        return CouplingManagerMaps::globalToLocalDomainIndices(domainI, domainJ);
    }

    //! Returns the coupling manager index for a given domain combination
    template<std::size_t i, std::size_t j>
    static constexpr auto getCouplingManagerIndex(Dune::index_constant<i> domainI, Dune::index_constant<j> domainJ)
    {
        static_assert(
            Detail::template isCoupled<CouplingManagerMaps>(domainI, domainJ),
            "Sub-coupling manager only exists for coupled domains."
        );
        return couplingManagerMap[i][j];
    }

    //! return the binary sub-coupling manager
    template<std::size_t i, std::size_t j>
    auto& subCouplingManager(Dune::index_constant<i> domainI,
                             Dune::index_constant<j> domainJ)
    {
        constexpr auto idx = getCouplingManagerIndex(domainI, domainJ);
        return *std::get<idx>(couplingManagers_);
    }

    //! return the binary sub-coupling manager
    template<std::size_t i, std::size_t j>
    const auto& subCouplingManager(Dune::index_constant<i> domainI,
                                   Dune::index_constant<j> domainJ) const
    {
        constexpr auto idx = getCouplingManagerIndex(domainI, domainJ);
        return *std::get<idx>(couplingManagers_);
    }

    //! Update the solution vector before assembly
    void updateSolution(const typename MDTraits::SolutionVector& curSol)
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(solutionVectors_)), [&](const auto id)
        {
            *std::get<id>(solutionVectors_) = curSol[id];
        });
    }

    /*!
     * \brief extend the jacobian pattern of the diagonal block of domain i
     *        by those entries that are not already in the uncoupled pattern
     * \note per default we do not add such additional dependencies
     * \note Such additional dependencies can arise from the coupling, e.g. if a coupling source
     *       term depends on a non-local average of a quantity of the same domain
     * \warning if you overload this also implement evalAdditionalDomainDerivatives
     */
    template<std::size_t id, class JacobianPattern>
    void extendJacobianPattern(Dune::index_constant<id> domainI, JacobianPattern& pattern) const
    {}

    /*!
     * \brief Return the coupling element stencil for a given bulk domain element
     */
    template<std::size_t i, class Entity, std::size_t j>
    const auto& couplingStencil(Dune::index_constant<i> domainI,
                                const Entity& entity,
                                Dune::index_constant<j> domainJ) const
    {
        // if the domains are coupled according to the map, forward to sub-coupling manager
        if constexpr (Detail::template isCoupled<CouplingManagerMaps>(domainI, domainJ))
            return subCouplingManager(domainI, domainJ).couplingStencil(
                globalToLocalDomainIndices(domainI, domainJ).domainI, entity,
                globalToLocalDomainIndices(domainI, domainJ).domainJ
            );
        else
            return emptyStencil_;
    }

    // ! evaluate coupling residual for the derivative low dim DOF with respect to bulk DOF
    // ! we only need to evaluate the part of the residual that will be influence by the bulk DOF
    template<std::size_t i, class LocalAssemblerI, std::size_t j>
    decltype(auto) evalCouplingResidual(Dune::index_constant<i> domainI,
                                        const SubControlVolumeFace<i>& scvfI,
                                        const LocalAssemblerI& localAssemblerI,
                                        Dune::index_constant<j> domainJ,
                                        std::size_t dofIdxGlobalJ) const
    {
        if constexpr (i == j || Detail::template isCoupled<CouplingManagerMaps>(domainI, domainJ))
            return subCouplingManager(domainI, domainJ).evalCouplingResidual(
                globalToLocalDomainIndices(domainI, domainJ).domainI,
                scvfI, localAssemblerI,
                globalToLocalDomainIndices(domainI, domainJ).domainJ,
                dofIdxGlobalJ
            );
        else
        {
            DUNE_THROW(Dune::InvalidStateException,
                "Calling evalCouplingResidual for uncoupled domains " << i << " and " << j
            );
            return localAssemblerI.evalLocalResidual();
        }
    }

    //! evaluate coupling residual for the derivative low dim DOF with respect to bulk DOF
    //! we only need to evaluate the part of the residual that will be influence by the bulk DOF
    template<std::size_t i, class LocalAssemblerI, std::size_t j>
    decltype(auto) evalCouplingResidual(Dune::index_constant<i> domainI,
                                        const LocalAssemblerI& localAssemblerI,
                                        Dune::index_constant<j> domainJ,
                                        std::size_t dofIdxGlobalJ) const
    {
        if constexpr (i == j || Detail::template isCoupled<CouplingManagerMaps>(domainI, domainJ))
            return subCouplingManager(domainI, domainJ).evalCouplingResidual(
                globalToLocalDomainIndices(domainI, domainJ).domainI,
                localAssemblerI,
                globalToLocalDomainIndices(domainI, domainJ).domainJ,
                dofIdxGlobalJ
            );
        else
        {
            DUNE_THROW(Dune::InvalidStateException,
                "Calling evalCouplingResidual for uncoupled domains " << i << " and " << j
            );
            return localAssemblerI.evalLocalResidual();
        }
    }

    //! evaluate coupling residual for the derivative low dim DOF with respect to bulk DOF
    //! we only need to evaluate the part of the residual that will be influence by the bulk DOF
    template<std::size_t i, class LocalAssemblerI, std::size_t j>
    decltype(auto) evalCouplingResidual(Dune::index_constant<i> domainI,
                                        const LocalAssemblerI& localAssemblerI,
                                        const SubControlVolume<i>& scvI,
                                        Dune::index_constant<j> domainJ,
                                        std::size_t dofIdxGlobalJ) const
    {
        if constexpr (i == j || Detail::template isCoupled<CouplingManagerMaps>(domainI, domainJ))
            return subCouplingManager(domainI, domainJ).evalCouplingResidual(
                globalToLocalDomainIndices(domainI, domainJ).domainI,
                localAssemblerI,
                scvI,
                globalToLocalDomainIndices(domainI, domainJ).domainJ,
                dofIdxGlobalJ
            );
        else
        {
            DUNE_THROW(Dune::InvalidStateException,
                "Calling evalCouplingResidual for uncoupled domains " << i << " and " << j
            );
            return localAssemblerI.evalLocalResidual();
        }
    }

    /*!
     * \brief Update the coupling context for the bulk face residual w.r.t to the lowDim dofs
     */
    template<std::size_t i, class LocalAssemblerI, std::size_t j, class PrimaryVariables>
    void updateCouplingContext(Dune::index_constant<i> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<j> domainJ,
                               const std::size_t dofIdxGlobalJ,
                               const PrimaryVariables& priVars,
                               int pvIdxJ)
    {
        // only one other manager needs to take action if i != j
        if constexpr (i != j)
        {
            if constexpr (Detail::template isCoupled<CouplingManagerMaps>(domainI, domainJ))
                subCouplingManager(domainI, domainJ).updateCouplingContext(
                    globalToLocalDomainIndices(domainI, domainJ).domainI,
                    localAssemblerI,
                    globalToLocalDomainIndices(domainI, domainJ).domainJ,
                    dofIdxGlobalJ,
                    priVars,
                    pvIdxJ
                );
        }
        else
        {
            // for i == j, we need to call all relevant managers where domainI is involved and make sure that the
            // context is updated with respect to its own domain
            using namespace Dune::Hybrid;
            forEach(CouplingManagerMaps::coupledDomains(domainI), [&](const auto domainJ)
            {
                static_assert(domainI != domainJ);
                subCouplingManager(domainI, domainJ).updateCouplingContext(
                    globalToLocalDomainIndices(domainI, domainJ).domainI,
                    localAssemblerI,
                    globalToLocalDomainIndices(domainI, domainJ).domainI,
                    dofIdxGlobalJ,
                    priVars,
                    pvIdxJ
                );
            });
        }
    }

    //! Bind the coupling context for a low dim element TODO remove Assembler
    template<std::size_t i, class Assembler = int>
    void bindCouplingContext(Dune::index_constant<i> domainI, const Element<i>& element, const Assembler& assembler = 0)
    {
        // for the coupling blocks
        using namespace Dune::Hybrid;
        forEach(CouplingManagerMaps::coupledDomains(domainI), [&](const auto domainJ)
        {
            this->subCouplingManager(domainI, domainJ).bindCouplingContext(
                globalToLocalDomainIndices(domainI, domainJ).domainI, element, assembler
            );
        });
    }

    /*!
     * \brief return the numeric epsilon used for deflecting primary variables of coupled domain i.
     * \note  specialization for free-flow schemes
     */
    template<std::size_t i>
    decltype(auto) numericEpsilon(Dune::index_constant<i> domainI,
                                  const std::string& paramGroup) const
    {
        static constexpr auto someOtherDomain = std::get<0>(CouplingManagerMaps::coupledDomains(domainI));
        return subCouplingManager(domainI, someOtherDomain).numericEpsilon(globalToLocalDomainIndices(domainI, someOtherDomain).domainI, paramGroup);
    }

    /*!
     * \brief evaluate additional derivatives of the element residual of a domain with respect
     *        to dofs in the same domain that are not in the regular stencil (see CouplingManager::extendJacobianPattern)
     * \note Such additional dependencies can arise from the coupling, e.g. if a coupling source
     *       term depends on a non-local average of a quantity of the same domain
     */
    template<std::size_t i, class LocalAssemblerI, class JacobianMatrixDiagBlock, class GridVariables>
    void evalAdditionalDomainDerivatives(Dune::index_constant<i> domainI,
                                         const LocalAssemblerI& localAssemblerI,
                                         const typename LocalAssemblerI::LocalResidual::ElementResidualVector& origResiduals,
                                         JacobianMatrixDiagBlock& A,
                                         GridVariables& gridVariables)
    {}

    template<std::size_t i, class LocalAssemblerI, class UpdatableElementVolVars, class UpdatableFluxVarCache>
    void updateCoupledVariables(Dune::index_constant<i> domainI,
                                const LocalAssemblerI& localAssemblerI,
                                UpdatableElementVolVars& elemVolVars,
                                UpdatableFluxVarCache& elemFluxVarsCache)
    {
        // for the coupling blocks
        using namespace Dune::Hybrid;
        forEach(CouplingManagerMaps::coupledDomains(domainI), [&](const auto domainJ)
        {
            subCouplingManager(domainI, domainJ).updateCoupledVariables(
                globalToLocalDomainIndices(domainI, domainJ).domainI,
                localAssemblerI,
                elemVolVars,
                elemFluxVarsCache
            );
        });
    }

protected:
    SolutionVectors& curSol()
    { return solutionVectors_; }

    const SolutionVectors& curSol() const
    { return solutionVectors_; }

private:
    CouplingManagers couplingManagers_;
    SolutionVectors solutionVectors_;

    CouplingStencil emptyStencil_;
};

} // end namespace Dumux

#endif
