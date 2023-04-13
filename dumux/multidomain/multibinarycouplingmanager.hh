// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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

} // end namespace Detail

/*!
 * \ingroup MultiDomain
 * \brief Coupling manager that combines an arbitrary number of binary coupling manager (coupling two domains each)
 * \tparam MDTraits the multidomain traits
 * \tparam CouplingMap a coupling policy class
 * \tparam CouplingMgrs the binary sub-coupling manager types
 *
 * The coupling policy has to provide the interfaces
 * - CouplingMap::coupledDomains(i): returns a tuple of Dune::index_constants with the coupled domains
 * - CouplingMap::globalToLocal(i, j): maps the indices i, j to the local index pair of the responsible sub coupling manager
 * - CouplingMap::managerMap(): returns a two-dimensional array mapping two indices to the coupling manager index
 */
template<class MDTraits, class CouplingMap, class ...CouplingMgrs>
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

    static constexpr auto couplingManagerMap_ = CouplingMap::managerMap();

    //! Returns the local domain indices used within the binary sub coupling managers.
    template<std::size_t i, std::size_t j>
    static constexpr auto globalToLocal_(Dune::index_constant<i> domainI, Dune::index_constant<j> domainJ)
    {
        static_assert(i <= MDTraits::numSubDomains && j <= MDTraits::numSubDomains);
        return CouplingMap::globalToLocal(domainI, domainJ);
    }

    //! If two domain are coupled
    template<class Map, std::size_t i, std::size_t j>
    static constexpr bool isCoupled_(Dune::index_constant<i> domainI, Dune::index_constant<j>)
    { return Detail::HasIndex<j, std::decay_t<decltype(Map::coupledDomains(domainI))>>::value; }

    //! Returns the coupling manager index for a given domain combination
    template<std::size_t i, std::size_t j>
    static constexpr auto subCouplingManagerIndex_(Dune::index_constant<i> domainI, Dune::index_constant<j> domainJ)
    {
        static_assert(
            isCoupled_<CouplingMap>(domainI, domainJ),
            "Sub-coupling manager only exists for coupled domains."
        );
        return couplingManagerMap_[i][j];
    }

public:
    template<std::size_t i, std::size_t j>
    using SubCouplingManager = SubCouplingManagerT<couplingManagerMap_[i][j]>;

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

    //! return the binary sub-coupling manager
    template<std::size_t i, std::size_t j>
    auto& subCouplingManager(Dune::index_constant<i> domainI,
                             Dune::index_constant<j> domainJ)
    {
        constexpr auto idx = subCouplingManagerIndex_(domainI, domainJ);
        return *std::get<idx>(couplingManagers_);
    }

    //! return the binary sub-coupling manager
    template<std::size_t i, std::size_t j>
    const auto& subCouplingManager(Dune::index_constant<i> domainI,
                                   Dune::index_constant<j> domainJ) const
    {
        constexpr auto idx = subCouplingManagerIndex_(domainI, domainJ);
        return *std::get<idx>(couplingManagers_);
    }

    //! apply a function to the domainI-domainJ sub coupling manager using its local indices
    template<std::size_t i, std::size_t j, class Apply>
    decltype(auto) subApply(Dune::index_constant<i> domainI,
                            Dune::index_constant<j> domainJ,
                            Apply&& apply)
    {
        constexpr auto localIndices = globalToLocal_(domainI, domainJ);
        return apply(subCouplingManager(domainI, domainJ), localIndices.first, localIndices.second);
    }

    //! apply a function to the domainI-domainJ sub coupling manager using its local indices
    template<std::size_t i, std::size_t j, class Apply>
    decltype(auto) subApply(Dune::index_constant<i> domainI,
                            Dune::index_constant<j> domainJ,
                            const Apply& apply) const
    {
        constexpr auto localIndices = globalToLocal_(domainI, domainJ);
        return apply(subCouplingManager(domainI, domainJ), localIndices.first, localIndices.second);
    }

    //! apply a function to a sub coupling manager containing this domain
    template<std::size_t i, class Apply>
    decltype(auto) subApply(Dune::index_constant<i> domainI, Apply&& apply)
    {
        constexpr auto dm = CouplingMap::coupledDomains(domainI);
        static_assert(std::tuple_size_v<std::decay_t<decltype(dm)>> != 0, "subApply on uncoupled domain!");
        constexpr auto domainJ = std::get<0>(dm);
        constexpr auto localIndices = globalToLocal_(domainI, domainJ);
        return apply(subCouplingManager(domainI, domainJ), localIndices.first);
    }

    //! apply a function to a sub coupling manager containing this domain
    template<std::size_t i, class Apply>
    decltype(auto) subApply(Dune::index_constant<i> domainI, const Apply& apply) const
    {
        constexpr auto dm = CouplingMap::coupledDomains(domainI);
        static_assert(std::tuple_size_v<std::decay_t<decltype(dm)>> != 0, "subApply on uncoupled domain!");
        constexpr auto domainJ = std::get<0>(dm);
        constexpr auto localIndices = globalToLocal_(domainI, domainJ);
        return apply(subCouplingManager(domainI, domainJ), localIndices.first);
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
        if constexpr (isCoupled_<CouplingMap>(domainI, domainJ))
            return subApply(domainI, domainJ, [&](const auto& cm, auto&& ii, auto&& jj) -> const auto& {
                return cm.couplingStencil(ii, entity, jj);
            });
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
        if constexpr (i == j || isCoupled_<CouplingMap>(domainI, domainJ))
            return subApply(domainI, domainJ, [&](const auto& cm, auto&& ii, auto&& jj) -> decltype(auto) {
                return cm.evalCouplingResidual(ii, scvfI, localAssemblerI, jj, dofIdxGlobalJ);
            });
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
        if constexpr (i == j || isCoupled_<CouplingMap>(domainI, domainJ))
            return subApply(domainI, domainJ, [&](const auto& cm, auto&& ii, auto&& jj) -> decltype(auto) {
                return cm.evalCouplingResidual(ii, localAssemblerI, jj, dofIdxGlobalJ);
            });
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
        if constexpr (i == j || isCoupled_<CouplingMap>(domainI, domainJ))
            return subApply(domainI, domainJ, [&](const auto& cm, auto&& ii, auto&& jj) -> decltype(auto) {
                return cm.evalCouplingResidual(ii, localAssemblerI, scvI, jj, dofIdxGlobalJ);
            });
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
            if constexpr (isCoupled_<CouplingMap>(domainI, domainJ))
                subApply(domainI, domainJ, [&](auto& cm, auto&& ii, auto&& jj){
                    cm.updateCouplingContext(ii, localAssemblerI, jj, dofIdxGlobalJ, priVars, pvIdxJ);
                });
        }
        else
        {
            // for i == j, we need to call all relevant managers where domainI is involved and make sure that the
            // context is updated with respect to its own domain (I-I coupling context)
            using namespace Dune::Hybrid;
            forEach(CouplingMap::coupledDomains(domainI), [&](const auto domainJ){
                static_assert(domainI != domainJ);
                subApply(domainI, domainJ, [&](auto& cm, auto&& ii, auto&& jj){
                    cm.updateCouplingContext(ii, localAssemblerI, ii, dofIdxGlobalJ, priVars, pvIdxJ);
                });
            });
        }
    }

    //! Bind the coupling context for a low dim element TODO remove Assembler
    template<std::size_t i, class Assembler = int>
    void bindCouplingContext(Dune::index_constant<i> domainI, const Element<i>& element, const Assembler& assembler = 0)
    {
        // for the coupling blocks
        using namespace Dune::Hybrid;
        forEach(CouplingMap::coupledDomains(domainI), [&](const auto domainJ){
            subApply(domainI, domainJ, [&](auto& cm, auto&& ii, auto&& jj) -> void {
                cm.bindCouplingContext(ii, element, assembler);
            });
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
        return subApply(domainI, [&](const auto& cm, auto&& ii) -> decltype(auto) {
            return cm.numericEpsilon(ii, paramGroup);
        });
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
        forEach(CouplingMap::coupledDomains(domainI), [&](const auto domainJ){
            subApply(domainI, domainJ, [&](auto& cm, auto&& ii, auto&& jj){
                cm.updateCoupledVariables(ii, localAssemblerI, elemVolVars, elemFluxVarsCache);
            });
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
