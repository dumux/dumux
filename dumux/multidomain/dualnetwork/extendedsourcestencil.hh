// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup DualNetworkCoupling
 * \ingroup PoreNetworkModels
 * \brief Extended source stencil helper class for coupling managers
 */

#ifndef DUMUX_DUAL_NETWORK_EXTENDEDSOURCESTENCIL_HH
#define DUMUX_DUAL_NETWORK_EXTENDEDSOURCESTENCIL_HH

#include <vector>

#include <dune/common/indices.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numericdifferentiation.hh>
#include <dumux/discretization/method.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup DualNetworkCoupling
 * \ingroup PoreNetworkModels
 * \brief A class managing an extended source stencil
 * \tparam CouplingManager the coupling manager type
 */
template<class CouplingManager>
class PNMHeatExtendedSourceStencil
{
    using MDTraits = typename CouplingManager::MultiDomainTraits;
    using Scalar = typename MDTraits::Scalar;

    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;

    static constexpr auto solidDomainIdx = CouplingManager::solidDomainIdx;
    static constexpr auto voidDomainIdx = CouplingManager::voidDomainIdx;

    template<std::size_t id>
    static constexpr bool isBox()
    { return GridGeometry<id>::discMethod == DiscretizationMethods::box; }
public:
    /*!
     * \brief extend the jacobian pattern of the diagonal block of domain i
     *        by those entries that are not already in the uncoupled pattern
     * \note Such additional dependencies can arise from the coupling, e.g. if a coupling source
     *       term depends on a non-local average of a quantity of the same domain
     */
    template<std::size_t id, class JacobianPattern>
    void extendJacobianPattern(const CouplingManager& couplingManager, Dune::index_constant<id> domainI, JacobianPattern& pattern) const
    {
        // add additional dof dependencies
        for (const auto& element : elements(couplingManager.gridView(domainI)))
        {
            const auto& dofs = extendedSourceStencil_(couplingManager, domainI, element);
            if constexpr (isBox<domainI>())
            {
                for (int i = 0; i < element.subEntities(GridView<domainI>::dimension); ++i)
                    for (const auto globalJ : dofs)
                        pattern.add(couplingManager.problem(domainI).gridGeometry().vertexMapper().subIndex(element, i, GridView<domainI>::dimension), globalJ);
            }
            else
            {
                const auto globalI = couplingManager.problem(domainI).gridGeometry().elementMapper().index(element);
                for (const auto globalJ : dofs)
                    pattern.add(globalI, globalJ);
            }
        }
    }

    /*!
     * \brief evaluate additional derivatives of the element residual of a domain with respect
     *        to dofs in the same domain that are not in the regular stencil (per default this is not the case)
     * \note Such additional dependencies can arise from the coupling, e.g. if a coupling source
     *       term depends on a non-local average of a quantity of the same domain
     * \note This is the same for box and cc
     */
    template<std::size_t i, class LocalAssemblerI, class SolutionVector, class JacobianMatrixDiagBlock, class GridVariables> //TODO edit function s.t. curSol from couplingmanager is used like in dumux/dumux/multidomain/embedded/extendedsourcestencil.hh (commit 3257876)
    void evalAdditionalDomainDerivatives(CouplingManager& couplingManager,
                                         Dune::index_constant<i> domainI,
                                         const LocalAssemblerI& localAssemblerI,
                                         const SolutionVector& curSol,
                                         JacobianMatrixDiagBlock& A,
                                         GridVariables& gridVariables) const
    {
        // const auto& curSolI = couplingManager.curSol(domainI); //TODO:curSol declared protected in .../dumux/dumux/multidomain/couplingmanager.hh
        const auto& curSolI = curSol;
        constexpr auto numEq = std::decay_t<decltype(curSolI[0])>::size();
        const auto& elementI = localAssemblerI.element();

        // only do something if we have an extended stencil
        if (extendedSourceStencil_(couplingManager, domainI, elementI).empty())
            return;

        // compute the undeflected residual (source only!)
        const auto origResidual = localAssemblerI.evalLocalSourceResidual(elementI);

        // compute derivate for all additional dofs in the circle stencil
        for (const auto dofIndex : extendedSourceStencil_(couplingManager, domainI, elementI))
        {
            // std::cout << "deflecting " << dofIndex << std::endl;
            auto partialDerivs = origResidual;
            const auto origPriVars = curSolI[dofIndex];

            // calculate derivatives w.r.t to the privars at the dof at hand
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
                // reset partial derivatives
                partialDerivs = 0.0;

                auto evalResiduals = [&](Scalar priVar)
                {
                    // update the coupling context (solution vector and recompute element residual)
                    auto priVars = origPriVars;
                    priVars[pvIdx] = priVar;
                    couplingManager.updateCouplingContext(domainI, localAssemblerI, domainI, dofIndex, priVars, pvIdx);
                    return localAssemblerI.evalLocalSourceResidual(elementI);
                };

                // derive the residuals numerically
                static const int numDiffMethod = getParam<int>("Assembly.NumericDifferenceMethod");
                NumericDifferentiation::partialDerivative(evalResiduals, curSolI[dofIndex][pvIdx],
                                                          partialDerivs, origResidual, numDiffMethod); //TODO: maybe change like in dumux/dumux/md/embedded/extendedsourcestencil.hh line 142

                // update the global stiffness matrix with the current partial derivatives
                for (const auto& scvJ : scvs(localAssemblerI.fvGeometry()))
                {
                    for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                    {
                        // A[i][col][eqIdx][pvIdx] is the rate of change of
                        // the residual of equation 'eqIdx' at dof 'i'
                        // depending on the primary variable 'pvIdx' at dof
                        // 'col'.
                        // std::cout << "Adding " << partialDerivs[scvJ.indexInElement()][eqIdx] << ", pvIdx " << pvIdx <<  std::endl;
                        A[scvJ.dofIndex()][dofIndex][eqIdx][pvIdx] += partialDerivs[scvJ.indexInElement()][eqIdx];
                    }
                }

                // restore the original coupling context
                couplingManager.updateCouplingContext(domainI, localAssemblerI, domainI, dofIndex, origPriVars, pvIdx);
            }
        }
    }

    //! clear the internal data
    void clear() { sourceStencils_.clear(); }

    //! return a reference to the stencil
    auto& stencil()
    { return sourceStencils_; }

private:
    //! the extended source stencil due to the source average (always empty for lowdim, but may be filled for bulk)
    template<std::size_t id>
    const auto& extendedSourceStencil_(const CouplingManager& couplingManager, Dune::index_constant<id> domainId, const Element<id>& element) const
    {
        if constexpr (domainId == voidDomainIdx)
        {
            const auto voidElementIdx = couplingManager.problem(voidDomainIdx).gridGeometry().elementMapper().index(element);
            if (sourceStencils_.count(voidElementIdx))
                return sourceStencils_.at(voidElementIdx);
        }

        return couplingManager.emptyStencil(domainId);
    }

    //! the additional stencil for the kernel evaluations / circle averages
    std::unordered_map<std::size_t, std::vector<std::size_t>> sourceStencils_;
};

} // end namespace Dumux::PoreNetwork

#endif
