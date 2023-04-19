// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedCoupling
 * \brief Extended source stencil helper class for coupling managers
 */

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_EXTENDEDSOURCESTENCIL_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_EXTENDEDSOURCESTENCIL_HH

#include <vector>

#include <dune/common/indices.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numericdifferentiation.hh>
#include <dumux/assembly/numericepsilon.hh>

#include <dumux/discretization/method.hh>

namespace Dumux::EmbeddedCoupling {

/*!
 * \ingroup EmbeddedCoupling
 * \brief A class managing an extended source stencil
 * \tparam CouplingManager the coupling manager type
 */
template<class CouplingManager>
class ExtendedSourceStencil
{
    using MDTraits = typename CouplingManager::MultiDomainTraits;
    using Scalar = typename MDTraits::Scalar;

    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;

    static constexpr auto bulkIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto lowDimIdx = typename MDTraits::template SubDomain<1>::Index();

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
    template<std::size_t i, class LocalAssemblerI, class JacobianMatrixDiagBlock, class GridVariables>
    void evalAdditionalDomainDerivatives(CouplingManager& couplingManager,
                                         Dune::index_constant<i> domainI,
                                         const LocalAssemblerI& localAssemblerI,
                                         JacobianMatrixDiagBlock& A,
                                         GridVariables& gridVariables) const
    {
        const auto& curSolI = couplingManager.curSol(domainI);
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
            auto partialDerivs = origResidual;
            const auto origPriVars = curSolI[dofIndex];
            auto priVars = origPriVars;

            // calculate derivatives w.r.t to the privars at the dof at hand
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
                // reset partial derivatives
                partialDerivs = 0.0;

                const auto evalResiduals = [&](const Scalar priVar)
                {
                    // update the coupling context (solution vector and recompute element residual)
                    priVars[pvIdx] = priVar;
                    couplingManager.updateCouplingContext(domainI, localAssemblerI, domainI, dofIndex, priVars, pvIdx);
                    return localAssemblerI.evalLocalSourceResidual(elementI);
                };

                // derive the residuals numerically
                static const NumericEpsilon<Scalar, numEq> eps_{localAssemblerI.problem().paramGroup()};
                static const int numDiffMethod = getParamFromGroup<int>(localAssemblerI.problem().paramGroup(), "Assembly.NumericDifferenceMethod");
                NumericDifferentiation::partialDerivative(
                    evalResiduals, priVars[pvIdx], partialDerivs, origResidual, eps_(priVars[pvIdx], pvIdx), numDiffMethod
                );

                // update the global stiffness matrix with the current partial derivatives
                for (const auto& scvJ : scvs(localAssemblerI.fvGeometry()))
                {
                    for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                    {
                        // A[i][col][eqIdx][pvIdx] is the rate of change of
                        // the residual of equation 'eqIdx' at dof 'i'
                        // depending on the primary variable 'pvIdx' at dof
                        // 'col'.
                        A[scvJ.dofIndex()][dofIndex][eqIdx][pvIdx] += partialDerivs[scvJ.indexInElement()][eqIdx];
                    }
                }

                // restore the current element solution
                priVars[pvIdx] = origPriVars[pvIdx];

                // restore the original coupling context
                couplingManager.updateCouplingContext(domainI, localAssemblerI, domainI, dofIndex, origPriVars, pvIdx);
            }
        }
    }

    //! clear the internal data
    void clear() { sourceStencils_.clear(); }

    //! return a reference to the stencil
    typename CouplingManager::template CouplingStencils<bulkIdx>& stencil()
    { return sourceStencils_; }

    //! return a const reference to the stencil
    const typename CouplingManager::template CouplingStencils<bulkIdx>& stencil() const
    { return sourceStencils_; }

private:
    //! the extended source stencil due to the source average (always empty for lowdim, but may be filled for bulk)
    template<std::size_t id>
    const auto& extendedSourceStencil_(const CouplingManager& couplingManager, Dune::index_constant<id> domainId, const Element<id>& element) const
    {
        if constexpr (domainId == bulkIdx)
        {
            const auto bulkElementIdx = couplingManager.problem(bulkIdx).gridGeometry().elementMapper().index(element);
            if (auto stencil = sourceStencils_.find(bulkElementIdx); stencil != sourceStencils_.end())
                return stencil->second;
        }

        return couplingManager.emptyStencil(domainId);
    }

    //! the additional stencil for the kernel evaluations / circle averages
    typename CouplingManager::template CouplingStencils<bulkIdx> sourceStencils_;
};

} // end namespace Dumux::EmbeddedCoupling

#endif
