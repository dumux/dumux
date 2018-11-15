// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup EmbeddedCoupling
 * \brief Extended source stencil helper class for coupling managers
 */

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_EXTENDEDSOURCESTENCIL_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_EXTENDEDSOURCESTENCIL_HH

#include <vector>
#include <dumux/common/properties.hh>

namespace Dumux {
namespace EmbeddedCoupling {

/*!
 * \ingroup MultiDomain
 * \ingroup EmbeddedCoupling
 * \brief A class managing an extended source stencil
 * \tparam CouplingManager the coupling manager type
 */
template<class CouplingManager>
class ExtendedSourceStencil
{
    using MDTraits = typename CouplingManager::MultiDomainTraits;
    using Scalar = typename MDTraits::Scalar;

    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomainTypeTag<id>;
    template<std::size_t id> using FVGridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::FVGridGeometry>;
    template<std::size_t id> using GridView = typename FVGridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;

    static constexpr auto bulkIdx = typename MDTraits::template DomainIdx<0>();
    static constexpr auto lowDimIdx = typename MDTraits::template DomainIdx<1>();

    template<std::size_t id>
    static constexpr bool isBox()
    { return FVGridGeometry<id>::discMethod == DiscretizationMethod::box; }
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

            if (isBox<domainI>())
            {
                for (int i = 0; i < element.subEntities(GridView<domainI>::dimension); ++i)
                    for (const auto globalJ : dofs)
                        pattern.add(couplingManager.problem(domainI).fvGridGeometry().vertexMapper().subIndex(element, i, GridView<domainI>::dimension), globalJ);
            }
            else
            {
                const auto globalI = couplingManager.problem(domainI).fvGridGeometry().elementMapper().index(element);
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
    template<std::size_t i, class LocalAssemblerI, class SolutionVector, class JacobianMatrixDiagBlock, class GridVariables>
    void evalAdditionalDomainDerivatives(CouplingManager& couplingManager,
                                         Dune::index_constant<i> domainI,
                                         const LocalAssemblerI& localAssemblerI,
                                         const SolutionVector& curSol,
                                         JacobianMatrixDiagBlock& A,
                                         GridVariables& gridVariables) const
    {
        constexpr auto numEq = std::decay_t<decltype(curSol[domainI][0])>::dimension;
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
            const auto origPriVars = curSol[domainI][dofIndex];

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
                NumericDifferentiation::partialDerivative(evalResiduals, curSol[domainI][dofIndex][pvIdx],
                                                          partialDerivs, origResidual, numDiffMethod);

                // update the global stiffness matrix with the current partial derivatives
                for (auto&& scvJ : scvs(localAssemblerI.fvGeometry()))
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

                // restore the original coupling context
                couplingManager.updateCouplingContext(domainI, localAssemblerI, domainI, dofIndex, origPriVars, pvIdx);
            }
        }
    }

    //! return a reference to the stencil
    typename CouplingManager::CouplingStencils& stencil()
    { return sourceStencils_; }

private:
    //! the extended source stencil for the bulk domain due to the source average
    const std::vector<std::size_t>& extendedSourceStencil_(const CouplingManager& couplingManager, Dune::index_constant<0> bulkDomain, const Element<0>& bulkElement) const
    {
        const auto bulkElementIdx = couplingManager.problem(bulkIdx).fvGridGeometry().elementMapper().index(bulkElement);
        if (sourceStencils_.count(bulkElementIdx))
            return sourceStencils_.at(bulkElementIdx);
        else
            return couplingManager.emptyStencil();
    }

    //! the extended source stencil for the low dim domain is empty
    const std::vector<std::size_t>& extendedSourceStencil_(const CouplingManager& couplingManager, Dune::index_constant<1> bulkDomain, const Element<1>& lowDimElement) const
    { return couplingManager.emptyStencil(); }

    //! the additional stencil for the kernel evaluations / circle averages
    typename CouplingManager::CouplingStencils sourceStencils_;
};

} // end namespace EmbeddedCoupling
} // end namespace Dumux

#endif
