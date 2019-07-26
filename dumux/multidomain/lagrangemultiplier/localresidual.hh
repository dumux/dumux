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
 * \brief Local residual for lagrange multipliers to be
 *        used within the multidomain framework.
 */
#ifndef DUMUX_MULTIDOMAIN_LAGRANGE_MULTIPLIER_LOCAL_RESIDUAL_HH
#define DUMUX_MULTIDOMAIN_LAGRANGE_MULTIPLIER_LOCAL_RESIDUAL_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/exceptions.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/discretization/fem/ipdata.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \brief Local residual for lagrange multipliers to be
 *        used within the multidomain framework.
 * \todo TODO: Extend docu on the interface expected in the problem
 */
template<class TypeTag>
class LagrangeMultiplierLocalResidual : public FELocalResidual<TypeTag>
{
    using ParentType = FELocalResidual<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using AnsatzSpaceBasis = typename GridGeometry::AnsatzSpaceBasis;
    using FEElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    static_assert(GridGeometry::discMethod == DiscretizationMethod::fem,
        "Lagrange multiplier can currently only used with fem discretization.");

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using SecondaryVariables = typename GridVariables::SecondaryVariables;
    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();

public:
    //! Pull up parents' constructor
    using ParentType::ParentType;

    //! the container storing the residual on all dofs of an element
    using typename ParentType::ElementResidualVector;

    /*!
     * \brief Compute the element local residual with
     *        a given local view on the trial space.
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param
     * \todo TODO Doc me.
     */
    template< class ElementSolution, class TrialSpaceBasisLocalView >
    ElementResidualVector eval(const Element& element,
                               const FEElementGeometry& feGeometry,
                               const ElementSolution& elemSol,
                               const TrialSpaceBasisLocalView& trialLocalView) const
    {
        using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

        using AnsatzLocalBasis = typename AnsatzSpaceBasis::LocalView::Tree::FiniteElement::Traits::LocalBasisType;
        using TrialLocalBasis = typename TrialSpaceBasisLocalView::Tree::FiniteElement::Traits::LocalBasisType;

        using IpDataAnsatz = FEIntegrationPointData<GlobalPosition, AnsatzLocalBasis>;
        using IpDataTrial = FEIntegrationPointData<GlobalPosition, TrialLocalBasis>;

        const auto& ansatzLocalView = feGeometry.feBasisLocalView();
        const auto& ansatzLocalBasis = ansatzLocalView.tree().finiteElement().localBasis();

        const auto& trialLocalBasis = trialLocalView.tree().finiteElement().localBasis();
        const auto trialNumLocalDofs = trialLocalBasis.size();

        ElementResidualVector residual(trialNumLocalDofs);

        const auto& geometry = element.geometry();
        const auto& quadRule = Dune::QuadratureRules<Scalar, GridView::dimension>::rule(geometry.type(), this->integrationOrder());
        for (const auto& quadPoint : quadRule)
        {
            // Obtain and store shape function values and gradients at the current quad point
            IpDataAnsatz ipDataAnsatz(geometry, quadPoint.position(), ansatzLocalBasis);
            IpDataTrial ipDataTrial(geometry, quadPoint.position(), trialLocalBasis);

            // calculate secondary variables for the previous and the current solution at the ip
            SecondaryVariables secVars;
            secVars.update(elemSol, this->problem(), element, ipDataAnsatz);

            // evaluate equation for lagrange multiplier
            const auto res = this->problem().evalConstraint(element, feGeometry, elemSol, ipDataAnsatz, secVars);

            // add entries to residual vector
            Scalar qWeight = quadPoint.weight()*geometry.integrationElement(quadPoint.position());
            for (unsigned int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                for (unsigned int i = 0; i < trialNumLocalDofs; ++i)
                    residual[i][eqIdx] += qWeight
                                          *secVars.extrusionFactor()
                                          *ipDataTrial.shapeValue(i)
                                          *res[eqIdx];
        }

        // add contribution from neumann segments
        residual += this->evalNeumannSegments_(element, feGeometry, elemSol, trialLocalView);

        return residual;
    }

    /*!
     * \brief Compute the element local residual for stationary problems
     *        with ansatz and trial space being identical.
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param
     * \todo TODO Doc me.
     */
    template<class ElementSolution>
    ElementResidualVector eval(const Element& element,
                               const FEElementGeometry& feGeometry,
                               const ElementSolution& elemSol) const
    {
        using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
        using AnsatzLocalBasis = typename AnsatzSpaceBasis::LocalView::Tree::FiniteElement::Traits::LocalBasisType;
        using IpDataAnsatz = FEIntegrationPointData<GlobalPosition, AnsatzLocalBasis>;

        const auto& ansatzLocalView = feGeometry.feBasisLocalView();
        const auto& ansatzLocalBasis = ansatzLocalView.tree().finiteElement().localBasis();
        const auto ansatzNumLocalDofs = ansatzLocalBasis.size();

        ElementResidualVector residual(ansatzNumLocalDofs);

        const auto& geometry = element.geometry();
        const auto& quadRule = Dune::QuadratureRules<Scalar, GridView::dimension>::rule(geometry.type(), this->integrationOrder());
        for (const auto& quadPoint : quadRule)
        {
            // Obtain and store shape function values and gradients at the current quad point
            IpDataAnsatz ipDataAnsatz(geometry, quadPoint.position(), ansatzLocalBasis);

            // calculate secondary variables for the previous and the current solution at the ip
            SecondaryVariables secVars;
            secVars.update(elemSol, this->problem(), element, ipDataAnsatz);

            // evaluate source term contribution
            const auto res = this->problem().evalConstraint(element, feGeometry, elemSol, ipDataAnsatz, secVars);

            // add entries to residual vector
            Scalar qWeight = quadPoint.weight()*geometry.integrationElement(quadPoint.position());
            for (unsigned int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                for (unsigned int i = 0; i < ansatzNumLocalDofs; ++i)
                    residual[i][eqIdx] += qWeight
                                          *secVars.extrusionFactor()
                                          *ipDataAnsatz.shapeValue(i)
                                          *res[eqIdx];
        }

        // add contribution from neumann segments
        residual += this->evalNeumannSegments_(element, feGeometry, elemSol);

        return residual;
    }
};

} // end namespace Dumux

#endif
