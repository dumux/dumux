// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \brief Element-wise calculation of the Navier-Stokes mass residual for single-phase flow when using fe discretization.
 */
#ifndef DUMUX_NAVIERSTOKES_MASS_ONEP_FE_LOCAL_RESIDUAL_HH
#define DUMUX_NAVIERSTOKES_MASS_ONEP_FE_LOCAL_RESIDUAL_HH

#include <vector>

#include <dune/common/exceptions.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/typetraits/griddiscretization.hh>
#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Helper class for evaluating FE-based local residuals
 */
template<class Scalar, class Indices>
class NavierStokesMassOnePFELocalResidualTerms
{
public:
    /*!
     * \brief Add the storage of non-CV local dofs at one time level
     *
     * \param storage The element storage vector to add to
     * \param problem The problem to solve
     * \param elemDisc The element discretization
     * \param elemVars The variables for all local dofs of the element at the given time level
     * \param isPreviousTimeLevel If the variables belong to the previous time level
     */
    template<class ResidualVector, class Problem, class ElementDiscretization, class ElementVariables>
    static void addStorageTerms(ResidualVector& storage,
                                const Problem& problem,
                                const ElementDiscretization& elemDisc,
                                const ElementVariables& elemVars,
                                const bool isPreviousTimeLevel)
    {
        if constexpr (Detail::LocalDofs::hasNonCVLocalDofsInterface<ElementDiscretization>())
        {
            // Make sure we don't iterate over quadrature points if there are no hybrid dofs
            if (nonCVLocalDofs(elemDisc).empty())
                return;

            std::vector<Scalar> integralShapeFunctions(Detail::LocalDofs::numLocalDofs(elemDisc), 0.0);

            // We apply mass lumping such that we only need to calculate the integral of basis functions
            const auto& element = elemDisc.element();
            for (const auto& qpData : CVFE::quadratureRule(elemDisc, element))
            {
                const auto& ipData = qpData.ipData();
                const auto& shapeValues = cache(elemVars, ipData).shapeValues();
                for (const auto& localDof : nonCVLocalDofs(elemDisc))
                    integralShapeFunctions[localDof.index()] += qpData.weight()*shapeValues[localDof.index()][0];
            }

            for (const auto& localDof : nonCVLocalDofs(elemDisc))
            {
                const auto localDofIdx = localDof.index();
                storage[localDofIdx][Indices::conti0EqIdx]
                    += integralShapeFunctions[localDofIdx]*elemVars[localDofIdx].density();
            }
        }
    }

    /*!
     * \brief Add flux and source residual contribution for non-CV local dofs
     *
     * \param residual The element residual vector to add to
     * \param problem The problem to solve
     * \param elemDisc The element discretization
     * \param elemVars The variables for all local dofs of the element
     */
    template<class ResidualVector, class Problem, class ElementDiscretization, class ElementVariables>
    static void addFluxAndSourceTerms(ResidualVector& residual,
                                      const Problem& problem,
                                      const ElementDiscretization& elemDisc,
                                      const ElementVariables& elemVars)
    {
        if constexpr (Detail::LocalDofs::hasNonCVLocalDofsInterface<ElementDiscretization>())
        {
            // Make sure we don't iterate over quadrature points if there are no hybrid dofs
            if (nonCVLocalDofs(elemDisc).empty())
                return;

            if constexpr (requires { problem.pointSources(); })
            {
                if (!problem.pointSources().empty())
                    DUNE_THROW(Dune::NotImplemented, "Point sources are not implemented for hybrid mass schemes.");
            }

            using GlobalPosition = typename GridDiscretization_t<ElementDiscretization>::GlobalCoordinate;

            const auto& element = elemDisc.element();
            for (const auto& qpData : CVFE::quadratureRule(elemDisc, element))
            {
                const auto& ipData = qpData.ipData();
                const auto& ipCache = cache(elemVars, ipData);
                const auto& shapeValues = ipCache.shapeValues();

                const Scalar divV = problem.velocityDivergence(elemDisc, ipData);

                Scalar density = 0.0;
                GlobalPosition gradDensity(0.0);
                for (const auto& localDof : localDofs(elemDisc))
                {
                    const auto localDensity = elemVars[localDof.index()].density();
                    density += localDensity*shapeValues[localDof.index()][0];
                    gradDensity.axpy(localDensity, ipCache.gradN(localDof.index()));
                }

                const auto velocity = problem.velocity(elemDisc, ipData);
                const auto source = problem.source(elemDisc, elemVars, ipData);
                const Scalar divMassFlux = density*divV + (velocity*gradDensity);

                // the strong form N_i (div(rho v) - q), the continuity equation is not integrated by parts
                for (const auto& localDof : nonCVLocalDofs(elemDisc))
                {
                    const auto localDofIdx = localDof.index();
                    residual[localDofIdx][Indices::conti0EqIdx]
                        += qpData.weight() * shapeValues[localDofIdx][0] * (divMassFlux - source[Indices::conti0EqIdx]);
                }
            }
        }
    }
};

/*!
 * \ingroup NavierStokesModel
 * \brief Element-wise calculation of the Navier-Stokes mass residual for single-phase flow when using fe discretization.
 */
template<class TypeTag>
class NavierStokesMassOnePFELocalResidual
: public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridVariablesCache = typename GridVariables::GridVariablesCache;
    using ElementVariables = typename GridVariablesCache::LocalView;

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FeResidual = NavierStokesMassOnePFELocalResidualTerms<Scalar, Indices>;

public:
    using ElementResidualVector = typename ParentType::ElementResidualVector;
    using ParentType::ParentType;

    /*!
     * \brief Adds the storage at one time level
     * \note See NavierStokesMassOnePFELocalResidualTerms::addStorageTerms
     */
    void addToElementStorage(ElementResidualVector& storage,
                             const Problem& problem,
                             const Element& element,
                             const ElementDiscretization& elemDisc,
                             const ElementVariables& elemVars,
                             bool isPreviousTimeLevel) const
    {
        FeResidual::addStorageTerms(storage, problem, elemDisc, elemVars, isPreviousTimeLevel);
    }

    /*!
     * \brief Adds the continuity residual
     * \note See NavierStokesMassOnePFELocalResidualTerms::addFluxAndSourceTerms
     */
    void addToElementFluxAndSourceResidual(ElementResidualVector& residual,
                                           const Problem& problem,
                                           const Element& element,
                                           const ElementDiscretization& elemDisc,
                                           const ElementVariables& elemVars) const
    {
        FeResidual::addFluxAndSourceTerms(residual, problem, elemDisc, elemVars);
    }

    /*!
     * \brief Add boundary flux contributions
     *
     * The continuity equation is not integrated by parts, so the mass flux over the boundary
     * follows from the velocity and the boundary fluxes of the problem do not enter the residual.
     */
    void addBoundaryFluxIntegral(ElementResidualVector& residual,
                                 const Problem& problem,
                                 const ElementDiscretization& elemDisc,
                                 const ElementVariables& elemVars) const
    {}
};

} // end namespace Dumux

#endif
