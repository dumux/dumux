// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \brief Pure-FE (continuous Galerkin) continuity residual for the one-phase Navier-Stokes
 *        mass subdomain.
 *
 * Implements the weak incompressible continuity: w*div(u)*dx
 */
#ifndef DUMUX_NAVIERSTOKES_MASS_ONEP_FE_LOCAL_RESIDUAL_HH
#define DUMUX_NAVIERSTOKES_MASS_ONEP_FE_LOCAL_RESIDUAL_HH

#include <vector>

#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Element-wise FE continuity residual for the one-phase mass subdomain
 */
template<class TypeTag>
class NavierStokesMassOnePFELocalResidual
: public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridVariablesCache = typename GridVariables::GridVariablesCache;
    using ElementVariables = typename GridVariablesCache::LocalView;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using ElementResidualVector = typename ParentType::ElementResidualVector;
    using ParentType::ParentType;

    /*!
     * \brief Adds the storage at one time level
     *
     * \f[ s_i = \int_E q_i \varrho \,dx \f]
     *
     * The term is mass lumped, so that only the integrals of the basis functions are needed.
     * The time discretization is applied by the caller, as for storageIntegral().
     * The storage drops out for an incompressible fluid.
     */
    void addToElementStorage(ElementResidualVector& storage,
                             const Problem& problem,
                             const Element& element,
                             const ElementDiscretization& elemDisc,
                             const ElementVariables& elemVars,
                             bool isPreviousTimeLevel) const
    {
        if constexpr (Detail::LocalDofs::hasNonCVLocalDofsInterface<ElementDiscretization>())
        {
            if (nonCVLocalDofs(elemDisc).empty())
                return;

            // mass lumping, so only the integrals of the basis functions are needed
            std::vector<Scalar> integralShapeFunctions(Detail::LocalDofs::numLocalDofs(elemDisc), 0.0);
            for (const auto& qpData : CVFE::quadratureRule(elemDisc, element))
            {
                const auto& shapeValues = cache(elemVars, qpData.ipData()).shapeValues();
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
     * \brief Adds the continuity residual
     *
     * \f[ r_i = \int_E q_i \left( \varrho\nabla\cdot\boldsymbol{u}
     *                             + \boldsymbol{u}\cdot\nabla\varrho - q \right) \,dx \f]
     *
     * The velocity and its divergence come from the coupling manager at each quadrature
     * point, the density and its gradient are interpolated from the local dofs. The density
     * gradient drops out for an incompressible fluid.
     */
    void addToElementFluxAndSourceResidual(ElementResidualVector& residual,
                                           const Problem& problem,
                                           const Element& element,
                                           const ElementDiscretization& elemDisc,
                                           const ElementVariables& elemVars) const
    {
        if constexpr (Detail::LocalDofs::hasNonCVLocalDofsInterface<ElementDiscretization>())
        {
            if (nonCVLocalDofs(elemDisc).empty())
                return;

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

                for (const auto& localDof : nonCVLocalDofs(elemDisc))
                {
                    const auto localDofIdx = localDof.index();
                    residual[localDofIdx][Indices::conti0EqIdx]
                        += qpData.weight() * shapeValues[localDofIdx]
                           * (divMassFlux - source[Indices::conti0EqIdx]);
                }
            }
        }
    }
};

} // end namespace Dumux

#endif
