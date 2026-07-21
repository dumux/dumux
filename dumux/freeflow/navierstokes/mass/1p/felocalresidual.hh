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
#ifndef DUMUX_NAVIERSTOKES_MASS_1P_FE_LOCAL_RESIDUAL_HH
#define DUMUX_NAVIERSTOKES_MASS_1P_FE_LOCAL_RESIDUAL_HH

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

public:
    using ElementResidualVector = typename ParentType::ElementResidualVector;
    using ParentType::ParentType;

    //! incompressible, stationary continuity has no storage term
    void addToElementStorageResidual(ElementResidualVector& residual,
                                     const Problem& problem,
                                     const Element& element,
                                     const ElementDiscretization& elemDisc,
                                     const ElementVariables& prevElemVars,
                                     const ElementVariables& curElemVars) const
    {}

    //! r_i = \int q_i (div u) dx, with div u from the coupling manager at each quadrature point
    void addToElementFluxAndSourceResidual(ElementResidualVector& residual,
                                           const Problem& problem,
                                           const Element& element,
                                           const ElementDiscretization& fvGeometry,
                                           const ElementVariables& elemVars) const
    {
        if constexpr (Detail::LocalDofs::hasNonCVLocalDofsInterface<ElementDiscretization>())
        {
            if (nonCVLocalDofs(fvGeometry).empty())
                return;

            for (const auto& qpData : CVFE::quadratureRule(fvGeometry, element))
            {
                const auto& ipData = qpData.ipData();
                const auto& ipCache = cache(elemVars, ipData);
                const auto& shapeValues = ipCache.shapeValues();

                const Scalar divV = problem.velocityDivergence(fvGeometry, ipData);

                for (const auto& localDof : nonCVLocalDofs(fvGeometry))
                {
                    const auto localDofIdx = localDof.index();
                    residual[localDofIdx][Indices::conti0EqIdx]
                        += qpData.weight() * shapeValues[localDofIdx] * divV;
                }
            }
        }
    }
};

} // end namespace Dumux

#endif
