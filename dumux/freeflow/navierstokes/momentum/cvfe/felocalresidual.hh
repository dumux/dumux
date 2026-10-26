// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \brief Helper functions for assembling FE-based local residuals
 */
#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_CVFE_FE_LOCAL_RESIDUAL_HELPER_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_CVFE_FE_LOCAL_RESIDUAL_HELPER_HH

#include <dune/geometry/quadraturerules.hh>
#include <dumux/common/typetraits/griddiscretization.hh>
#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/common/boundaryflag.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/fem/interpolationpointdata.hh>

#include "flux.hh"

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Helper class for evaluating FE-based local residuals
 */
template<class Scalar, class NumEqVector, class LocalBasis, class Extrusion>
class NavierStokesMomentumFELocalResidualTerms
{
    using RangeType = typename LocalBasis::Traits::RangeType;

public:
    /*!
     * \brief Add the storage of non-CV local dofs at one time level
     *
     * \param storage The element storage vector to add to
     * \param problem The problem to solve
     * \param elemDisc The finite-volume geometry of the element
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

            const auto& localBasis = elemDisc.feLocalBasis();
            std::vector<RangeType> integralShapeFunctions(localBasis.size(), RangeType(0.0));

            // We apply mass lumping such that we only need to calculate the integral of basis functions
            const auto& geometry = elemDisc.elementGeometry();
            const auto& element = elemDisc.element();
            using GlobalPosition = typename GridDiscretization_t<ElementDiscretization>::GlobalCoordinate;
            using FeIpData = FEInterpolationPointData<GlobalPosition, LocalBasis>;

            for (const auto& qpData : CVFE::quadratureRule(elemDisc, element))
            {
                const auto& ipData = qpData.ipData();
                // Obtain and store shape function values and gradients at the current quad point
                FeIpData feIpData(geometry, ipData.local(), ipData.global(), localBasis);

                // get density from the problem
                for (const auto& localDof : nonCVLocalDofs(elemDisc))
                    integralShapeFunctions[localDof.index()] += qpData.weight() * feIpData.shapeValue(localDof.index());
            }

            for (const auto& localDof : nonCVLocalDofs(elemDisc))
            {
                const auto localDofIdx = localDof.index();
                const auto& data = ipData(elemDisc, localDof);
                const auto density = problem.density(element, elemDisc, data, isPreviousTimeLevel);
                const auto momentum = density*elemVars[localDofIdx].velocity();

                for (int eqIdx = 0; eqIdx < NumEqVector::dimension; ++eqIdx)
                    storage[localDofIdx][eqIdx] += integralShapeFunctions[localDofIdx][0]*momentum[eqIdx];
            }
        }
    }

    /*!
     * \brief Add flux and source residual contribution for non-CV local dofs
     *
     * \param residual The element residual vector to add to
     * \param problem The problem to solve
     * \param elemDisc The finite-volume geometry of the element
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
                    DUNE_THROW(Dune::NotImplemented, "Point sources are not implemented for hybrid momentum schemes.");
            }

            using GridDiscretization = GridDiscretization_t<ElementDiscretization>;
            using FluxTerms = Experimental::CVFENavierStokesMomentumFluxAtIp<Scalar, GridDiscretization>;

            const auto& element = elemDisc.element();
            for (const auto& qpData : CVFE::quadratureRule(elemDisc, element))
            {
                const auto& ipData = qpData.ipData();
                const auto& ipCache = cache(elemVars, ipData);
                const auto& shapeValues = ipCache.shapeValues();

                auto fluxTerm = FluxTerms::diffusiveFluxTerm(problem, element, elemDisc, elemVars, ipData, ipCache);

                // the density enters both the advective and the gravity term
                const Scalar density = problem.density(element, elemDisc, ipData);

                // there is no face to upwind across, so the advective term is interpolated as well
                if (problem.enableInertiaTerms())
                {
                    const auto velocity = FluxTerms::velocity(elemDisc, elemVars, ipCache);
                    const auto momentum = [&](const auto&, const auto&, const auto&) { return density*velocity; };
                    const auto transportingVelocity = [&](const auto&, const auto&, const auto&) { return velocity; };
                    fluxTerm += FluxTerms::advectiveFluxTerm(elemDisc, elemVars, ipData, momentum, transportingVelocity);
                }

                // add gravity term rho*g (note that gravity might be zero in case it's disabled in the problem)
                auto source = problem.source(elemDisc, elemVars, ipData);
                source.axpy(density, problem.gravity());

                const auto pressureFluxTerm = FluxTerms::pressureFluxTerm(problem, element, elemDisc, ipData);

                // the weak form -F grad N_i - N_i q
                for (const auto& localDof : nonCVLocalDofs(elemDisc))
                {
                    const auto localDofIdx = localDof.index();
                    NumEqVector fluxAndSourceTerm(0.0);
                    fluxTerm.mmv(ipCache.gradN(localDofIdx), fluxAndSourceTerm);
                    pressureFluxTerm.mmv(ipCache.gradN(localDofIdx), fluxAndSourceTerm);
                    fluxAndSourceTerm.axpy(-shapeValues[localDofIdx][0], source);

                    for (int eqIdx = 0; eqIdx < NumEqVector::dimension; ++eqIdx)
                        residual[localDofIdx][eqIdx] += qpData.weight()*fluxAndSourceTerm[eqIdx];
                }
            }
        }
    }

};

/*!
 * \ingroup NavierStokesModel
 * \brief Element-wise calculation of the Navier-Stokes residual for models using FE discretizations
 */
template<class TypeTag>
class NavierStokesMomentumFELocalResidual
: public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridVariablesCache = typename GridVariables::GridVariablesCache;
    using ElementVariables = typename GridVariablesCache::LocalView;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridDiscretization::LocalView;
    using GridView = typename GridDiscretization::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    using Extrusion = Extrusion_t<GridDiscretization>;

    using LocalBasis = typename GridDiscretization::FeCache::FiniteElementType::Traits::LocalBasisType;
    using FeResidual = NavierStokesMomentumFELocalResidualTerms<Scalar, NumEqVector, LocalBasis, Extrusion>;

public:
    using ElementResidualVector = typename ParentType::ElementResidualVector;
    using ParentType::ParentType;

    void addToElementStorage(ElementResidualVector& storage,
                             const Problem& problem,
                             const Element& element,
                             const ElementDiscretization& elemDisc,
                             const ElementVariables& elemVars,
                             bool isPreviousTimeLevel) const
    {
        FeResidual::addStorageTerms(storage, problem, elemDisc, elemVars, isPreviousTimeLevel);
    }

    void addToElementFluxAndSourceResidual(ElementResidualVector& residual,
                                           const Problem& problem,
                                           const Element& element,
                                           const ElementDiscretization& elemDisc,
                                           const ElementVariables& elemVars) const
    {
        FeResidual::addFluxAndSourceTerms(
            residual, problem, elemDisc, elemVars
        );
    }
};

} // end namespace Dumux

#endif
