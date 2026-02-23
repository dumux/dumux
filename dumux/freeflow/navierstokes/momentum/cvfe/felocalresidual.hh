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
#include <dumux/common/concepts/datacache_.hh>
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
     * \brief Add storage residual contribution for non-CV local dofs
     *
     * \param residual The element residual vector to add to
     * \param problem The problem to solve
     * \param fvGeometry The finite-volume geometry of the element
     * \param prevElemVars The variables for all local dofs of the element at the previous time level
     * \param curElemVars The variables for all local dofs of the element at the current  time level
     * \param timeStepSize The current time step size
     */
    template<class ResidualVector, class Problem, class FVElementGeometry, class ElementVariables>
    static void addStorageTerms(ResidualVector& residual,
                                const Problem& problem,
                                const FVElementGeometry& fvGeometry,
                                const ElementVariables& prevElemVars,
                                const ElementVariables& curElemVars,
                                const Scalar timeStepSize)
    {
        if constexpr (Detail::LocalDofs::hasNonCVLocalDofsInterface<FVElementGeometry>())
        {
            // Make sure we don't iterate over quadrature points if there are no hybrid dofs
            if (nonCVLocalDofs(fvGeometry).empty())
                return;

            const auto& localBasis = fvGeometry.feLocalBasis();
            std::vector<RangeType> integralShapeFunctions(localBasis.size(), RangeType(0.0));

            // We apply mass lumping such that we only need to calculate the integral of basis functions
            const auto& geometry = fvGeometry.elementGeometry();
            const auto& element = fvGeometry.element();
            using GlobalPosition = typename FVElementGeometry::GridGeometry::GlobalCoordinate;
            using FeIpData = FEInterpolationPointData<GlobalPosition, LocalBasis>;

            for (const auto& qpData : CVFE::quadratureRule(fvGeometry, element))
            {
                const auto& ipData = qpData.ipData();
                // Obtain and store shape function values and gradients at the current quad point
                FeIpData feIpData(geometry, ipData.local(), ipData.global(), localBasis);

                // get density from the problem
                for (const auto& localDof : nonCVLocalDofs(fvGeometry))
                    integralShapeFunctions[localDof.index()] += qpData.weight() * feIpData.shapeValue(localDof.index());
            }

            for (const auto& localDof : nonCVLocalDofs(fvGeometry))
            {
                const auto localDofIdx = localDof.index();
                const auto& data = ipData(fvGeometry, localDof);
                const auto curDensity = problem.density(element, fvGeometry, data, false);
                const auto prevDensity = problem.density(element, fvGeometry, data, true);
                const auto curVelocity = curElemVars[localDofIdx].velocity();
                const auto prevVelocity = prevElemVars[localDofIdx].velocity();
                auto timeDeriv = (curDensity*curVelocity - prevDensity*prevVelocity);
                timeDeriv /= timeStepSize;

                // add storage to residual
                for (int eqIdx = 0; eqIdx < NumEqVector::dimension; ++eqIdx)
                    residual[localDofIdx][eqIdx] += integralShapeFunctions[localDofIdx]*timeDeriv[eqIdx];
            }
        }
    }

    /*!
     * \brief Add flux and source residual contribution for non-CV local dofs
     *
     * \param residual The element residual vector to add to
     * \param problem The problem to solve
     * \param fvGeometry The finite-volume geometry of the element
     * \param elemVars The variables for all local dofs of the element
     * \param elemDataCache The element data cache
     * \param elemBcTypes The element boundary types
     */
    template<class ResidualVector, class Problem, class FVElementGeometry,
             class ElementVariables, class ElementDataCache, class ElementBoundaryTypes>
    static void addFluxAndSourceTerms(ResidualVector& residual,
                                      const Problem& problem,
                                      const FVElementGeometry& fvGeometry,
                                      const ElementVariables& elemVars,
                                      const ElementDataCache& elemDataCache,
                                      const ElementBoundaryTypes& elemBcTypes)
    {
        if constexpr (Detail::LocalDofs::hasNonCVLocalDofsInterface<FVElementGeometry>())
        {
            // Make sure we don't iterate over quadrature points if there are no hybrid dofs
            if (nonCVLocalDofs(fvGeometry).empty())
                return;

            if (!problem.pointSourceMap().empty())
                DUNE_THROW(Dune::NotImplemented, "Point sources are not implemented for hybrid momentum schemes.");

            static const bool enableUnsymmetrizedVelocityGradient
                = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.EnableUnsymmetrizedVelocityGradient", false);

            const auto& element = fvGeometry.element();
            using FluxVariablesCache = Concept::DataCache_t<ElementDataCache>;
            using FluxFunctionContext = NavierStokesMomentumFluxFunctionContext<Problem, FVElementGeometry, ElementVariables, FluxVariablesCache>;
            for (const auto& qpData : CVFE::quadratureRule(fvGeometry, element))
            {
                const auto& ipData = qpData.ipData();
                // Obtain and store shape function values and gradients at the current quad point
                const auto& cache = elemDataCache[ipData];
                FluxFunctionContext context(problem, fvGeometry, elemVars, cache);
                const auto& v = context.velocity();
                const auto& gradV = context.gradVelocity();

                // get viscosity from the problem
                const Scalar mu = problem.effectiveViscosity(element, fvGeometry, ipData);
                // get density from the problem
                const Scalar density = problem.density(element, fvGeometry, ipData);

                for (const auto& localDof : nonCVLocalDofs(fvGeometry))
                {
                    const auto localDofIdx = localDof.index();
                    NumEqVector fluxAndSourceTerm(0.0);
                    // add advection term
                    if (problem.enableInertiaTerms())
                        fluxAndSourceTerm -= density*(v*cache.gradN(localDofIdx))*v;

                    // add diffusion term
                    fluxAndSourceTerm += enableUnsymmetrizedVelocityGradient ?
                                            mu*mv(gradV, cache.gradN(localDofIdx))
                                            : mu*mv(gradV + getTransposed(gradV), cache.gradN(localDofIdx));

                    // add pressure term
                    fluxAndSourceTerm -= problem.pressure(element, fvGeometry, ipData) * cache.gradN(localDofIdx);

                    // finally add source and Neumann term and add everything to residual
                    auto sourceAtIp = problem.source(fvGeometry, elemVars, ipData);
                    // add gravity term rho*g (note that gravity might be zero in case it's disabled in the problem)
                    sourceAtIp += density * problem.gravity();

                    const auto& shapeValues = cache.shapeValues();
                    for (int eqIdx = 0; eqIdx < NumEqVector::dimension; ++eqIdx)
                    {
                        fluxAndSourceTerm[eqIdx] -= shapeValues[localDofIdx] * sourceAtIp[eqIdx];
                        residual[localDofIdx][eqIdx] += qpData.weight()*fluxAndSourceTerm[eqIdx];
                    }
                }
            }

            if (elemBcTypes.hasNeumann())
                addBoundaryFluxes(residual, problem, fvGeometry, elemVars, elemDataCache, elemBcTypes);
        }
    }

    /*!
     * \brief Evaluate Neumann boundary contributions
     *
     * \param residual The element residual vector to add to
     * \param problem The problem to solve
     * \param fvGeometry The finite-volume geometry of the element
     * \param elemVars The variables for all local dofs of the element
     * \param elemDataCache The element data cache
     * \param elemBcTypes The element boundary types
     */
    template<class ResidualVector, class Problem, class FVElementGeometry,
             class ElementVariables, class ElementDataCache, class ElementBoundaryTypes>
    static void addBoundaryFluxes(ResidualVector& residual,
                                  const Problem& problem,
                                  const FVElementGeometry& fvGeometry,
                                  const ElementVariables& elemVars,
                                  const ElementDataCache& elemDataCache,
                                  const ElementBoundaryTypes& elemBcTypes)
    {
        ResidualVector flux(0.0);

        const auto& element = fvGeometry.element();
        for (const auto& intersection : intersections(fvGeometry.gridGeometry().gridView(), element))
        {
            if (!intersection.boundary())
                continue;

            const auto& bcTypes = elemBcTypes.get(fvGeometry, intersection);
            if (!bcTypes.hasNeumann())
                continue;

            problem.addBoundaryFluxIntegrals(flux, fvGeometry, elemVars, elemDataCache, intersection, bcTypes);
        }
        residual += flux;
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

    using GridVariablesCache = typename GridVariables::GridVolumeVariables;
    using ElementVariables = typename GridVariablesCache::LocalView;

    using GridDataCache = Concept::GridCache_t<GridVariables>;
    using ElementDataCache = typename GridDataCache::LocalView;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    using Extrusion = Extrusion_t<GridGeometry>;

    using LocalBasis = typename GridGeometry::FeCache::FiniteElementType::Traits::LocalBasisType;
    using FeResidual = NavierStokesMomentumFELocalResidualTerms<Scalar, NumEqVector, LocalBasis, Extrusion>;

public:
    //! Use the parent type's constructor
    using ElementResidualVector = typename ParentType::ElementResidualVector;
    using ParentType::ParentType;

    void addToElementStorageResidual(ElementResidualVector& residual,
                                     const Problem& problem,
                                     const Element& element,
                                     const FVElementGeometry& fvGeometry,
                                     const ElementVariables& prevElemVolVars,
                                     const ElementVariables& curElemVolVars) const
    {
        FeResidual::addStorageTerms(
            residual, problem, fvGeometry, prevElemVolVars, curElemVolVars, this->timeLoop().timeStepSize()
        );
    }

    void addToElementFluxAndSourceResidual(ElementResidualVector& residual,
                                           const Problem& problem,
                                           const Element& element,
                                           const FVElementGeometry& fvGeometry,
                                           const ElementVariables& elemVars,
                                           const ElementDataCache& elemDataCache,
                                           const ElementBoundaryTypes &elemBcTypes) const
    {
        FeResidual::addFluxAndSourceTerms(
            residual, problem, fvGeometry, elemVars, elemDataCache, elemBcTypes
        );
    }
};

} // end namespace Dumux

#endif
