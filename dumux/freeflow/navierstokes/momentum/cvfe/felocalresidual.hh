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
class NavierStokesMomentumFELocalResidual
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

            static const auto intOrder
                = getParamFromGroup<int>(problem.paramGroup(), "Assembly.FEIntegrationOrderStorage", 4);

            const auto& localBasis = fvGeometry.feLocalBasis();
            std::vector<RangeType> integralShapeFunctions(localBasis.size(), RangeType(0.0));

            // We apply mass lumping such that we only need to calculate the integral of basis functions
            const auto& geometry = fvGeometry.elementGeometry();
            const auto& element = fvGeometry.element();
            using GlobalPosition = typename FVElementGeometry::GridGeometry::GlobalCoordinate;
            using IpData = FEInterpolationPointData<GlobalPosition, LocalBasis>;
            const auto& quadRule = Dune::QuadratureRules<Scalar, FVElementGeometry::GridGeometry::GridView::dimension>::rule(geometry.type(), intOrder);
            for (const auto& quadPoint : quadRule)
            {
                const Scalar qWeight = quadPoint.weight()*Extrusion::integrationElement(geometry, quadPoint.position());
                // Obtain and store shape function values and gradients at the current quad point
                IpData ipData(geometry, quadPoint.position(), localBasis);

                // get density from the problem
                for (const auto& localDof : nonCVLocalDofs(fvGeometry))
                    integralShapeFunctions[localDof.index()] += ipData.shapeValue(localDof.index())*qWeight;
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
     * \param elemFluxVarsCache The element flux variables cache
     * \param elemBcTypes The element boundary types
     */
    template<class ResidualVector, class Problem, class FVElementGeometry,
             class ElementVariables, class ElementFluxVariablesCache, class ElementBoundaryTypes>
    static void addFluxAndSourceTerms(ResidualVector& residual,
                                      const Problem& problem,
                                      const FVElementGeometry& fvGeometry,
                                      const ElementVariables& elemVars,
                                      const ElementFluxVariablesCache& elemFluxVarsCache,
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
            static const auto intOrder
                = getParamFromGroup<int>(problem.paramGroup(), "Assembly.FEIntegrationOrderFluxAndSource", 4);

            const auto& localBasis = fvGeometry.feLocalBasis();

            const auto& geometry = fvGeometry.elementGeometry();
            const auto& element = fvGeometry.element();
            using GlobalPosition = typename FVElementGeometry::GridGeometry::GlobalCoordinate;
            using IpData = FEInterpolationPointData<GlobalPosition, LocalBasis>;
            using FluxFunctionContext = NavierStokesMomentumFluxFunctionContext<Problem, FVElementGeometry, ElementVariables, IpData>;
            const auto& quadRule = Dune::QuadratureRules<Scalar, FVElementGeometry::GridGeometry::GridView::dimension>::rule(geometry.type(), intOrder);
            for (const auto& quadPoint : quadRule)
            {
                const Scalar qWeight = quadPoint.weight()*Extrusion::integrationElement(geometry, quadPoint.position());

                // Obtain and store shape function values and gradients at the current quad point
                IpData ipData(geometry, quadPoint.position(), localBasis);
                FluxFunctionContext context(problem, fvGeometry, elemVars, ipData);
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
                        fluxAndSourceTerm -= density*(v*ipData.gradN(localDofIdx))*v;

                    // add diffusion term
                    fluxAndSourceTerm += enableUnsymmetrizedVelocityGradient ?
                                            mu*mv(gradV, ipData.gradN(localDofIdx))
                                            : mu*mv(gradV + getTransposed(gradV), ipData.gradN(localDofIdx));

                    // add pressure term
                    fluxAndSourceTerm -= problem.pressure(element, fvGeometry, ipData) * ipData.gradN(localDofIdx);

                    // finally add source and Neumann term and add everything to residual
                    auto sourceAtIp = problem.source(fvGeometry, elemVars, ipData);
                    // add gravity term rho*g (note that gravity might be zero in case it's disabled in the problem)
                    sourceAtIp += density * problem.gravity();

                    for (int eqIdx = 0; eqIdx < NumEqVector::dimension; ++eqIdx)
                    {
                        fluxAndSourceTerm[eqIdx] -= ipData.shapeValue(localDofIdx) * sourceAtIp[eqIdx];
                        residual[localDofIdx][eqIdx] += qWeight*fluxAndSourceTerm[eqIdx];
                    }
                }
            }

            if (elemBcTypes.hasNeumann())
                addBoundaryFluxes(residual, problem, fvGeometry, elemVars, elemFluxVarsCache, elemBcTypes);
        }
    }

    /*!
     * \brief Evaluate Neumann boundary contributions
     *
     * \param residual The element residual vector to add to
     * \param problem The problem to solve
     * \param fvGeometry The finite-volume geometry of the element
     * \param elemVars The variables for all local dofs of the element
     * \param elemFluxVarsCache The element flux variables cache
     * \param elemBcTypes The element boundary types
     */
    template<class ResidualVector, class Problem, class FVElementGeometry,
             class ElementVariables, class ElementFluxVariablesCache, class ElementBoundaryTypes>
    static void addBoundaryFluxes(ResidualVector& residual,
                                  const Problem& problem,
                                  const FVElementGeometry& fvGeometry,
                                  const ElementVariables& elemVars,
                                  const ElementFluxVariablesCache& elemFluxVarsCache,
                                  const ElementBoundaryTypes& elemBcTypes)
    {
        ResidualVector flux(0.0);

        static const auto intOrder
            = getParamFromGroup<int>(problem.paramGroup(), "Assembly.FEIntegrationOrderBoundary", 4);

        const auto& localBasis = fvGeometry.feLocalBasis();

        const auto& geometry = fvGeometry.elementGeometry();
        using GlobalPosition = typename FVElementGeometry::GridGeometry::GlobalCoordinate;
        using GridView = typename FVElementGeometry::GridGeometry::GridView;
        using BoundaryFlag = Dumux::BoundaryFlag<typename GridView::Grid>;
        using FaceIpData = FEFaceInterpolationPointData<GlobalPosition, LocalBasis, BoundaryFlag>;

        const auto& element = fvGeometry.element();
        for (const auto& intersection : intersections(fvGeometry.gridGeometry().gridView(), element))
        {
            if (!intersection.boundary())
                continue;

            const auto& bcTypes = elemBcTypes.get(fvGeometry, intersection);
            if (!bcTypes.hasNeumann())
                continue;

            // select quadrature rule for intersection faces (dim-1)
            auto isGeometry = intersection.geometry();
            constexpr auto faceDim = FVElementGeometry::GridGeometry::GridView::dimension - 1;
            const auto& faceRule = Dune::QuadratureRules<Scalar, faceDim>::rule(isGeometry.type(), intOrder);
            for (const auto& quadPoint : faceRule)
            {
                // position of quadrature point in local coordinates of element
                auto local = geometry.local(isGeometry.global(quadPoint.position()));

                // get quadrature rule weight for intersection
                Scalar qWeight = quadPoint.weight() * Extrusion::integrationElement(isGeometry, quadPoint.position());
                FaceIpData faceIpData(geometry, local, localBasis, intersection.centerUnitOuterNormal(), BoundaryFlag{ intersection });

                for (const auto& localDof : nonCVLocalDofs(fvGeometry))
                {
                    const auto& boundaryFlux = qWeight*problem.boundaryFlux(fvGeometry, elemVars, elemFluxVarsCache, faceIpData);
                    // only add fluxes to equations for which Neumann is set
                    for (int eqIdx = 0; eqIdx < NumEqVector::dimension; ++eqIdx)
                        if (bcTypes.isNeumann(eqIdx))
                            flux[localDof.index()][eqIdx] += faceIpData.shapeValue(localDof.index()) * boundaryFlux[eqIdx];
                }
            }
        }
        residual += flux;
    }
};

} // end namespace Dumux

#endif
