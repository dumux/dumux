// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \brief FE (Galerkin) weak-form contributions on the non-CV (P2 edge) dofs for the co-located
 *        momentum + Cahn-Hilliard CVFE model.
 *
 * On the hybrid PQ2 scheme the vertex dofs are control volumes (assembled by the box/FV
 * compute{Storage,Flux,Source}) while the edge dofs are FE-only and get their residual from these
 * helpers. This is a fork of NavierStokesMomentumFELocalResidualTerms that (a) reads the mixture
 * material rho(c)/nu(c) from the LOCAL volvars (correct monolithic Jacobian), and (b) adds the
 * Cahn-Hilliard weak form for the c-transport and mu-definition equations.
 *
 * Weak forms (test function N_i), sign-matched to the mass/2p CV residual so vertex and edge dofs
 * share one convention:
 *   momentum [0..dim-1]:  -rho (v.gradN_i) v  +  nu (gradV+gradV^T).gradN_i  -  p gradN_i
 *                         -  N_i (source_mom + rho g)          [storage: mass-lumped rho dv/dt]
 *   phase field [dim]:    -c (v.gradN_i)  +  M (gradMu.gradN_i)  [storage: mass-lumped dc/dt]
 *   chem. pot.  [dim+1]:   sigma~eps (gradC.gradN_i)  -  N_i (mu - energyScale f'(c))
 * where source_mom is the capillary force (well-balanced -c gradMu) from problem.source() and the
 * chem.pot. source (mu - energyScale f'(c)) is problem.source()[chemPotEqIdx].
 */
#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_CHNS_CVFE_FE_LOCAL_RESIDUAL_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_CHNS_CVFE_FE_LOCAL_RESIDUAL_HH

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/fem/interpolationpointdata.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

#include <dumux/freeflow/navierstokes/momentum/chns/cvfe/flux.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief FE weak-form terms on non-CV dofs for the co-located momentum + Cahn-Hilliard model.
 * \tparam Scalar the scalar type
 * \tparam NumEqVector the full (dim+2) per-dof residual vector type
 * \tparam MomVec the dim-sized momentum/velocity vector type
 * \tparam Indices the model indices (phaseFieldEqIdx, chemicalPotentialEqIdx, ...)
 * \tparam LocalBasis the P2 local basis
 * \tparam Extrusion the extrusion policy
 */
template<class Scalar, class NumEqVector, class MomVec, class Indices, class LocalBasis, class Extrusion>
class NavierStokesMomentumCHNSFELocalResidualTerms
{
    using RangeType = typename LocalBasis::Traits::RangeType;
    static constexpr int dim = MomVec::dimension;
    static constexpr int phaseFieldEqIdx = Indices::phaseFieldEqIdx;
    static constexpr int chemicalPotentialEqIdx = Indices::chemicalPotentialEqIdx;

public:
    /*!
     * \brief Element-local monotonicity limiter for the ADVECTED phase field.
     *
     * Returns local phase-field coefficients in which each vertex (P1) dof keeps its value and each
     * P2 edge (midpoint) dof's hierarchical bump d = c_m - 0.5(c_a + c_b) is scaled by
     *     alpha = min(1, 0.25 |c_a - c_b| / |d|)  in [0,1]
     * so the quadratic profile along that edge stays monotone (|d| <= 0.25|c_a-c_b| is the exact 1D
     * no-interior-extremum bound). alpha -> 1 where the interface is smooth/well-resolved (keep full
     * P2), alpha -> 0 where the P2 reconstruction overshoots (degrade that edge to P1). It reads ONLY
     * the element's own dofs (the edge's two endpoints + midpoint), so it adds NO Jacobian stencil
     * entries beyond the already-present intra-element coupling. Because the edge midpoint dof and its
     * endpoints are SHARED by both incident elements, cLim is identical from either side -> the
     * reconstruction is continuous across edges and, used only in the divergence-form advective flux
     * (tested against the full partition of unity), conserves int c exactly.
     */
    template<class FVElementGeometry, class ElementVariables>
    static void limitedPhaseFieldCoeffs(const FVElementGeometry& fvGeometry,
                                        const ElementVariables& elemVars,
                                        std::vector<Scalar>& cLim,
                                        const Scalar boundScale = 0.25)
    {
        const auto& fe = fvGeometry.gridGeometry().feCache().get(fvGeometry.element().type());
        const auto& lc = fe.localCoefficients();
        const auto numLocalDofs = fe.localBasis().size();
        cLim.resize(numLocalDofs);
        for (std::size_t j = 0; j < numLocalDofs; ++j)
            cLim[j] = elemVars[j].phaseField();

        using ctype = typename FVElementGeometry::GridGeometry::GridView::ctype;
        const auto refElem = Dune::referenceElement<ctype, dim>(fvGeometry.element().type());
        constexpr int edgeCodim = dim - 1; // P2 edge dofs: codim 1 in 2D, codim 2 in 3D

        // map vertex sub-entity index -> local basis index (codim == dim are the P1 vertex nodes)
        std::array<int, 8> vertexBasisIdx;
        vertexBasisIdx.fill(-1);
        for (std::size_t j = 0; j < numLocalDofs; ++j)
            if (lc.localKey(j).codim() == dim)
                vertexBasisIdx[lc.localKey(j).subEntity()] = int(j);

        for (std::size_t j = 0; j < numLocalDofs; ++j)
        {
            const auto& key = lc.localKey(j);
            if (key.codim() != edgeCodim) continue;
            const int vA = vertexBasisIdx[refElem.subEntity(key.subEntity(), edgeCodim, 0, dim)];
            const int vB = vertexBasisIdx[refElem.subEntity(key.subEntity(), edgeCodim, 1, dim)];
            if (vA < 0 || vB < 0) continue;
            const Scalar cA = elemVars[vA].phaseField();
            const Scalar cB = elemVars[vB].phaseField();
            const Scalar mid = 0.5*(cA + cB);
            const Scalar d = elemVars[j].phaseField() - mid;
            const Scalar bound = boundScale*std::abs(cB - cA);
            const Scalar alpha = std::min(Scalar(1.0), bound/(std::abs(d) + 1e-12));
            cLim[j] = mid + alpha*d;
        }
    }

    //! standard consistent-mass Galerkin storage for the momentum block and the phase-field
    //! transport equation on the non-CV (edge) dofs
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
            if (nonCVLocalDofs(fvGeometry).empty())
                return;

            const auto& localBasis = fvGeometry.feLocalBasis();
            const auto& geometry = fvGeometry.elementGeometry();
            const auto& element = fvGeometry.element();
            using GlobalPosition = typename FVElementGeometry::GridGeometry::GlobalCoordinate;
            using FeIpData = FEInterpolationPointData<GlobalPosition, LocalBasis>;

            // Convective-form momentum storage rho^{n+1}(u^{n+1}-u^n)/dt when enabled (see
            // FreeFlow.ConvectiveFormMomentumStorage in the base momentum residual).
            static const bool convectiveFormStorage = getParamFromGroup<bool>(
                problem.paramGroup(), "FreeFlow.ConvectiveFormMomentumStorage", false);
            // The FE (edge-dof) part is standard consistent-mass Galerkin: R_i = int N_i d/dt(.).
            // Row-sum lumping (R_i = (int N_i) d/dt(.)_i) is only kept for A/B comparison via
            // FreeFlow.LumpFEStorage; it is NOT standard Galerkin and is questionable for the
            // co-located Cahn-Hilliard phase field (it was inherited from the base momentum model).
            static const bool lump = getParamFromGroup<bool>(
                problem.paramGroup(), "FreeFlow.LumpFEStorage", false);
            // Eyre convex splitting of the double-well f'(c)=E(c^3-c): treat the concave -E c part
            // explicitly (-E c^n) for unconditional gradient stability. Implemented here (the storage
            // hook has both time levels) as the correction +E(c^{n+1}-c^n) added to the mu-equation,
            // which turns the flux/source term's implicit -E c^{n+1} into -E c^n. Default off.
            static const bool eyre = getParamFromGroup<bool>(
                problem.paramGroup(), "FreeFlow.EyreConvexSplitting", false);
            const Scalar energyScale = eyre ? problem.energyScale() : Scalar(0.0);
            // Semi-implicit FULL Taylor linearization of the double-well (Aland & Voigt 2011):
            // linearize f'(c^{n+1})=c^3-c about c^n so the CH sub-problem is LINEAR (unconditionally
            // solvable). The correction on the mu-row E*(-2(c^n)^3 + 3(c^n)^2 c^{n+1} - (c^{n+1})^3)
            // exactly cancels the source term's implicit E*(c^{n+1})^3 (SAME consistent-mass Galerkin
            // interpolation at the quadrature point), leaving the linear part. Use INSTEAD of Eyre
            // (which only linearizes the concave -Ec part, leaving c^3 nonlinear). Gated (default off).
            static const bool semiImplicitDW = getParamFromGroup<bool>(
                problem.paramGroup(), "FreeFlow.SemiImplicitDoubleWell", false);
            const Scalar dwScale = semiImplicitDW ? problem.energyScale() : Scalar(0.0);

            // per-dof storage time rates d/dt(rho u) and d/dt(c)
            const auto numLocalDofs = localBasis.size();
            std::vector<MomVec> momRate(numLocalDofs, MomVec(0.0));
            std::vector<Scalar> cRate(numLocalDofs, Scalar(0.0));
            std::vector<Scalar> cOldArr(numLocalDofs, Scalar(0.0)), cNewArr(numLocalDofs, Scalar(0.0));
            for (const auto& localDof : localDofs(fvGeometry))
            {
                const auto j = localDof.index();
                const auto curDensity = curElemVars[j].density();
                const auto prevDensity = convectiveFormStorage ? curDensity : prevElemVars[j].density();
                momRate[j] = (curDensity*curElemVars[j].velocity() - prevDensity*prevElemVars[j].velocity());
                momRate[j] /= timeStepSize;
                cOldArr[j] = prevElemVars[j].phaseField();
                cNewArr[j] = curElemVars[j].phaseField();
                cRate[j] = (cNewArr[j] - cOldArr[j])/timeStepSize;
            }

            for (const auto& qpData : CVFE::quadratureRule(fvGeometry, element))
            {
                const auto& ipData = qpData.ipData();
                FeIpData feIpData(geometry, ipData.local(), ipData.global(), localBasis);
                const auto w = qpData.weight();

                // consistent mass: interpolate the storage rate at the quadrature point
                MomVec momRateIp(0.0); Scalar cRateIp(0.0);
                if (!lump)
                    for (const auto& localDof : localDofs(fvGeometry))
                    {
                        const auto j = localDof.index();
                        const Scalar Nj = feIpData.shapeValue(j)[0];
                        momRateIp.axpy(Nj, momRate[j]);
                        cRateIp += Nj*cRate[j];
                    }

                // consistent-mass interpolation of c^n, c^{n+1} for the semi-implicit double-well
                // correction (must match the source term's interpolation for exact cubic cancellation)
                Scalar cOldIp(0.0), cNewIp(0.0);
                if (semiImplicitDW)
                    for (const auto& localDof : localDofs(fvGeometry))
                    {
                        const auto j = localDof.index();
                        const Scalar Nj = feIpData.shapeValue(j)[0];
                        cOldIp += Nj*cOldArr[j];
                        cNewIp += Nj*cNewArr[j];
                    }
                const Scalar dwCorr = semiImplicitDW
                    ? dwScale*(-2.0*cOldIp*cOldIp*cOldIp + 3.0*cOldIp*cOldIp*cNewIp - cNewIp*cNewIp*cNewIp)
                    : Scalar(0.0);

                for (const auto& localDof : nonCVLocalDofs(fvGeometry))
                {
                    const auto i = localDof.index();
                    const Scalar Ni = feIpData.shapeValue(i)[0];
                    const auto& mr = lump ? momRate[i] : momRateIp; // lumped: nodal rate at i
                    const Scalar cr = lump ? cRate[i]  : cRateIp;
                    for (int d = 0; d < dim; ++d)
                        residual[i][d] += w*Ni*mr[d];
                    residual[i][phaseFieldEqIdx] += w*Ni*cr;
                    // Eyre: +E (c^{n+1}-c^n) = E*cr*dt (same consistent/lumped mass as the storage)
                    if (eyre)
                        residual[i][chemicalPotentialEqIdx] += w*Ni*energyScale*cr*timeStepSize;
                    // Semi-implicit full linearization of the double-well (alternative to Eyre)
                    if (semiImplicitDW)
                        residual[i][chemicalPotentialEqIdx] += w*Ni*dwCorr;
                }
            }
        }
    }

    //! flux + source weak form for momentum and Cahn-Hilliard on non-CV dofs
    template<class ResidualVector, class Problem, class FVElementGeometry, class ElementVariables>
    static void addFluxAndSourceTerms(ResidualVector& residual,
                                      const Problem& problem,
                                      const FVElementGeometry& fvGeometry,
                                      const ElementVariables& elemVars)
    {
        if constexpr (Detail::LocalDofs::hasNonCVLocalDofsInterface<FVElementGeometry>())
        {
            if (nonCVLocalDofs(fvGeometry).empty())
                return;

            static const bool enableUnsymmetrizedVelocityGradient
                = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.EnableUnsymmetrizedVelocityGradient", false);
            static const bool degenerateMobility
                = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.DegenerateMobility", false);
            // Skew-symmetric (Temam) forms: the "plain" terms below (-rho(v.gradN)v for momentum,
            // -c(v.gradN) for the phase field) are the weak form of the CONSERVATIVE/divergence
            // advection (integrated by parts against the zero-flux boundary), matching computeFlux's
            // divergence-form CV flux. Averaging that with the advective-form weak term
            // N_i*(v.grad(.)) cancels half the spurious div(v)!=0 force at the discretely
            // non-mass-conservative density/phase-field interface -- see the CV-side correction in
            // localresidual.hh for the full rationale. Both default off, matching the CV/base models.
            static const bool skewMomentum
                = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.SkewSymmetricMomentumAdvection", false);
            static const bool skewPhaseField
                = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.SkewSymmetricPhaseFieldAdvection", false);

            const auto& element = fvGeometry.element();
            using Cache = typename ElementVariables::InterpolationPointData;
            using FluxFunctionContext = NavierStokesMomentumFluxFunctionContext<Problem, FVElementGeometry, ElementVariables, Cache>;
            using GlobalPosition = typename FVElementGeometry::GridGeometry::GlobalCoordinate;
            using SurfaceTensionForm = typename Problem::SurfaceTensionForm;

            const Scalar sigmaEps = problem.surfaceTension(); // = sigma~ eps (Laplacian coefficient)
            const Scalar energyScale = problem.energyScale(); // = sigma~ / eps (double-well scale)
            const Scalar refDensity = problem.referenceDensity(); // Boussinesq reference (heavy phase)
            const auto& gravity = problem.gravity();
            const auto stForm = problem.surfaceTensionForm();
            const bool prescribeMu = problem.prescribeChemicalPotential(); // diagnostic mu-pin
            // Abels-Garcke-Grün diffusive momentum flux (edge-dof weak form of the CV interfaceFlux);
            // variable-density consistency, only with inertia. Default off.
            static const bool enableInterfaceFlux
                = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.EnableInterfaceFlux", false);
            const bool enableAGG = enableInterfaceFlux && problem.enableInertiaTerms();
            const Scalar dRhodC = problem.mixtureDensityDerivative();
            // Quadrature test: integrate the cubic double-well N_i E c(c^2-1) (degree 8 on P2, badly
            // under-integrated by the order-4 element cache) with a manual Dune rule of this order
            // instead (evaluating shape functions directly). 0 = use the cache-based integration.
            static const int dwQuadOrder
                = getParamFromGroup<int>(problem.paramGroup(), "FreeFlow.DoubleWellQuadOrder", 0);
            const bool manualDW = dwQuadOrder > 0 && !prescribeMu;
            // Consistent-Galerkin chemical-potential equation on ALL dofs (vertex + edge) instead of the
            // hybrid box(vertex)/Galerkin(edge) split. The mu-definition mu = -sigma~eps*Lap(c) + E f'(c)
            // is algebraic & elliptic (no advection), so there is no monotonicity reason to box-integrate
            // it on the CV vertices. The split projects the (non-constant) discrete g=-sigma~eps*Lap(c)+Ef'(c)
            // two different ways on vertices (box) vs edges (Galerkin), so a constant mu is NOT an exact
            // discrete equilibrium -> a high-frequency vertex-vs-edge grad(mu) ripple that drives the
            // solenoidal (pressure-irremovable) parasitic capillary current -c grad(mu). With this on, the
            // box path drops the mu-row (see localresidual fluxIntegral/sourceIntegral) and it is assembled
            // here by consistent Galerkin for every local dof. Default off.
            static const bool consistentChemPot
                = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.ConsistentChemPotGalerkin", false);
            // Balanced-force normal projection of the wellBalanced capillary force: replace grad(mu) by
            // its component along the interface normal n=gradC/|gradC|, dropping the TANGENTIAL grad(mu)
            // that is the whole source of the force's discrete curl (grad(mu) x grad(c)) and hence of the
            // solenoidal parasitic current. Reg floor keeps 1/|gradC|^2 finite -> force ->0 in the bulk.
            static const bool normalProject
                = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.CapillaryNormalProjection", false);
            static const Scalar npRegEps
                = getParamFromGroup<Scalar>(problem.paramGroup(), "FreeFlow.CapillaryNormalProjectionRegEps", 1e-8);
            // pressure-robust capillary force is assembled separately (localresidual.hh
            // addPressureRobustCapillary_) via an H(div) test-function reconstruction; skip the raw
            // wellbalanced/potential volumetric force here when it is on.
            static const bool pressureRobustCap
                = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.PressureRobustCapillary", false);
            // Boussinesq buoyancy (rho(c)-rho_ref) g is likewise gradient-dominated at the density jump;
            // when on it is reconstructed into H(div) (addPressureRobustCapillary_) instead of tested
            // against the raw N_i below -> removes the density-contrast+gravity parasitic current.
            static const bool pressureRobustGrav
                = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.PressureRobustGravity", false);
            // grad-div stabilization gamma*(div u, div v): penalises the discretely non-solenoidal
            // velocity modes the P1 pressure cannot control, damping the pressure-robustness parasitic
            // current at the density/CH interface. A convergent, symmetric-positive-definite alternative
            // to the (non-converging) RT1 reconstruction; gamma ~ O(viscosity..1). Default off.
            static const Scalar gradDivGamma
                = getParamFromGroup<Scalar>(problem.paramGroup(), "FreeFlow.GradDivStabilization", 0.0);
            // Optional Peclet-limited artificial diffusion c_stab*h*|v| for the phase-field advection
            // (stabilizes the unupwinded central-Galerkin edge dofs). Simple conservative baseline;
            // see the P2->P1 limiter note. Default off.
            static const Scalar phaseFieldArtDiff
                = getParamFromGroup<Scalar>(problem.paramGroup(), "FreeFlow.PhaseFieldArtificialDiffusion", 0.0);
            const Scalar hElem = std::pow(element.geometry().volume(), 1.0/Scalar(dim));

            // Element-local P2->P1 monotonicity limiter for the phase-field ADVECTION (see
            // limitedPhaseFieldCoeffs); reconstructs a limited cHat used only in the advective flux.
            static const bool phaseFieldLimiter
                = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.PhaseFieldMonotonicityLimiter", false);
            static const Scalar phaseFieldLimiterBound
                = getParamFromGroup<Scalar>(problem.paramGroup(), "FreeFlow.PhaseFieldLimiterBound", 0.25);
            std::vector<Scalar> cLim;
            if (phaseFieldLimiter)
                limitedPhaseFieldCoeffs(fvGeometry, elemVars, cLim, phaseFieldLimiterBound);

            for (const auto& qpData : CVFE::quadratureRule(fvGeometry, element))
            {
                const auto& ipData = qpData.ipData();
                const auto& ipCache = cache(elemVars, ipData);
                FluxFunctionContext context(problem, fvGeometry, elemVars, ipCache);
                const auto& v = context.velocity();
                const auto& gradV = context.gradVelocity();
                const auto& shapeValues = ipCache.shapeValues();

                // interpolate CH quantities and mixture material at the quadrature point
                Scalar c(0.0), mu(0.0), density(0.0), viscosity(0.0);
                GlobalPosition gradC(0.0), gradMu(0.0);
                // limited reconstruction of the ADVECTED phase field (cHat == c if limiter off)
                Scalar cHat(0.0);
                GlobalPosition gradCHat(0.0);
                for (const auto& localDof : localDofs(fvGeometry))
                {
                    const auto j = localDof.index();
                    const auto& vv = elemVars[j];
                    c += shapeValues[j][0]*vv.phaseField();
                    mu += shapeValues[j][0]*vv.chemicalPotential();
                    density += shapeValues[j][0]*vv.density();
                    viscosity += shapeValues[j][0]*vv.effectiveViscosity();
                    gradC.axpy(vv.phaseField(), ipCache.gradN(j));
                    gradMu.axpy(vv.chemicalPotential(), ipCache.gradN(j));
                    const Scalar cAdv = phaseFieldLimiter ? cLim[j] : vv.phaseField();
                    cHat += shapeValues[j][0]*cAdv;
                    gradCHat.axpy(cAdv, ipCache.gradN(j));
                }

                // Cahn-Hilliard mobility M (optionally degenerate M0 (c^2-1)^2, clamped)
                Scalar mobility = problem.mobility();
                if (degenerateMobility)
                {
                    const Scalar g = std::max(0.0, 1.0 - c*c);
                    mobility *= g*g;
                }

                // double-well derivative f'(c) = c(c^2-1); dropped here when integrated manually below
                const Scalar dwSource = manualDW ? mu : (mu - energyScale*c*(c*c - 1.0));

                // velocity divergence at this qp (for grad-div stabilization)
                Scalar divU = 0.0;
                for (int d = 0; d < dim; ++d) divU += gradV[d][d];

                // grad-div stabilization on ALL velocity dofs (vertex CV + edge FE): gamma (div u) *
                // div(N_k e_d) = gamma (div u) dN_k/dx_d. Assembled FE-style over the element for every
                // velocity dof; the vertex dofs' momentum is otherwise finite-volume, but this added
                // stabilization must penalise their divergence too or the parasitic current re-grows in
                // the vertex modes (observed with edge-only: mu bounded ~6 steps then climbs again).
                if (gradDivGamma > 0.0)
                    for (const auto& ld : localDofs(fvGeometry))
                    {
                        const auto k = ld.index();
                        const auto gNk = ipCache.gradN(k);
                        for (int d = 0; d < dim; ++d)
                            residual[k][d] += qpData.weight()*gradDivGamma*divU*gNk[d];
                    }

                for (const auto& localDof : nonCVLocalDofs(fvGeometry))
                {
                    const auto i = localDof.index();
                    const auto gradNi = ipCache.gradN(i);
                    const auto Ni = shapeValues[i][0];

                    // ---- momentum block [0..dim-1] ----
                    // (grad-div stabilization for these edge dofs is added in the all-dofs loop above)
                    MomVec mom(0.0);
                    if (problem.enableInertiaTerms())
                    {
                        const Scalar vGradNi = v*gradNi; // reused below
                        if (skewMomentum)
                        {
                            // 0.5*[ N_i*rho*(v.gradV_d) - rho*(v.gradN_i)*v_d ]
                            for (int d = 0; d < dim; ++d)
                                mom[d] += 0.5*(Ni*density*(v*gradV[d]) - density*vGradNi*v[d]);
                        }
                        else
                            mom -= density*vGradNi*v; // advection -rho(v.gradN)v (conservative form)
                    }
                    // AGG diffusive momentum flux (weak form of the CV interfaceFlux): -(J.gradNi) v
                    // with J = -M (drho/dc) grad mu, i.e. momentum advected with the CH mass flux.
                    if (enableAGG)
                    {
                        const Scalar JgradNi = -mobility*dRhodC*(gradMu*gradNi);
                        for (int d = 0; d < dim; ++d)
                            mom[d] -= JgradNi*v[d];
                    }
                    mom += enableUnsymmetrizedVelocityGradient
                             ? viscosity*mv(gradV, gradNi)
                             : viscosity*mv(gradV + getTransposed(gradV), gradNi); // viscous
                    mom -= problem.pressure(element, fvGeometry, ipData) * gradNi;  // pressure
                    // Capillary force. stress: weak form of +sigma~eps (gradC x gradC).n (a flux
                    // term, tested with gradN_i); wellbalanced/potential: volumetric body forces
                    // (tested with N_i, like the Boussinesq term below).
                    if (stForm == SurfaceTensionForm::stress)
                    {
                        const Scalar gcgn = gradC*gradNi;
                        for (int d = 0; d < dim; ++d)
                            mom[d] -= sigmaEps*gcgn*gradC[d];
                    }
                    else if (pressureRobustCap)
                        ; // wellbalanced/potential volumetric force assembled separately (H(div) recon.)
                    else if (stForm == SurfaceTensionForm::wellBalanced)
                    {
                        // grad(mu), optionally projected onto the interface normal (balanced force)
                        GlobalPosition gmu = gradMu;
                        if (normalProject)
                        {
                            const Scalar gmn = (gradMu*gradC)/(gradC*gradC + npRegEps);
                            gmu = gradC; gmu *= gmn;
                        }
                        for (int d = 0; d < dim; ++d)
                            mom[d] += Ni*c*gmu[d]; // body force -c*gradMu => residual +c*gradMu
                    }
                    else // potential
                        for (int d = 0; d < dim; ++d)
                            mom[d] -= Ni*mu*gradC[d]; // body force +mu*gradC => residual -mu*gradC
                    // Boussinesq body force (rho(c) - rho_ref) g. Dropped here when pressure-robust
                    // gravity is on (then assembled via the H(div) reconstruction, no double-count).
                    const Scalar buoyRho = pressureRobustGrav ? 0.0 : (density - refDensity);
                    for (int d = 0; d < dim; ++d)
                        residual[i][d] += qpData.weight()*(mom[d] - Ni*buoyRho*gravity[d]);

                    // ---- phase-field transport [dim] ----
                    // advection uses the limited reconstruction cHat/gradCHat (== c/gradC if the
                    // limiter is off); mobility keeps the true gradMu.
                    Scalar pf;
                    if (skewPhaseField)
                        pf = 0.5*(Ni*(v*gradCHat) - cHat*(v*gradNi)) + mobility*(gradMu*gradNi);
                    else
                        pf = -cHat*(v*gradNi) + mobility*(gradMu*gradNi); // conservative (divergence) form
                    // conservative -D_art*gradC diffusive flux (mirrors the mobility term; CV/vertex
                    // counterpart in fluxIntegral), vanishing as |v|->0 so the static ST balance stands
                    if (phaseFieldArtDiff > 0.0)
                        pf += phaseFieldArtDiff*hElem*v.two_norm()*(gradC*gradNi);
                    residual[i][phaseFieldEqIdx] += qpData.weight()*pf;

                    // ---- chemical-potential definition [dim+1] ----
                    // diagnostic: pin mu to the analytic equilibrium value (int Ni (mu - mu_eq))
                    // instead of solving mu = -sigma~eps*Lap(c) + energyScale*f'(c).
                    // (when consistentChemPot, the mu-row is assembled below for ALL dofs instead)
                    if (!consistentChemPot)
                    {
                        const Scalar cp = prescribeMu
                            ? Ni*(mu - problem.analyticChemicalPotential(ipData.global()))
                            : sigmaEps*(gradC*gradNi) - Ni*dwSource;
                        residual[i][chemicalPotentialEqIdx] += qpData.weight()*cp;
                    }
                }

                // consistent-Galerkin mu-row on ALL local dofs (vertex + edge): same integrand as the
                // edge-dof branch above, but for every dof, so the discrete mu-definition uses a single
                // (element Galerkin) projection -> constant mu becomes an exact discrete equilibrium and
                // the vertex-vs-edge grad(mu) ripple driving the parasitic current is removed.
                if (consistentChemPot)
                    for (const auto& localDof : localDofs(fvGeometry))
                    {
                        const auto i = localDof.index();
                        const auto gradNi = ipCache.gradN(i);
                        const auto Ni = shapeValues[i][0];
                        const Scalar cp = prescribeMu
                            ? Ni*(mu - problem.analyticChemicalPotential(ipData.global()))
                            : sigmaEps*(gradC*gradNi) - Ni*dwSource;
                        residual[i][chemicalPotentialEqIdx] += qpData.weight()*cp;
                    }
            }

            // Manual high-order Dune integration of the double-well +N_i E c(c^2-1) on the mu-equation
            // (evaluate the P2 basis directly, bypassing the order-4 element cache).
            if (manualDW)
            {
                const auto& localBasis = fvGeometry.feLocalBasis();
                const auto elemGeo = fvGeometry.elementGeometry();
                std::vector<RangeType> shapeVals;
                using CT = typename FVElementGeometry::GridGeometry::GridView::ctype;
                const auto& dwQuad = Dune::QuadratureRules<CT, dim>::rule(elemGeo.type(), dwQuadOrder);
                for (const auto& qp : dwQuad)
                {
                    localBasis.evaluateFunction(qp.position(), shapeVals);
                    Scalar cc(0.0);
                    for (const auto& localDof : localDofs(fvGeometry))
                        cc += shapeVals[localDof.index()][0]*elemVars[localDof.index()].phaseField();
                    const auto w = qp.weight()*elemGeo.integrationElement(qp.position());
                    const Scalar dw = energyScale*cc*(cc*cc - 1.0);
                    // consistentChemPot: the vertex mu-rows come from the FE Galerkin path too, so the
                    // manual double-well must feed ALL dofs (the box path adds nothing to the mu-row).
                    if (consistentChemPot)
                        for (const auto& localDof : localDofs(fvGeometry))
                            residual[localDof.index()][chemicalPotentialEqIdx] += w*shapeVals[localDof.index()][0]*dw;
                    else
                        for (const auto& localDof : nonCVLocalDofs(fvGeometry))
                            residual[localDof.index()][chemicalPotentialEqIdx] += w*shapeVals[localDof.index()][0]*dw;
                }
            }
        }
    }
};

} // end namespace Dumux

#endif
