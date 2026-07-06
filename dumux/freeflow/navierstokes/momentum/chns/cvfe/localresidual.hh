// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \brief Local residual for the co-located momentum + Cahn-Hilliard CVFE model (P2-P1-P2-P2).
 *
 * numEq = dim + 2: [ velocity (dim), phase field c, chemical potential mu ]. The momentum block
 * is assembled by the volvars-local momentum flux helper (NavierStokesMomentumCHNSFluxCVFE) as a
 * dim-vector and embedded into components [0..dim-1]; the Cahn-Hilliard c-transport and
 * mu-definition are ported from mass/2p/localresidual.hh with the advecting velocity taken LOCALLY
 * (co-located, no coupling manager). Non-CV (P2 edge) dofs are handled by the FE weak-form helper
 * (NavierStokesMomentumCHNSFELocalResidualTerms) via the addToElement* hooks. Pressure is the only
 * inter-domain coupling (from the P1 pressure subdomain).
 */
#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_CHNS_CVFE_LOCAL_RESIDUAL_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_CHNS_CVFE_LOCAL_RESIDUAL_HH

#include <algorithm>
#include <cmath>
#include <vector>
#include <dune/common/fvector.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/concepts/variables_.hh>
#include <dumux/common/typetraits/localdofs_.hh>

#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

#include <dumux/freeflow/navierstokes/momentum/chns/cvfe/flux.hh>
#include <dumux/freeflow/navierstokes/momentum/chns/cvfe/felocalresidual.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Element-wise residual for the co-located momentum + Cahn-Hilliard CVFE model.
 */
template<class TypeTag>
class NavierStokesMomentumCHNSCVFELocalResidual
: public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridVariablesCache = Concept::GridVariablesCache_t<GridVariables>;
    using ElementVariables = typename GridVariablesCache::LocalView;
    using Variables = Concept::Variables_t<GridVariables>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    using Extrusion = Extrusion_t<GridGeometry>;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr int phaseFieldEqIdx = Indices::phaseFieldEqIdx;
    static constexpr int chemicalPotentialEqIdx = Indices::chemicalPotentialEqIdx;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using MomVec = Dune::FieldVector<Scalar, dim>;

    using LocalBasis = typename GridGeometry::FeCache::FiniteElementType::Traits::LocalBasisType;
    using FluxHelper = NavierStokesMomentumCHNSFluxCVFE<GridGeometry, MomVec>;
    using FluxFunctionHelper = NavierStokesMomentumCHNSFluxFunctionCVFE<GridGeometry, MomVec>;
    using FeResidual = NavierStokesMomentumCHNSFELocalResidualTerms<Scalar, NumEqVector, MomVec, Indices, LocalBasis, Extrusion>;
    using SurfaceTensionForm = typename Problem::SurfaceTensionForm;

public:
    using ElementResidualVector = typename ParentType::ElementResidualVector;
    using ParentType::ParentType;

    //! storage: momentum rho u (conservative), phase field c; mu is algebraic (0)
    NumEqVector computeStorage(const Problem& problem,
                               const FVElementGeometry& fvGeometry,
                               const SubControlVolume& scv,
                               const Variables& vars,
                               const bool isPreviousStorage) const
    {
        NumEqVector storage(0.0);
        const auto mom = vars.density()*vars.velocity();
        for (int d = 0; d < dim; ++d)
            storage[d] = mom[d];
        storage[phaseFieldEqIdx] = vars.phaseField();
        // storage[chemicalPotentialEqIdx] = 0 (algebraic mu definition)
        return storage;
    }

    //! flux across a sub-control-volume face: momentum (embedded) + Cahn-Hilliard scalars
    template<class ElementFluxVariablesCache>
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVariables& elemVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        // momentum flux as a dim-vector (rho/nu from volvars, pressure from coupled subdomain)
        using FluxContext = NavierStokesMomentumFluxContext<Problem, FVElementGeometry, ElementVariables, ElementFluxVariablesCache>;
        FluxContext context(problem, fvGeometry, elemVars, elemFluxVarsCache, scvf);
        FluxHelper fluxHelper;
        MomVec mflux(0.0);
        mflux += fluxHelper.advectiveMomentumFlux(context);
        mflux += fluxHelper.diffusiveMomentumFlux(context);
        mflux += fluxHelper.pressureContribution(context);

        // Cahn-Hilliard fluxes (ported from mass/2p) with LOCAL advecting velocity
        const auto& fluxVarCache = elemFluxVarsCache[scvf];
        const auto& shapeValues = fluxVarCache.shapeValues();
        MomVec v(0.0);
        for (const auto& localDof : localDofs(fvGeometry))
            v.axpy(shapeValues[localDof.index()][0], elemVars[localDof.index()].velocity());
        const Scalar volumeFlux = (v*scvf.unitOuterNormal())*scvf.area();

        GlobalPosition gradC(0.0), gradMu(0.0);
        for (const auto& localDof : localDofs(fvGeometry))
        {
            const auto& vv = elemVars[localDof.index()];
            gradC.axpy(vv.phaseField(), fluxVarCache.gradN(localDof.index()));
            gradMu.axpy(vv.chemicalPotential(), fluxVarCache.gradN(localDof.index()));
        }

        // Korteweg capillary stress face flux +sigma~eps (gradC.n) gradC on the momentum block.
        // Only for SurfaceTensionForm::stress; the wellbalanced/potential forms are volumetric
        // (see computeSource) since they have no natural divergence/flux structure.
        if (problem.surfaceTensionForm() == SurfaceTensionForm::stress)
        {
            const Scalar gcn = gradC*scvf.unitOuterNormal();
            const Scalar stressScale = problem.surfaceTension()*gcn
                * Extrusion::area(fvGeometry, scvf) * elemVars[scvf.insideScvIdx()].extrusionFactor();
            for (int d = 0; d < dim; ++d)
                mflux[d] += stressScale*gradC[d];
        }

        // mobility (optionally degenerate M0 (c^2-1)^2, clamped)
        Scalar mobility = problem.mobility();
        static const bool degenerateMobility = getParamFromGroup<bool>(
            problem.paramGroup(), "FreeFlow.DegenerateMobility", false);
        if (degenerateMobility)
        {
            const Scalar cFace = 0.5*(elemVars[scvf.insideScvIdx()].phaseField()
                                    + elemVars[scvf.outsideScvIdx()].phaseField());
            const Scalar g = std::max(0.0, 1.0 - cFace*cFace);
            mobility *= g*g;
        }

        static const Scalar phiUpwindWeight = getParamFromGroup<Scalar>(
            problem.paramGroup(), "FreeFlow.PhaseFieldUpwindWeight", 1.0);

        NumEqVector flux(0.0);
        for (int d = 0; d < dim; ++d)
            flux[d] = mflux[d];

        // phase-field transport: upwind(c) volumeFlux  - M grad(mu).n dA
        flux[phaseFieldEqIdx] = upwindSchemeMultiplier_(
            elemVars, scvf, volumeFlux,
            [](const auto& vv) { return vv.phaseField(); }, phiUpwindWeight
        )*volumeFlux;
        flux[phaseFieldEqIdx] += -1.0*vtmv(scvf.unitOuterNormal(), mobility, gradMu)*scvf.area();

        // chemical potential: -sigma~eps grad(c).n dA
        flux[chemicalPotentialEqIdx] = -1.0*vtmv(scvf.unitOuterNormal(), problem.surfaceTension(), gradC)*scvf.area();

        return flux;
    }

    //! source at a control volume (nodal values at the vertex): Boussinesq body force on momentum
    //! and the double-well mu-source; the capillary stress is a face flux (see computeFlux).
    //! NOTE: this classic (non-quadrature) interface has no gradient access, so it cannot add the
    //! wellbalanced/potential volumetric capillary force (see sourceIntegral() below, which is the
    //! one actually called by the hybrid assembler and does add it).
    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVariables& elemVars,
                              const SubControlVolume& scv) const
    {
        NumEqVector source(0.0);
        const auto& vv = elemVars[scv];
        const Scalar c = vv.phaseField();
        const auto& g = problem.gravity();
        // Boussinesq body force (rho(c) - rho_ref) g on the momentum block
        const Scalar buoyancy = vv.density() - problem.referenceDensity();
        for (int d = 0; d < dim; ++d)
            source[d] = buoyancy*g[d];
        // chemical-potential double-well source mu - energyScale f'(c), f'(c) = c(c^2-1)
        source[chemicalPotentialEqIdx] = vv.chemicalPotential() - problem.energyScale()*c*(c*c - 1.0);
        return source;
    }

    //! quadrature-based storage integral (used by the hybrid assembler): momentum rho u + phase field c
    NumEqVector storageIntegral(const FVElementGeometry& fvGeometry,
                                const ElementVariables& elemVars,
                                const SubControlVolume& scv,
                                bool isPreviousTimeLevel) const
    {
        const auto& vars = elemVars[scv];
        NumEqVector storage(0.0);
        const auto mom = vars.density()*vars.velocity();
        for (int d = 0; d < dim; ++d)
            storage[d] = mom[d];
        storage[phaseFieldEqIdx] = vars.phaseField();
        storage *= Extrusion::volume(fvGeometry, scv) * vars.extrusionFactor();
        return storage;
    }

    //! quadrature-based source integral: Boussinesq body force + double-well mu-source, plus the
    //! volumetric capillary force (wellbalanced -c*gradMu / potential mu*gradC); the stress form
    //! is a face flux instead (see fluxIntegral) and contributes nothing here.
    NumEqVector sourceIntegral(const FVElementGeometry& fvGeometry,
                               const ElementVariables& elemVars,
                               const SubControlVolume& scv) const
    {
        const auto& problem = this->asImp().problem();
        const auto& g = problem.gravity();
        const Scalar refDensity = problem.referenceDensity();
        const Scalar energyScale = problem.energyScale();
        const auto stForm = problem.surfaceTensionForm();
        const bool prescribeMu = problem.prescribeChemicalPotential(); // diagnostic mu-pin
        // Quadrature test: integrate the cubic double-well E c(c^2-1) with a manual Dune rule over the
        // scv (bypassing the Midpoint scv rule) instead of the low-order cache; matches the FE side.
        static const int dwQuadOrder = getParamFromGroup<int>(problem.paramGroup(), "FreeFlow.DoubleWellQuadOrder", 0);
        const bool manualDW = dwQuadOrder > 0 && !prescribeMu;
        // consistent-Galerkin mu-equation on all dofs -> the box path contributes no mu-row source
        static const bool consistentChemPot = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.ConsistentChemPotGalerkin", false);
        // balanced-force normal projection of the wellBalanced capillary force (see felocalresidual)
        static const bool normalProject = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.CapillaryNormalProjection", false);
        static const Scalar npRegEps = getParamFromGroup<Scalar>(problem.paramGroup(), "FreeFlow.CapillaryNormalProjectionRegEps", 1e-8);

        NumEqVector source(0.0);
        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scv))
        {
            const auto& ipCache = cache(elemVars, qpData.ipData());
            const auto& sv = ipCache.shapeValues();
            Scalar c(0.0), mu(0.0), density(0.0);
            GlobalPosition gradC(0.0), gradMu(0.0);
            for (const auto& localDof : localDofs(fvGeometry))
            {
                const auto j = localDof.index();
                c += sv[j][0]*elemVars[j].phaseField();
                mu += sv[j][0]*elemVars[j].chemicalPotential();
                density += sv[j][0]*elemVars[j].density();
                gradC.axpy(elemVars[j].phaseField(), ipCache.gradN(j));
                gradMu.axpy(elemVars[j].chemicalPotential(), ipCache.gradN(j));
            }
            // NOTE sign: evalSource() does residual -= sourceIntegral(), so this returns the
            // physical body force S (the residual then gets -N_i S). This matches the validated
            // pq1bubble baseline's problem.source() exactly: buoyancy S=(rho-rho_ref)*g, capillary
            // wellbalanced S=-c*gradMu, potential S=+mu*gradC. (felocalresidual.hh assembles the
            // residual contribution -N_i S directly, i.e. +N_i c gradMu / -N_i mu gradC, which is
            // the same physics.)
            NumEqVector s(0.0);
            for (int d = 0; d < dim; ++d)
                s[d] = (density - refDensity)*g[d];
            if (stForm == SurfaceTensionForm::wellBalanced)
            {
                GlobalPosition gmu = gradMu;
                if (normalProject) // project grad(mu) onto the interface normal n=gradC/|gradC|
                {
                    const Scalar gmn = (gradMu*gradC)/(gradC*gradC + npRegEps);
                    gmu = gradC; gmu *= gmn;
                }
                for (int d = 0; d < dim; ++d)
                    s[d] -= c*gmu[d];
            }
            else if (stForm == SurfaceTensionForm::potential)
                for (int d = 0; d < dim; ++d)
                    s[d] += mu*gradC[d];
            // chem-pot equation: normally mu - energyScale f'(c) (the -sigma~eps*Lap(c) part is
            // the face flux in fluxIntegral); diagnostic pin makes residual -= s enforce mu = mu_eq.
            // box mu-row source; zeroed when the mu-eq is assembled by consistent Galerkin on all dofs
            // (the FE path then carries the entire vertex mu-row: mass, Laplacian and double-well).
            s[chemicalPotentialEqIdx] = consistentChemPot ? Scalar(0.0)
                : (prescribeMu
                    ? problem.analyticChemicalPotential(qpData.ipData().global()) - mu
                    : (manualDW ? mu : mu - energyScale*c*(c*c - 1.0)));
            source += qpData.weight()*s;
        }

        // manual high-order Dune integration of the double-well over the scv (see FE side): the
        // chem-pot source carries -E c(c^2-1), so residual -= source adds +int E c(c^2-1).
        // Skipped under consistentChemPot: the FE Galerkin path integrates the double-well for the
        // vertex mu-rows too (over the element), so the box must not add it again.
        if (manualDW && !consistentChemPot)
        {
            const auto& localBasis = fvGeometry.feLocalBasis();
            const auto scvGeo = fvGeometry.geometry(scv);
            const auto elemGeo = fvGeometry.elementGeometry();
            using CT = typename GridView::ctype;
            std::vector<typename LocalBasis::Traits::RangeType> shapeVals;
            Scalar dwInt(0.0);
            for (const auto& qp : Dune::QuadratureRules<CT, dim>::rule(scvGeo.type(), dwQuadOrder))
            {
                const auto elemLocal = elemGeo.local(scvGeo.global(qp.position()));
                localBasis.evaluateFunction(elemLocal, shapeVals);
                Scalar cc(0.0);
                for (const auto& localDof : localDofs(fvGeometry))
                    cc += shapeVals[localDof.index()][0]*elemVars[localDof.index()].phaseField();
                dwInt += qp.weight()*scvGeo.integrationElement(qp.position())*energyScale*cc*(cc*cc - 1.0);
            }
            source[chemicalPotentialEqIdx] -= dwInt;
        }

        source *= elemVars[scv].extrusionFactor();
        return source;
    }

    //! quadrature-based flux integral over a face: momentum (embedded) + Cahn-Hilliard face fluxes
    NumEqVector fluxIntegral(const FVElementGeometry& fvGeometry,
                             const ElementVariables& elemVars,
                             const SubControlVolumeFace& scvf) const
    {
        const auto& problem = this->asImp().problem();
        const Scalar sigmaEps = problem.surfaceTension();
        const bool addStressForm = problem.surfaceTensionForm() == SurfaceTensionForm::stress;
        const bool prescribeMu = problem.prescribeChemicalPotential(); // diagnostic mu-pin
        static const bool degenerateMobility = getParamFromGroup<bool>(
            problem.paramGroup(), "FreeFlow.DegenerateMobility", false);
        static const Scalar phaseFieldArtDiff = getParamFromGroup<Scalar>(
            problem.paramGroup(), "FreeFlow.PhaseFieldArtificialDiffusion", 0.0);
        static const bool phaseFieldLimiter = getParamFromGroup<bool>(
            problem.paramGroup(), "FreeFlow.PhaseFieldMonotonicityLimiter", false);
        static const Scalar phaseFieldLimiterBound = getParamFromGroup<Scalar>(
            problem.paramGroup(), "FreeFlow.PhaseFieldLimiterBound", 0.25);
        // Abels-Garcke-Grün diffusive momentum flux (variable-density consistency); mirrors the
        // pq1bubble baseline's interfaceFlux. Only meaningful with inertia. Default off.
        static const bool enableInterfaceFlux = getParamFromGroup<bool>(
            problem.paramGroup(), "FreeFlow.EnableInterfaceFlux", false);
        const bool enableAGG = enableInterfaceFlux && problem.enableInertiaTerms();
        const Scalar dRhodC = problem.mixtureDensityDerivative();
        // consistent-Galerkin mu-equation on all dofs -> the box path adds no mu-row (see felocalresidual)
        static const bool consistentChemPot = getParamFromGroup<bool>(
            problem.paramGroup(), "FreeFlow.ConsistentChemPotGalerkin", false);

        FluxFunctionHelper fluxFunctionHelper;
        using FluxFunctionContext = NavierStokesMomentumFluxFunctionContext<Problem, FVElementGeometry, ElementVariables, typename GridVariablesCache::InterpolationPointData>;

        const auto& element = fvGeometry.element();
        const Scalar hElem = std::pow(element.geometry().volume(), 1.0/Scalar(dim));
        // element-local monotonicity-limited phase-field coeffs for the advective reconstruction
        // (same cLim as the FE edge terms -> consistent, conservative); cHat below == c if off.
        std::vector<Scalar> cLim;
        if (phaseFieldLimiter)
            FeResidual::limitedPhaseFieldCoeffs(fvGeometry, elemVars, cLim, phaseFieldLimiterBound);
        MomVec mflux(0.0);
        MomVec velIntegral(0.0);
        Scalar pfFlux(0.0), cpFlux(0.0);
        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scvf))
        {
            const auto& ipCache = cache(elemVars, qpData.ipData());
            FluxFunctionContext context(problem, fvGeometry, elemVars, ipCache);
            velIntegral += context.velocity() * qpData.weight();

            const auto& sv = ipCache.shapeValues();
            const auto n = qpData.ipData().unitOuterNormal();
            const auto& u = context.velocity();
            const auto& gradV = context.gradVelocity();

            // interpolate c, viscosity, grad c, grad mu at the ip from the local volvars
            Scalar c(0.0), viscosity(0.0); GlobalPosition gradC(0.0), gradMu(0.0);
            Scalar cHat(0.0); // limited reconstruction of the advected phase field (== c if off)
            for (const auto& localDof : localDofs(fvGeometry))
            {
                const auto j = localDof.index();
                c += sv[j][0]*elemVars[j].phaseField();
                viscosity += sv[j][0]*elemVars[j].effectiveViscosity();
                gradC.axpy(elemVars[j].phaseField(), ipCache.gradN(j));
                gradMu.axpy(elemVars[j].chemicalPotential(), ipCache.gradN(j));
                cHat += sv[j][0]*(phaseFieldLimiter ? cLim[j] : elemVars[j].phaseField());
            }

            // momentum: viscous + pressure integrands (nu from volvars, pressure from coupled p domain)
            MomVec diffusive = mv(gradV + getTransposed(gradV), n);
            diffusive *= -viscosity;
            const Scalar pressure = problem.pressure(element, fvGeometry, qpData.ipData()) - problem.referencePressure();
            MomVec pn(n); pn *= pressure;
            mflux += qpData.weight()*(diffusive + pn);
            // Korteweg capillary stress on momentum (SurfaceTensionForm::stress only; the
            // wellbalanced/potential volumetric forms are added in sourceIntegral instead)
            if (addStressForm)
            {
                const Scalar gcn = gradC*n;
                for (int d = 0; d < dim; ++d)
                    mflux[d] += qpData.weight()*sigmaEps*gcn*gradC[d];
            }

            // Cahn-Hilliard face fluxes (central)
            Scalar mobility = problem.mobility();
            if (degenerateMobility)
            {
                const Scalar gg = std::max(0.0, 1.0 - c*c);
                mobility *= gg*gg;
            }
            const auto un = u*n;
            // AGG diffusive momentum flux: momentum advected with the CH mass flux J=(drho/dc)(-M grad mu),
            // i.e. add u (J.n) = -M (drho/dc)(grad mu . n) u to the momentum flux (baseline interfaceFlux).
            if (enableAGG)
            {
                const Scalar JdotN = -mobility*dRhodC*(gradMu*n);
                for (int d = 0; d < dim; ++d)
                    mflux[d] += qpData.weight()*u[d]*JdotN;
            }
            pfFlux += qpData.weight()*(cHat*un - mobility*(gradMu*n)); // advect the limited cHat
            // conservative artificial-diffusion flux -D_art*gradC.n (D_art=c_stab*h*|u|); FE (edge)
            // counterpart in felocalresidual, so total is a divergence -> no source in the c-balance
            if (phaseFieldArtDiff > 0.0)
                pfFlux += qpData.weight()*(-phaseFieldArtDiff*hElem*u.two_norm()*(gradC*n));
            // box mu-row Laplacian face flux; dropped when the mu-eq is assembled by consistent Galerkin
            // on all dofs (vertex mu-row then comes entirely from the FE path) or under the mu-pin.
            if (!prescribeMu && !consistentChemPot)
                cpFlux += qpData.weight()*(-sigmaEps*(gradC*n));
        }
        mflux += fluxFunctionHelper.advectiveMomentumFluxIntegral(problem, fvGeometry, elemVars, scvf, velIntegral);

        const auto extrusion = elemVars[fvGeometry.scv(scvf.insideScvIdx())].extrusionFactor();
        NumEqVector flux(0.0);
        for (int d = 0; d < dim; ++d)
            flux[d] = mflux[d]*extrusion;
        flux[phaseFieldEqIdx] = pfFlux*extrusion;
        flux[chemicalPotentialEqIdx] = cpFlux*extrusion;
        return flux;
    }

    void addToElementStorageResidual(ElementResidualVector& residual,
                                     const Problem& problem,
                                     const Element& element,
                                     const FVElementGeometry& fvGeometry,
                                     const ElementVariables& prevElemVars,
                                     const ElementVariables& curElemVars) const
    {
        FeResidual::addStorageTerms(residual, problem, fvGeometry, prevElemVars, curElemVars,
                                    this->timeLoop().timeStepSize());
        addConvectiveMomentumStorageCorrection_(residual, problem, fvGeometry, prevElemVars, curElemVars);
        addEyreConvexSplittingCorrection_(residual, problem, fvGeometry, prevElemVars, curElemVars);
    }

    void addToElementFluxAndSourceResidual(ElementResidualVector& residual,
                                           const Problem& problem,
                                           const Element& element,
                                           const FVElementGeometry& fvGeometry,
                                           const ElementVariables& elemVars) const
    {
        FeResidual::addFluxAndSourceTerms(residual, problem, fvGeometry, elemVars);
        addSkewSymmetricAdvectionCorrections_(residual, problem, fvGeometry, elemVars);
    }

    void addToElementFluxAndSourceResidual(ElementResidualVector& residual,
                                           const Problem& problem,
                                           const Element& element,
                                           const FVElementGeometry& fvGeometry,
                                           const ElementVariables& elemVars,
                                           const ElementBoundaryTypes& elemBcTypes) const
    {
        FeResidual::addFluxAndSourceTerms(residual, problem, fvGeometry, elemVars);
        addSkewSymmetricAdvectionCorrections_(residual, problem, fvGeometry, elemVars);
    }

private:
    template<class ElemVolVars, class UpwindTermFunction>
    Scalar upwindSchemeMultiplier_(const ElemVolVars& elemVolVars,
                                   const SubControlVolumeFace& scvf,
                                   const Scalar flux,
                                   const UpwindTermFunction& upwindTerm,
                                   const Scalar upwindWeight) const
    {
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        using std::signbit;
        if (signbit(flux))
            return upwindWeight*upwindTerm(outsideVolVars) + (1.0 - upwindWeight)*upwindTerm(insideVolVars);
        else
            return upwindWeight*upwindTerm(insideVolVars) + (1.0 - upwindWeight)*upwindTerm(outsideVolVars);
    }

    /*!
     * \brief Skew-symmetric (Temam) corrections of the CV (vertex-dof) momentum and phase-field
     * advection, ported from momentum/cvfe/localresidual.hh and mass/2p/localresidual.hh. Both
     * fluxIntegral()'s momentum advection and the CH phase-field advection above are assembled in
     * CONSERVATIVE (divergence) form; this is energy-/monotonicity-stable only if the advecting
     * velocity is discretely mass-conservative on the momentum control volumes, which it is not
     * here (div v = 0 is enforced on a DIFFERENT stencil than the Cahn-Hilliard-transported
     * density/phase field lives on). Converting to convective form removes the resulting parasitic
     * interfacial energy injection (momentum) / overshoot (phase field). Both corrections vanish
     * when the discrete divergence is exactly zero, so they are consistent, not just stabilizing.
     * Gated by FreeFlow.SkewSymmetricMomentumAdvection / FreeFlow.SkewSymmetricPhaseFieldAdvection
     * (both default off, matching the base momentum/mass models).
     *
     * Computed here (rather than via the elemFluxVarsCache-based face loop the base models use)
     * because the hybrid CVFE assembler drives this model through the quadrature-based *Integral
     * methods above, which recompute shape values/gradients per scvf via CVFE::quadratureRule +
     * cache(), not through a separately-cached elemFluxVarsCache.
     */
    void addSkewSymmetricAdvectionCorrections_(ElementResidualVector& residual,
                                               const Problem& problem,
                                               const FVElementGeometry& fvGeometry,
                                               const ElementVariables& elemVars) const
    {
        static const bool enableMomentumSkew = getParamFromGroup<bool>(
            problem.paramGroup(), "FreeFlow.SkewSymmetricMomentumAdvection", false);
        static const bool enablePhaseFieldSkew = getParamFromGroup<bool>(
            problem.paramGroup(), "FreeFlow.SkewSymmetricPhaseFieldAdvection", false);
        const bool enableMomentum = enableMomentumSkew && problem.enableInertiaTerms();
        if (!enableMomentum && !enablePhaseFieldSkew)
            return;

        static const auto upwindWeight = getParamFromGroup<Scalar>(problem.paramGroup(), "Flux.UpwindWeight");

        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (scvf.boundary())
                continue;

            // face-integrated velocity int_face v dA (CVFE::quadratureRule weights already
            // include the physical face measure, matching fluxIntegral's velIntegral above)
            MomVec velIntegral(0.0);
            for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scvf))
            {
                const auto& ipCache = cache(elemVars, qpData.ipData());
                const auto& sv = ipCache.shapeValues();
                MomVec v(0.0);
                for (const auto& localDof : localDofs(fvGeometry))
                    v.axpy(sv[localDof.index()][0], elemVars[localDof.index()].velocity());
                velIntegral += v*qpData.weight();
            }
            const auto n = scvf.unitOuterNormal();
            const auto volumeFlux = velIntegral*n;

            const auto insideIdx = scvf.insideScvIdx();
            const auto outsideIdx = scvf.outsideScvIdx();
            const auto& insideVars = elemVars[insideIdx];
            const auto& outsideVars = elemVars[outsideIdx];
            const auto extrusion = insideVars.extrusionFactor();

            if (enableMomentum)
            {
                // mass flux m_f = rho_f (v.n) dA, same upwind-blended density as the advective
                // flux (NavierStokesMomentumCHNSFluxFunctionCVFE::advectiveMomentumFluxIntegral)
                const auto density = volumeFlux > 0
                    ? upwindWeight*insideVars.density() + (1.0-upwindWeight)*outsideVars.density()
                    : upwindWeight*outsideVars.density() + (1.0-upwindWeight)*insideVars.density();
                const auto massFlux = density*volumeFlux*extrusion;

                const auto& uIn = insideVars.velocity();
                const auto& uOut = outsideVars.velocity();
                for (int d = 0; d < dim; ++d)
                {
                    residual[insideIdx][d] -= massFlux*uIn[d];
                    residual[outsideIdx][d] += massFlux*uOut[d];
                }
            }

            if (enablePhaseFieldSkew)
            {
                const auto pfFlux = 0.5*volumeFlux*extrusion;
                residual[insideIdx][phaseFieldEqIdx]  -= pfFlux*insideVars.phaseField();
                residual[outsideIdx][phaseFieldEqIdx] += pfFlux*outsideVars.phaseField();
            }
        }
    }

    /*!
     * \brief Convective-form momentum storage correction for CV (vertex) dofs. storageIntegral()
     * above always uses the conservative storage rho^{n+1} u^{n+1} - rho^n u^n (the framework
     * evaluates it once per time level and subtracts, so freezing density at the current level
     * within storageIntegral itself is not possible without access to the other time level). This
     * adds the correction conservative -> convective directly, since addToElementStorageResidual
     * uniquely has both prevElemVars and curElemVars at once:
     *     convective - conservative = rho^{n+1}(u^{n+1}-u^n)/dt - (rho^{n+1}u^{n+1}-rho^n u^n)/dt
     *                               = -(rho^{n+1}-rho^n) u^n / dt
     * matching the FE (P2 edge dof) treatment in NavierStokesMomentumCHNSFELocalResidualTerms.
     * Gated by FreeFlow.ConvectiveFormMomentumStorage (default off).
     */
    /*!
     * \brief Eyre convex splitting of the CH double-well on the CV (vertex) dofs: makes the concave
     * -E c part of f'(c)=E(c^3-c) explicit (-E c^n) for unconditional gradient stability. Added here
     * (the storage hook has both time levels) as +E(c^{n+1}-c^n) box-integrated onto the mu-equation,
     * turning the sourceIntegral's implicit -E c^{n+1} into -E c^n. Gated by FreeFlow.EyreConvexSplitting.
     */
    void addEyreConvexSplittingCorrection_(ElementResidualVector& residual,
                                           const Problem& problem,
                                           const FVElementGeometry& fvGeometry,
                                           const ElementVariables& prevElemVars,
                                           const ElementVariables& curElemVars) const
    {
        static const bool eyre = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.EyreConvexSplitting", false);
        if (!eyre)
            return;
        const Scalar energyScale = problem.energyScale();
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto i = scv.localDofIndex();
            const auto cDiff = curElemVars[scv].phaseField() - prevElemVars[scv].phaseField();
            const auto vol = Extrusion::volume(fvGeometry, scv) * curElemVars[scv].extrusionFactor();
            residual[i][chemicalPotentialEqIdx] += energyScale*cDiff*vol;
        }
    }

    void addConvectiveMomentumStorageCorrection_(ElementResidualVector& residual,
                                                 const Problem& problem,
                                                 const FVElementGeometry& fvGeometry,
                                                 const ElementVariables& prevElemVars,
                                                 const ElementVariables& curElemVars) const
    {
        static const bool convectiveFormStorage = getParamFromGroup<bool>(
            problem.paramGroup(), "FreeFlow.ConvectiveFormMomentumStorage", false);
        if (!convectiveFormStorage)
            return;

        const auto dt = this->timeLoop().timeStepSize();
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& curVars = curElemVars[scv];
            const auto& prevVars = prevElemVars[scv];
            const auto densityDiff = curVars.density() - prevVars.density();
            const auto vol = Extrusion::volume(fvGeometry, scv) * curVars.extrusionFactor();
            const auto& prevVelocity = prevVars.velocity();
            const auto i = scv.localDofIndex();
            for (int d = 0; d < dim; ++d)
                residual[i][d] -= densityDiff*prevVelocity[d]/dt*vol;
        }
    }
};

} // end namespace Dumux

#endif
