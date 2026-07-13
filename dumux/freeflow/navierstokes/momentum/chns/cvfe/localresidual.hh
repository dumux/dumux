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
#include <dune/common/fmatrix.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

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
        // pressure-robust capillary force is assembled separately (addPressureRobustCapillary_) via
        // an H(div) test-function reconstruction; when on, skip the raw volumetric force here.
        static const bool pressureRobust = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.PressureRobustCapillary", false);

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
            if (!pressureRobust && stForm == SurfaceTensionForm::wellBalanced)
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
            else if (!pressureRobust && stForm == SurfaceTensionForm::potential)
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
        addSemiImplicitDoubleWellCorrection_(residual, problem, fvGeometry, prevElemVars, curElemVars);
    }

    void addToElementFluxAndSourceResidual(ElementResidualVector& residual,
                                           const Problem& problem,
                                           const Element& element,
                                           const FVElementGeometry& fvGeometry,
                                           const ElementVariables& elemVars) const
    {
        FeResidual::addFluxAndSourceTerms(residual, problem, fvGeometry, elemVars);
        addSkewSymmetricAdvectionCorrections_(residual, problem, fvGeometry, elemVars);
        addPressureRobustCapillary_(residual, problem, element, fvGeometry, elemVars);
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
        addPressureRobustCapillary_(residual, problem, element, fvGeometry, elemVars);
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

    /*!
     * \brief Semi-implicit full Taylor linearization of the CH double-well f'(c)=c^3-c on the CV
     * (vertex) dofs (Aland & Voigt 2011, eq. after (30)): linearize about the previous time level,
     * f'(c^{n+1}) ~ (c^n)^3 - c^n + (3(c^n)^2-1)(c^{n+1}-c^n) = -2(c^n)^3 + (3(c^n)^2-1) c^{n+1},
     * which is LINEAR in c^{n+1} -> the whole CH sub-problem is linear -> unconditionally solvable
     * (removes the monolithic-Newton divergence on the stiff double-well at physical M/eps that Eyre,
     * which only makes the CONCAVE -Ec part explicit while leaving the convex c^3 implicit/nonlinear,
     * does NOT cure). Added here (the storage hook has both time levels) as the box-integrated
     * correction E*(f'_lin - f'_impl) = E*(-2(c^n)^3 + 3(c^n)^2 c^{n+1} - (c^{n+1})^3) on the mu-row,
     * exactly cancelling the sourceIntegral's implicit E*(c^{n+1})^3 (both use the vertex value +
     * box volume) and leaving the linear part. Use INSTEAD of EyreConvexSplitting and with the
     * default DoubleWellQuadOrder=0 (the box/cache double-well this correction is derived against).
     * Gated by FreeFlow.SemiImplicitDoubleWell (default off).
     */
    void addSemiImplicitDoubleWellCorrection_(ElementResidualVector& residual,
                                              const Problem& problem,
                                              const FVElementGeometry& fvGeometry,
                                              const ElementVariables& prevElemVars,
                                              const ElementVariables& curElemVars) const
    {
        static const bool semiImplicit = getParamFromGroup<bool>(
            problem.paramGroup(), "FreeFlow.SemiImplicitDoubleWell", false);
        if (!semiImplicit)
            return;
        const Scalar energyScale = problem.energyScale();
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto i = scv.localDofIndex();
            const Scalar cOld = prevElemVars[scv].phaseField();
            const Scalar cNew = curElemVars[scv].phaseField();
            const auto vol = Extrusion::volume(fvGeometry, scv) * curElemVars[scv].extrusionFactor();
            const Scalar corr = energyScale*(-2.0*cOld*cOld*cOld + 3.0*cOld*cOld*cNew - cNew*cNew*cNew);
            residual[i][chemicalPotentialEqIdx] += corr*vol;
        }
    }

    /*!
     * \brief Pressure-robust capillary body force (Linke reconstruction). Tests the capillary force
     * against the H(div)-conforming RT0 reconstruction Pi(v_h) of the velocity test functions instead
     * of the raw P2 test functions. Since Pi(v_h) is (edge-flux) divergence-controlled, the GRADIENT
     * part of the force - which the P1 pressure cannot represent, producing the parasitic current -
     * becomes L2-orthogonal to Pi(v_h): (grad psi, Pi v_h) = -(psi, div Pi v_h) ~ 0. Only the genuine
     * solenoidal capillary response survives. Element-local RT0 basis psi_e = (x - a_opp)/(2|T|) has
     * unit flux through its own (outward) edge and zero through the others; because the reconstruction
     * flux int_e v.n and psi_e both flip with the edge orientation, the product is orientation-
     * invariant -> H(div) conformity holds with NO global bookkeeping and the numeric Jacobian stays
     * element-stencil-local. Assembly: residual[k][d] += sum_e (n_e.e_d)(int_e N_k) int_T g.psi_e,
     * with the tested field g = c*gradMu (wellBalanced) / -mu*gradC (potential), int_e N_k = Le/6
     * (P2 vertex endpoint) or 2Le/3 (P2 edge midpoint). Replaces the raw volumetric capillary source
     * (gated off in sourceIntegral/felocalresidual when this is on). Gated FreeFlow.PressureRobustCapillary.
     */
    void addPressureRobustCapillary_(ElementResidualVector& residual,
                                     const Problem& problem,
                                     const Element& element,
                                     const FVElementGeometry& fvGeometry,
                                     const ElementVariables& elemVars) const
    {
        static const bool pr = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.PressureRobustCapillary", false);
        // Optionally also reconstruct the Boussinesq buoyancy (rho(c)-rho_ref) g into H(div): it is
        // gradient-dominated at the sharp density interface and, tested against the raw P2 N_i, is the
        // source of the density-contrast+gravity parasitic current that capillary-only PR does NOT
        // remove (confirmed by physics isolation: contrast-only and gravity-only both stay bounded at
        // mu~+-700; only their product blows up mu to +-2e5 -> the interfacial buoyancy force).
        static const bool prGravity = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.PressureRobustGravity", false);
        if (!pr && !prGravity)
            return;
        const auto stForm = problem.surfaceTensionForm();
        // capillary is reconstructed only as a body force (wellBalanced/potential); in the stress form
        // it is a face flux handled in felocalresidual, so then only the buoyancy (if on) is done here.
        const bool doCap = pr && (stForm != SurfaceTensionForm::stress);
        if (!doCap && !prGravity)
            return;
        const auto& gravityVec = problem.gravity();
        const Scalar refDensityPR = problem.referenceDensity();

        const auto& geo = element.geometry();
        const auto refEl = Dune::referenceElement(geo);
        const Scalar area = geo.volume();
        const int numV = refEl.size(dim);
        const int numE = refEl.size(dim-1);

        std::vector<GlobalPosition> corner(numV);
        for (int v = 0; v < numV; ++v)
            corner[v] = geo.corner(v);
        std::vector<int> oppV(numE, 0), evA(numE, 0), evB(numE, 0);
        for (int e = 0; e < numE; ++e)
        {
            evA[e] = refEl.subEntity(e, dim-1, 0, dim);
            evB[e] = refEl.subEntity(e, dim-1, 1, dim);
            for (int v = 0; v < numV; ++v)
                if (v != evA[e] && v != evB[e]) oppV[e] = v;
        }

        // reconstruction order: 0 = RT0 (constant edge flux, div=P0); >=1 = BDM1 (linear edge flux,
        // full [P1]^2, removes the P1-part of the gradient force the P1 pressure can't absorb).
        static const int prOrder = getParamFromGroup<int>(problem.paramGroup(), "FreeFlow.PressureRobustOrder", 1);

        // per-edge geometry: outward normal nE[e], length L[e]
        std::vector<GlobalPosition> nE(numE);
        std::vector<Scalar> L(numE);
        for (int e = 0; e < numE; ++e)
        {
            GlobalPosition tang = corner[evB[e]]; tang -= corner[evA[e]];
            L[e] = tang.two_norm();
            GlobalPosition n; n[0] = tang[1]; n[1] = -tang[0]; n /= L[e];
            GlobalPosition mid = corner[evA[e]]; mid += corner[evB[e]]; mid *= 0.5;
            GlobalPosition toOpp = corner[oppV[e]]; toOpp -= mid;
            if (n*toOpp > 0.0) n *= -1.0; // outward (away from opposite vertex)
            nE[e] = n;
        }

        // force moments: RT0 uses GM_e = int_T g.psi_e (psi_e=(x-a_opp)/2|T|); BDM1 uses the P1-vertex
        // moments G_v = int_T g N_v^{P1} (lambda_v). Accumulate both in one quadrature sweep.
        std::vector<Scalar> GM(numE, 0.0);
        std::vector<GlobalPosition> Gv(numV, GlobalPosition(0.0));
        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, element))
        {
            const auto& qpIp = qpData.ipData();
            const auto& ipCache = cache(elemVars, qpIp);
            const auto& sv = ipCache.shapeValues();
            Scalar c(0.0), mu(0.0), rho(0.0);
            GlobalPosition gradC(0.0), gradMu(0.0);
            for (const auto& localDof : localDofs(fvGeometry))
            {
                const auto j = localDof.index();
                c += sv[j][0]*elemVars[j].phaseField();
                mu += sv[j][0]*elemVars[j].chemicalPotential();
                rho += sv[j][0]*elemVars[j].density();
                gradC.axpy(elemVars[j].phaseField(), ipCache.gradN(j));
                gradMu.axpy(elemVars[j].chemicalPotential(), ipCache.gradN(j));
            }
            GlobalPosition gfield(0.0); // tested field g (capillary + optional buoyancy)
            if (doCap)
            {
                if (stForm == SurfaceTensionForm::wellBalanced)
                    for (int d = 0; d < dim; ++d) gfield[d] += c*gradMu[d];
                else // potential
                    for (int d = 0; d < dim; ++d) gfield[d] += -mu*gradC[d];
            }
            if (prGravity)
                for (int d = 0; d < dim; ++d) gfield[d] -= (rho - refDensityPR)*gravityVec[d];
            const auto xg = qpIp.global();
            const auto xl = qpIp.local();
            const Scalar w = qpData.weight();
            // P1 barycentric shapes on the Dune reference triangle
            const std::array<Scalar,3> lam = { Scalar(1.0)-xl[0]-xl[1], xl[0], xl[1] };
            for (int e = 0; e < numE; ++e)
            {
                GlobalPosition psi = xg; psi -= corner[oppV[e]]; psi /= (2.0*area);
                GM[e] += w*(gfield*psi);
            }
            for (int v = 0; v < numV; ++v)
                Gv[v].axpy(w*lam[v], gfield);
        }

        if (prOrder <= 0)
        {
            // RT0 distribution over the P2 velocity dofs on each edge (Le/6 vertex, 2Le/3 midpoint)
            for (int e = 0; e < numE; ++e)
                for (const auto& localDof : localDofs(fvGeometry))
                {
                    const auto k = localDof.index();
                    const auto pos = ipData(fvGeometry, localDof).global();
                    GlobalPosition mid = corner[evA[e]]; mid += corner[evB[e]]; mid *= 0.5;
                    Scalar intN = 0.0;
                    if ((pos-corner[evA[e]]).two_norm() < 1e-8*L[e] || (pos-corner[evB[e]]).two_norm() < 1e-8*L[e])
                        intN = L[e]/6.0;
                    else if ((pos-mid).two_norm() < 1e-8*L[e])
                        intN = 2.0*L[e]/3.0;
                    if (intN != 0.0)
                        for (int d = 0; d < dim; ++d)
                            residual[k][d] += nE[e][d]*intN*GM[e];
                }
            return;
        }

        // RT1 (order>=2): reconstruct into RT1 = [P1]^2 (+) {xtilde*xi, xtilde*eta}, div=P1, matching
        // the 6 edge P1-flux moments + 2 interior [P0]^2 moments -> annihilates the P1-gradient part of
        // the capillary force the P1 pressure cannot represent (the O(h) residual that RT0/BDM1, both
        // div=P0, leave). 8-dof monomial basis (centroid-shifted for conditioning), per-element 8x8
        // moment solve: residual[k][d] += (F^T A^{-1}) . b(N_k e_d), F_beta=int_T g.b_beta.
        if (prOrder >= 2)
        {
            const auto xc = geo.center();
            // normalise the monomials by the element size so the 8x8 moment matrix stays well-
            // conditioned: with raw xi~h the columns span O(1)..O(h^2) -> cond(A)~h^-3, and A.invert()
            // then returns garbage (observed: RT1 removed LESS buoyancy current than RT0, impossible if
            // correct). Scaling is a pure change of basis (same RT1 space) so the reconstruction is
            // mathematically invariant; it only fixes the floating-point conditioning.
            const Scalar hscale = std::sqrt(area);
            auto bvec = [&](const GlobalPosition& x, std::array<GlobalPosition,8>& b){
                const Scalar xi = (x[0]-xc[0])/hscale, et = (x[1]-xc[1])/hscale;
                b[0]={1.0,0.0}; b[1]={xi,0.0}; b[2]={et,0.0};
                b[3]={0.0,1.0}; b[4]={0.0,xi}; b[5]={0.0,et};
                b[6]={xi*xi, xi*et}; b[7]={xi*et, et*et};
            };
            Dune::FieldMatrix<Scalar,8,8> A(0.0);
            std::array<Scalar,8> F{};
            std::array<GlobalPosition,8> bb;
            for (const auto& qpData : CVFE::quadratureRule(fvGeometry, element)) // interior moments + F
            {
                const auto& qpIp = qpData.ipData();
                const auto& ipCache = cache(elemVars, qpIp);
                const auto& sv = ipCache.shapeValues();
                Scalar c(0.0), mu(0.0), rho(0.0); GlobalPosition gradC(0.0), gradMu(0.0);
                for (const auto& ld : localDofs(fvGeometry)) {
                    const auto j = ld.index();
                    c += sv[j][0]*elemVars[j].phaseField(); mu += sv[j][0]*elemVars[j].chemicalPotential();
                    rho += sv[j][0]*elemVars[j].density();
                    gradC.axpy(elemVars[j].phaseField(), ipCache.gradN(j));
                    gradMu.axpy(elemVars[j].chemicalPotential(), ipCache.gradN(j));
                }
                GlobalPosition gf(0.0);
                if (doCap) {
                    if (stForm==SurfaceTensionForm::wellBalanced) for (int d=0;d<dim;++d) gf[d]+=c*gradMu[d];
                    else for (int d=0;d<dim;++d) gf[d]+=-mu*gradC[d];
                }
                if (prGravity) for (int d=0;d<dim;++d) gf[d] -= (rho - refDensityPR)*gravityVec[d];
                bvec(qpIp.global(), bb);
                const Scalar w = qpData.weight();
                for (int be=0; be<8; ++be) { F[be]+=w*(gf*bb[be]); A[6][be]+=w*bb[be][0]; A[7][be]+=w*bb[be][1]; }
            }
            const auto& gauss = Dune::QuadratureRules<Scalar,1>::rule(Dune::GeometryTypes::line, 4);
            for (int e=0;e<numE;++e) { // edge P1-flux moments (rows 2e, 2e+1)
                GlobalPosition tang = corner[evB[e]]; tang -= corner[evA[e]];
                for (const auto& gp : gauss) {
                    const Scalar t = gp.position()[0];
                    GlobalPosition p = corner[evA[e]]; p.axpy(t, tang);
                    bvec(p, bb);
                    const Scalar wL = gp.weight()*L[e];
                    for (int be=0;be<8;++be) {
                        const Scalar bn = bb[be]*nE[e];
                        A[2*e][be]   += wL*bn*(1.0-t);
                        A[2*e+1][be] += wL*bn*t;
                    }
                }
            }
            A.invert();
            std::array<Scalar,8> FA{}; // F^T A^{-1}
            for (int a=0;a<8;++a){ Scalar s=0.0; for(int b2=0;b2<8;++b2) s+=F[b2]*A[b2][a]; FA[a]=s; }
            for (const auto& ld : localDofs(fvGeometry)) {
                const auto k = ld.index();
                const auto posk = ipData(fvGeometry, ld).global();
                bool isMid=false;
                for (int e=0;e<numE;++e){ GlobalPosition mid=corner[evA[e]]; mid+=corner[evB[e]]; mid*=0.5; if((posk-mid).two_norm()<1e-8*L[e]) isMid=true; }
                for (int d=0; d<dim; ++d) {
                    std::array<Scalar,8> bm{};
                    for (int e=0;e<numE;++e) {
                        GlobalPosition mid=corner[evA[e]]; mid+=corner[evB[e]]; mid*=0.5;
                        const Scalar edn = nE[e][d];
                        Scalar iA=0.0,iB=0.0;
                        if ((posk-corner[evA[e]]).two_norm()<1e-8*L[e]) iA=L[e]/6.0;
                        else if ((posk-corner[evB[e]]).two_norm()<1e-8*L[e]) iB=L[e]/6.0;
                        else if ((posk-mid).two_norm()<1e-8*L[e]) { iA=L[e]/3.0; iB=L[e]/3.0; }
                        bm[2*e]=edn*iA; bm[2*e+1]=edn*iB;
                    }
                    const Scalar intNk = isMid ? area/3.0 : 0.0;
                    bm[6] = (d==0)?intNk:0.0; bm[7] = (d==1)?intNk:0.0;
                    Scalar contrib=0.0; for(int a=0;a<8;++a) contrib += FA[a]*bm[a];
                    residual[k][d] += contrib;
                }
            }
            return;
        }

        // BDM1: reconstruct Pi(N_k e_d) as a [P1]^2 field matching the two P1 edge-flux moments per
        // edge, then residual[k][d] += sum_v Gv[v].U_v. U_v is found per vertex from its two
        // edge-normal projections (U_v.n_e), which come from the per-edge 2x2 P1 edge-mass inverse.
        // precompute per-vertex inverse of [n_ea; n_eb] (its two incident edges)
        std::vector<std::array<int,2>> vEdges(numV); // the (up to 2) edges incident to vertex v
        for (int v = 0; v < numV; ++v)
        {
            int cnt = 0;
            for (int e = 0; e < numE && cnt < 2; ++e)
                if (evA[e] == v || evB[e] == v) vEdges[v][cnt++] = e;
        }
        for (const auto& localDof : localDofs(fvGeometry))
        {
            const auto k = localDof.index();
            const auto posk = ipData(fvGeometry, localDof).global();
            for (int d = 0; d < dim; ++d)
            {
                // U_v.n_e per (vertex,edge); nonzero only for edges carrying dof k
                std::array<std::array<Scalar,3>,3> Un{}; // Un[v][e]
                for (int e = 0; e < numE; ++e)
                {
                    GlobalPosition mid = corner[evA[e]]; mid += corner[evB[e]]; mid *= 0.5;
                    // edge moments b_i,b_j of the test function N_k e_d on edge e
                    Scalar bA = 0.0, bB = 0.0;
                    const Scalar edn = nE[e][d]; // e_d . n_e
                    if ((posk-corner[evA[e]]).two_norm() < 1e-8*L[e])       bA = edn*L[e]/6.0; // k = endpoint A
                    else if ((posk-corner[evB[e]]).two_norm() < 1e-8*L[e])  bB = edn*L[e]/6.0; // k = endpoint B
                    else if ((posk-mid).two_norm() < 1e-8*L[e]) { bA = edn*L[e]/3.0; bB = edn*L[e]/3.0; } // k = midpoint
                    if (bA == 0.0 && bB == 0.0) continue;
                    // 2x2 P1 edge-mass inverse: (U_A.n, U_B.n) = (1/L)[[4,-2],[-2,4]](bA,bB)
                    Un[evA[e]][e] = (4.0*bA - 2.0*bB)/L[e];
                    Un[evB[e]][e] = (-2.0*bA + 4.0*bB)/L[e];
                }
                // per vertex: solve [n_e1; n_e2] U_v = (U_v.n_e1, U_v.n_e2)
                for (int v = 0; v < numV; ++v)
                {
                    const int e1 = vEdges[v][0], e2 = vEdges[v][1];
                    const Scalar a11=nE[e1][0], a12=nE[e1][1], a21=nE[e2][0], a22=nE[e2][1];
                    const Scalar det = a11*a22 - a12*a21;
                    if (std::abs(det) < 1e-30) continue;
                    const Scalar r1 = Un[v][e1], r2 = Un[v][e2];
                    GlobalPosition Uv;
                    Uv[0] = ( a22*r1 - a12*r2)/det;
                    Uv[1] = (-a21*r1 + a11*r2)/det;
                    residual[k][d] += Gv[v]*Uv;
                }
            }
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
