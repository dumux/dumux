// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \brief Momentum flux helpers for the co-located momentum + Cahn-Hilliard CVFE model.
 *
 * These are forks of NavierStokesMomentumFluxCVFE / NavierStokesMomentumFluxFunctionCVFE that read
 * the mixture density rho(c) and viscosity nu(c) from the (local) VOLUME VARIABLES instead of from
 * problem.density()/effectiveViscosity(). This is REQUIRED for a correct monolithic Jacobian: in the
 * co-located model the phase field c is a primary variable of the SAME subdomain, so the numeric
 * Jacobian deflects c per dof and rebuilds the volvars; only a volvars-local rho(c)/nu(c) responds to
 * that deflection, whereas a coupling-manager/problem lookup of a stored solution would not.
 *
 * The velocity FluxContext / FluxFunctionContext are unchanged and reused from the base flux.hh.
 * The template parameter MomVec must be a dim-sized vector (the momentum block), satisfying the
 * static_assert MomVec::dimension == dimWorld inherited from the base helpers' design.
 */
#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_CHNS_CVFE_FLUXVARIABLES_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_CHNS_CVFE_FLUXVARIABLES_HH

#include <dune/common/fmatrix.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/extrusion.hh>

// reuse the velocity contexts (NavierStokesMomentumFluxContext / ...FluxFunctionContext)
#include <dumux/freeflow/navierstokes/momentum/cvfe/flux.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief CVFE momentum flux helper with volvars-local mixture material (co-located CHNS).
 * \tparam MomVec a dim-sized momentum/velocity vector type
 */
template<class GridGeometry, class MomVec>
class NavierStokesMomentumCHNSFluxCVFE
{
    using GridView = typename GridGeometry::GridView;
    using Scalar = typename MomVec::value_type;
    using Extrusion = Extrusion_t<GridGeometry>;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using Tensor = Dune::FieldMatrix<Scalar, dim, dimWorld>;
    static_assert(MomVec::dimension == dimWorld, "MomVec must have dimension == dimWorld");

    //! interpolate the mixture viscosity nu(c) at a face from the local volvars
    template<class ElemVolVars, class FVElementGeometry, class ShapeValues>
    Scalar faceViscosity_(const ElemVolVars& elemVolVars,
                          const FVElementGeometry& fvGeometry,
                          const ShapeValues& shapeValues) const
    {
        Scalar mu(0.0);
        for (const auto& localDof : localDofs(fvGeometry))
            mu += shapeValues[localDof.index()][0] * elemVolVars[localDof.index()].effectiveViscosity();
        return mu;
    }

public:
    //! advective momentum flux sum_f m_f u_up with m_f = rho_f (v.n) dA, rho from volvars
    template<class Context>
    MomVec advectiveMomentumFlux(const Context& context) const
    {
        if (!context.problem().enableInertiaTerms())
            return MomVec(0.0);

        const auto& fvGeometry = context.fvGeometry();
        const auto& elemVolVars = context.elemVolVars();
        const auto& scvf = context.scvFace();
        const auto& fluxVarCache = context.elemFluxVarsCache()[scvf];
        const auto& shapeValues = fluxVarCache.shapeValues();

        MomVec v(0.0);
        for (const auto& localDof : localDofs(fvGeometry))
            v.axpy(shapeValues[localDof.index()][0], elemVolVars[localDof.index()].velocity());

        const auto vn = v*scvf.unitOuterNormal();
        const auto& insideVolVars = elemVolVars[fvGeometry.scv(scvf.insideScvIdx())];
        const auto& outsideVolVars = elemVolVars[fvGeometry.scv(scvf.outsideScvIdx())];
        // rho(c) from the local (deflected) volvars -> correct monolithic Jacobian
        const auto insideDensity = insideVolVars.density();
        const auto outsideDensity = outsideVolVars.density();
        const auto upwindMom = vn > 0 ? insideDensity*insideVolVars.velocity() : outsideDensity*outsideVolVars.velocity();
        const auto downwindMom = vn > 0 ? outsideDensity*outsideVolVars.velocity() : insideDensity*insideVolVars.velocity();
        static const auto upwindWeight = getParamFromGroup<Scalar>(context.problem().paramGroup(), "Flux.UpwindWeight");
        const auto advectiveTermIntegrand = vn * (upwindWeight * upwindMom + (1.0-upwindWeight)*downwindMom);

        return advectiveTermIntegrand * Extrusion::area(fvGeometry, scvf) * insideVolVars.extrusionFactor();
    }

    //! diffusive (viscous) momentum flux -nu(c) (gradV + gradV^T) . n dA, nu from volvars
    template<class Context>
    MomVec diffusiveMomentumFlux(const Context& context) const
    {
        const auto& fvGeometry = context.fvGeometry();
        const auto& elemVolVars = context.elemVolVars();
        const auto& scvf = context.scvFace();
        const auto& fluxVarCache = context.elemFluxVarsCache()[scvf];

        Tensor gradV(0.0);
        for (const auto& localDof : localDofs(fvGeometry))
        {
            const auto& volVars = elemVolVars[localDof];
            for (int dir = 0; dir < dim; ++dir)
                gradV[dir].axpy(volVars.velocity(dir), fluxVarCache.gradN(localDof.index()));
        }

        const auto mu = faceViscosity_(elemVolVars, fvGeometry, fluxVarCache.shapeValues());

        static const bool enableUnsymmetrizedVelocityGradient
            = getParamFromGroup<bool>(context.problem().paramGroup(), "FreeFlow.EnableUnsymmetrizedVelocityGradient", false);

        MomVec diffusiveFlux = enableUnsymmetrizedVelocityGradient ?
                mv(gradV, scvf.unitOuterNormal())
                : mv(gradV + getTransposed(gradV), scvf.unitOuterNormal());

        diffusiveFlux *= -mu;

        static const bool enableDilatationTerm = getParamFromGroup<bool>(context.problem().paramGroup(), "FreeFlow.EnableDilatationTerm", false);
        if (enableDilatationTerm)
            diffusiveFlux += 2.0/3.0 * mu * trace(gradV) * scvf.unitOuterNormal();

        diffusiveFlux *= Extrusion::area(fvGeometry, scvf) * elemVolVars[fvGeometry.scv(scvf.insideScvIdx())].extrusionFactor();
        return diffusiveFlux;
    }

    //! pressure contribution (p-p_ref) n dA; pressure comes from the coupled P1 pressure subdomain
    template<class Context>
    MomVec pressureContribution(const Context& context) const
    {
        const auto& element = context.element();
        const auto& fvGeometry = context.fvGeometry();
        const auto& elemVolVars = context.elemVolVars();
        const auto& scvf = context.scvFace();
        const auto& fluxVarCache = context.elemFluxVarsCache()[scvf];

        const auto pressure = context.problem().pressure(element, fvGeometry, fluxVarCache.ipData());
        const auto referencePressure = context.problem().referencePressure();

        MomVec pn(scvf.unitOuterNormal());
        pn *= (pressure-referencePressure)*Extrusion::area(fvGeometry, scvf)*elemVolVars[fvGeometry.scv(scvf.insideScvIdx())].extrusionFactor();
        return pn;
    }
};

/*!
 * \ingroup NavierStokesModel
 * \brief CVFE momentum flux-function helper (integration faces) with volvars-local material.
 */
template<class GridGeometry, class MomVec>
class NavierStokesMomentumCHNSFluxFunctionCVFE
{
    using GridView = typename GridGeometry::GridView;
    using Scalar = typename MomVec::value_type;
    using Extrusion = Extrusion_t<GridGeometry>;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using Tensor = Dune::FieldMatrix<Scalar, dim, dimWorld>;
    static_assert(MomVec::dimension == dimWorld, "MomVec must have dimension == dimWorld");

public:
    //! advective momentum flux for a given integrated face velocity; rho(c) upwind-blended from volvars
    template<class Problem, class FVElementGeometry, class ElementVariables,
             class SubControlVolumeFace, class VelocityVector>
    MomVec advectiveMomentumFluxIntegral(const Problem& problem,
                                         const FVElementGeometry& fvGeometry,
                                         const ElementVariables& elemVars,
                                         const SubControlVolumeFace& scvf,
                                         const VelocityVector& integratedVelocity) const
    {
        if (!problem.enableInertiaTerms())
            return MomVec(0.0);

        const auto vn_integral = integratedVelocity*scvf.unitOuterNormal();
        const auto& insideVolVars = elemVars[fvGeometry.scv(scvf.insideScvIdx())];
        const auto& outsideVolVars = elemVars[fvGeometry.scv(scvf.outsideScvIdx())];
        // blend the momentum rho(c) u with the SAME upwind weighting as the CV advective flux
        const auto upwindMom = vn_integral > 0 ? insideVolVars.density()*insideVolVars.velocity()
                                               : outsideVolVars.density()*outsideVolVars.velocity();
        const auto downwindMom = vn_integral > 0 ? outsideVolVars.density()*outsideVolVars.velocity()
                                                 : insideVolVars.density()*insideVolVars.velocity();
        static const auto upwindWeight = getParamFromGroup<Scalar>(problem.paramGroup(), "Flux.UpwindWeight");
        return vn_integral * (upwindWeight*upwindMom + (1.0-upwindWeight)*downwindMom);
    }

    //! diffusive momentum flux integrand; nu(c) interpolated from volvars at the ip
    template<class Context, class IpData>
    MomVec diffusiveMomentumFluxIntegrand(const Context& context, const IpData& ipData) const
    {
        const auto& fvGeometry = context.fvGeometry();
        const auto& elemVolVars = context.elemVolVars();

        Scalar mu(0.0);
        const auto& shapeValues = ipData.shapeValues();
        for (const auto& localDof : localDofs(fvGeometry))
            mu += shapeValues[localDof.index()][0] * elemVolVars[localDof.index()].effectiveViscosity();

        static const bool enableUnsymmetrizedVelocityGradient
            = getParamFromGroup<bool>(context.problem().paramGroup(), "FreeFlow.EnableUnsymmetrizedVelocityGradient", false);

        const auto& gradV = context.gradVelocity();
        MomVec diffusiveFluxIntegrand = enableUnsymmetrizedVelocityGradient ?
                mv(gradV, ipData.unitOuterNormal())
                : mv(gradV + getTransposed(gradV), ipData.unitOuterNormal());

        diffusiveFluxIntegrand *= -mu;

        static const bool enableDilatationTerm = getParamFromGroup<bool>(context.problem().paramGroup(), "FreeFlow.EnableDilatationTerm", false);
        if (enableDilatationTerm)
            diffusiveFluxIntegrand += 2.0/3.0 * mu * trace(gradV) * ipData.unitOuterNormal();

        return diffusiveFluxIntegrand;
    }

    //! pressure flux integrand; pressure from the coupled P1 pressure subdomain
    template<class Context, class IpData>
    MomVec pressureFluxIntegrand(const Context& context, const IpData& ipData) const
    {
        const auto& element = context.element();
        const auto& fvGeometry = context.fvGeometry();
        const auto pressure = context.problem().pressure(element, fvGeometry, ipData);
        const auto referencePressure = context.problem().referencePressure();

        MomVec pn(ipData.unitOuterNormal());
        pn *= (pressure-referencePressure);
        return pn;
    }
};

} // end namespace Dumux

#endif
