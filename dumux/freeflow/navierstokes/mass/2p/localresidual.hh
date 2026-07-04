// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesMassTwoPLocalResidual
 */
#ifndef DUMUX_FREEFLOW_NAVIERSTOKES_MASS_2P_LOCAL_RESIDUAL_HH
#define DUMUX_FREEFLOW_NAVIERSTOKES_MASS_2P_LOCAL_RESIDUAL_HH

#include <algorithm>
#include <type_traits>

#include <dumux/common/numeqvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/common/typetraits/problem.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Element-wise calculation of the Navier-Stokes residual for two-phase flow.
 */
template<class TypeTag>
class NavierStokesMassTwoPLocalResidual
: public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    using GridFluxVariablesCache = typename GridVariables::GridFluxVariablesCache;
    using ElementFluxVariablesCache = typename GridFluxVariablesCache::LocalView;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    // NB: ElementResidualVector is inherited (public) from the base local residual - do not
    // re-alias it here, that would shadow it as a private member and break the assemblers.

    using Extrusion = Extrusion_t<GridGeometry>;

    static constexpr int dimWorld = GridGeometry::GridView::dimensionworld;

public:
    //! Use the parent type's constructor
    using ParentType::ParentType;

    /*!
     * \brief Calculate the storage term of the equation
     */
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        NumEqVector storage(0.0);
        // Incompressible (Boussinesq) continuity ∇·v = 0 has no density storage and is
        // well-conditioned; the density contrast then only enters momentum (buoyancy).
        // The compressible/quasi-incompressible form ∂ρ/∂t + ∇·(ρv) = 0 keeps the density
        // but introduces a ρ'(φ)/dt stiffness in the φ-coupling that is hard on the solver.
        static const bool incompressible = getParamFromGroup<bool>(
            problem.paramGroup(), "FreeFlow.IncompressibleContinuity", true);
        storage[Indices::conti0EqIdx] = incompressible ? 0.0 : volVars.density();
        storage[Indices::phaseFieldEqIdx] = volVars.phaseField();
        storage[Indices::chemicalPotentialEqIdx] = 0.0;
        return storage;
    }

    /*!
     * \brief Evaluate the mass flux over a face of a sub control volume.
     *
     * \param problem The problem
     * \param element The element
     * \param fvGeometry The finite volume geometry context
     * \param elemVolVars The volume variables for all flux stencil elements
     * \param scvf The sub control volume face to compute the flux on
     * \param elemFluxVarsCache The cache related to flux computation
     */
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        NumEqVector flux(0.0);

        // TODO variable extrusion factor?
        const auto velocity = problem.faceVelocity(element, fvGeometry, scvf);
        const Scalar volumeFlux = velocity*scvf.unitOuterNormal()*scvf.area();
        static const bool incompressible = getParamFromGroup<bool>(
            problem.paramGroup(), "FreeFlow.IncompressibleContinuity", true);
        flux[Indices::conti0EqIdx] = incompressible
            ? volumeFlux
            : upwindSchemeMultiplier_(
                  elemVolVars, scvf, volumeFlux,
                  [](const auto& volVars) { return volVars.density(); }
              )*volumeFlux;

        const auto& fluxVarCache = elemFluxVarsCache[scvf];
        Dune::FieldVector<Scalar, dimWorld> gradPhaseField(0.0);
        Dune::FieldVector<Scalar, dimWorld> gradChemicalPotential(0.0);
        for (const auto& localDof : localDofs(fvGeometry))
        {
            const auto& volVars = elemVolVars[localDof];
            // v.axpy(a, w) means v += a*w
            gradPhaseField.axpy(
                volVars.phaseField(),
                fluxVarCache.gradN(localDof.index())
            );
            gradChemicalPotential.axpy(
                volVars.chemicalPotential(),
                fluxVarCache.gradN(localDof.index())
            );
        }

        // advection of the phase field. Default FULL UPWIND (weight 1): upwind's dissipation
        // keeps phi monotone (no overshoot), which is what suppresses the single-node interfacial
        // pressure spikes here; the diffusion does not smear the interface because the
        // Cahn-Hilliard dynamics restore its width. The energy-stable alternative is CENTRAL
        // (weight 0.5) together with FreeFlow.SkewSymmetricPhaseFieldAdvection=true (the skew form
        // 1/2[u.grad phi + div(phi u)] that is robust to discrete div u); it works but showed no
        // benefit over upwind in these benchmarks (Case-1 spike slightly worse, Case-2 identical).
        static const Scalar phiUpwindWeight = getParamFromGroup<Scalar>(
            problem.paramGroup(), "FreeFlow.PhaseFieldUpwindWeight", 1.0);
        flux[Indices::phaseFieldEqIdx] = upwindSchemeMultiplier_(
            elemVolVars, scvf, volumeFlux,
            [](const auto& volVars) { return volVars.phaseField(); },
            phiUpwindWeight
        )*volumeFlux;

        // Cahn-Hilliard diffusion flux -M grad(mu). Optionally use a DEGENERATE mobility
        // M(c) = M0 (c^2-1)^2 which vanishes in the bulk (c=+-1) and concentrates the diffusion at
        // the interface (c=0). This gives the correct interface-controlled sharp-interface limit
        // and is well-behaved even at relatively large eps (cf. Aland & Voigt 2011, who use
        // M0 = gamma ~ 1e-3*eps); a CONSTANT mobility (default) adds spurious bulk diffusion that
        // contaminates the dynamics and converges more slowly.
        Scalar mobility = problem.mobility();
        static const bool degenerateMobility = getParamFromGroup<bool>(
            problem.paramGroup(), "FreeFlow.DegenerateMobility", false);
        if (degenerateMobility)
        {
            const Scalar cFace = 0.5*(elemVolVars[scvf.insideScvIdx()].phaseField()
                                    + elemVolVars[scvf.outsideScvIdx()].phaseField());
            // clamp: for |c| > 1 (overshoot nodes) the mobility must stay 0; the raw
            // (c^2-1)^2 grows again there and would feed diffusion exactly where the
            // degenerate law is supposed to switch it off
            const Scalar g = std::max(0.0, 1.0 - cFace*cFace);
            mobility *= g*g;
        }
        flux[Indices::phaseFieldEqIdx] += -1.0*vtmv(
            scvf.unitOuterNormal(), mobility, gradChemicalPotential
        )*scvf.area();

        flux[Indices::chemicalPotentialEqIdx] = -1.0*vtmv(
            scvf.unitOuterNormal(), problem.surfaceTension(), gradPhaseField
        )*scvf.area();

        return flux;
    }

    /*!
     * \brief Compute the source term of the equation
     * \param problem The problem
     * \param element The element
     * \param fvGeometry The finite volume geometry context
     * \param elemVolVars The volume variables for all flux stencil elements
     * \param scv The sub control volume
     */
    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume& scv) const
    {
        NumEqVector source(0.0);

        // add contributions from problem (e.g. double well potential)
        source += problem.source(element, fvGeometry, elemVolVars, scv);

        return source;
    }

    /*!
     * \brief Element-level hook: optional skew-symmetric correction of the phase-field advection.
     */
    void addToElementFluxAndSourceResidual(typename ParentType::ElementResidualVector& residual,
                                           const Problem& problem,
                                           const Element& element,
                                           const FVElementGeometry& fvGeometry,
                                           const ElementVolumeVariables& elemVolVars,
                                           const ElementFluxVariablesCache& elemFluxVarsCache,
                                           const ElementBoundaryTypes& elemBcTypes) const
    {
        addSkewSymmetricPhaseFieldAdvectionCorrection_(residual, problem, element, fvGeometry, elemVolVars);
    }

private:
    /*!
     * \brief Skew-symmetric (Temam) form of the Cahn-Hilliard phase-field convection.
     *
     * computeFlux assembles phi-advection in CONSERVATIVE (divergence) form sum_f phi_up m_f
     * (m_f = (v.n) dA the volume flux through the sub-control-volume face). That equals the
     * advective form u.grad(phi) only if sum_f m_f = 0 per control volume; with incompressible
     * continuity div v = 0 enforced only discretely, sum_f m_f != 0 at the density interface, and
     * the spurious phi*(div u) term makes phi OVERSHOOT (>1) there -> mu = (gamma/eps)phi(phi^2-1)
     * - gamma*eps*Laplacian(phi) spikes -> single-node interfacial pressure spikes. The skew form
     *   1/2[u.grad phi + div(phi u)] = div(phi u) - 1/2 phi (div u)
     * cancels half that term (subtract 1/2 phi_cv m_f from each adjacent CV); it vanishes when
     * sum_f m_f = 0 (consistent) and is robust to the small discrete divergence at the interface.
     * Gated by FreeFlow.SkewSymmetricPhaseFieldAdvection (default off).
     */
    void addSkewSymmetricPhaseFieldAdvectionCorrection_(typename ParentType::ElementResidualVector& residual,
                                                        const Problem& problem,
                                                        const Element& element,
                                                        const FVElementGeometry& fvGeometry,
                                                        const ElementVolumeVariables& elemVolVars) const
    {
        static const bool enable = getParamFromGroup<bool>(
            problem.paramGroup(), "FreeFlow.SkewSymmetricPhaseFieldAdvection", false);
        if (!enable)
            return;

        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (scvf.boundary())
                continue;

            const auto velocity = problem.faceVelocity(element, fvGeometry, scvf);
            const Scalar volumeFlux = velocity*scvf.unitOuterNormal()*scvf.area();

            residual[scvf.insideScvIdx()][Indices::phaseFieldEqIdx]
                += -0.5*volumeFlux*elemVolVars[scvf.insideScvIdx()].phaseField();
            residual[scvf.outsideScvIdx()][Indices::phaseFieldEqIdx]
                += 0.5*volumeFlux*elemVolVars[scvf.outsideScvIdx()].phaseField();
        }
    }



    template<class ElemVolVars, class SubControlVolumeFace, class UpwindTermFunction, class Scalar>
    Scalar upwindSchemeMultiplier_(const ElemVolVars& elemVolVars,
                                   const SubControlVolumeFace& scvf,
                                   const Scalar flux,
                                   const UpwindTermFunction& upwindTerm,
                                   const Scalar upwindWeight = 1.0) const
    {
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        using std::signbit;
        if (signbit(flux)) // if sign of flux is negative
            return (upwindWeight*upwindTerm(outsideVolVars)
                    + (1.0 - upwindWeight)*upwindTerm(insideVolVars));
        else
            return (upwindWeight*upwindTerm(insideVolVars)
                    + (1.0 - upwindWeight)*upwindTerm(outsideVolVars));
    }

};

} // end namespace Dumux

#endif
