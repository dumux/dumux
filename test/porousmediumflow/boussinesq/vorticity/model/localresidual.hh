// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Local residual for the vector-potential Boussinesq model (dimensional, 2D/3D).
 *
 * Primary variables: A_k (vector-potential components) and C_i (concentrations).
 *
 * The Darcy velocity is recovered as u = ∇ × A, which satisfies ∇·u = 0 identically.
 * Each component A_k is governed by the divergence-form equation
 *
 *   ∇·F_Ak = 0,   F_Ak = (μ/K) ∇A_k − ρ₀(Σᵢ βᵢ Cᵢ)(ê_k × g)
 *
 * derived by taking the curl of Darcy's law under the Boussinesq approximation
 * (∇·u = 0).  This is the correct form for spatially variable K(x) and
 * concentration-dependent μ(C): K appears as a resistance weight on the
 * diffusive gradient, not as a prefactor on the buoyancy term.
 * Incompressibility (∇·u = 0) is the only assumption; full compressibility
 * would reintroduce a baroclinic pressure coupling that cannot be eliminated.
 *
 * In 2D only the z-component (ψ) is non-trivial, reducing to:
 *   F_ψ = (μ/K) ∇ψ − ρ₀ β C (ê_z × g)   (ê_z × g = (−g_y, g_x))
 *
 * Each scalar transport C_i satisfies:
 *   φ ∂Cᵢ/∂t + ∇·(u Cᵢ − Dᵢ(C) ∇Cᵢ) = 0
 *
 * Permeability K(x) and gravity g come from the spatial parameters.
 * Fluid properties (μ(C), ρ₀, βᵢ, Dᵢ(C)) come from the FluidSystem;
 * μ and D may depend on the local concentration vector.
 */
#ifndef DUMUX_BOUSSINESQ_VORTICITY_LOCAL_RESIDUAL_HH
#define DUMUX_BOUSSINESQ_VORTICITY_LOCAL_RESIDUAL_HH

#include <array>

#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/assembly/cvfelocalresidual.hh>

namespace Dumux {

template<class TypeTag>
class BoussinesqVorticityLocalResidual
    : public CVFELocalResidual<TypeTag>
{
    using ParentType = CVFELocalResidual<TypeTag>;

    using Scalar         = GetPropType<TypeTag, Properties::Scalar>;
    using Problem        = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry   = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView       = typename GridGeometry::GridView;
    using Element        = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry    = typename GridGeometry::LocalView;
    using SubControlVolume     = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ElementVolumeVariables =
        typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables  = GetPropType<TypeTag, Properties::VolumeVariables>;
    using ElementFluxVariablesCache =
        typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using ModelTraits  = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices      = typename ModelTraits::Indices;
    using FluidSystem  = GetPropType<TypeTag, Properties::FluidSystem>;
    using NumEqVector  = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using Extrusion    = Extrusion_t<GridGeometry>;
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimensionworld>;

    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr int nPot  = ModelTraits::numPotentialEqs;
    static constexpr int nComp = ModelTraits::numFluidComponents();

public:
    using ParentType::ParentType;

    /*!
     * \brief Storage term.
     *
     * The vector-potential equations are quasi-static (no time derivative).
     * Each transport component contributes φ · Cᵢ.
     */
    NumEqVector computeStorage(const Problem&,
                               const SubControlVolume&,
                               const VolumeVariables& volVars) const
    {
        NumEqVector storage(0.0);
        for (int i = 0; i < nComp; ++i)
            storage[Indices::transportEqIdx(i)] = volVars.porosity() * volVars.concentration(i);
        return storage;
    }

    /*!
     * \brief Flux term.
     *
     * Vector-potential equation k:
     *   F_Ak · n = ((μ/K) ∇A_k − ρ₀ S (ê_k × g)) · n · |σ|
     *   with  S = Σᵢ βᵢ C̄ᵢ  (face-averaged buoyancy concentration)
     *         μ = μ(C̄)       (viscosity evaluated at face concentrations)
     *         K = K(x)        (permeability at face centre)
     *
     * Transport equation i:
     *   (u C̃ᵢ − Dᵢ(C̄) ∇Cᵢ) · n · |σ|   (upwind advection + Fickian diffusion)
     *   with  u = ∇ × A  (curl of the vector potential)
     */
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        const auto& fluxCache = elemFluxVarsCache[scvf];
        const auto& n    = scvf.unitOuterNormal();
        const Scalar area = Extrusion::area(fvGeometry, scvf);

        // Accumulate shape-function weighted gradients and face values
        std::array<GlobalPosition, nPot>  gradA;
        std::array<GlobalPosition, nComp> gradC;
        std::array<Scalar, nComp>         Cface;
        for (auto& v : gradA) v  = GlobalPosition(0.0);
        for (auto& v : gradC) v  = GlobalPosition(0.0);
        Cface.fill(Scalar(0.0));

        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& vv    = elemVolVars[scv];
            const auto& gradN = fluxCache.gradN(scv.indexInElement());
            const Scalar N    = fluxCache.shapeValues()[scv.indexInElement()][0];

            for (int k = 0; k < nPot; ++k)
                gradA[k].axpy(vv.vectorPotential(k), gradN);
            for (int i = 0; i < nComp; ++i)
            {
                gradC[i].axpy(vv.concentration(i), gradN);
                Cface[i] += N * vv.concentration(i);
            }
        }

        // Medium properties from spatial params (may vary in space)
        const auto& sp         = problem.spatialParams();
        const Scalar K         = sp.permeabilityAtPos(scvf.center());
        const auto&  g         = sp.gravity(scvf.center());
        // Fluid properties from FluidSystem (μ and D may depend on concentration)
        const Scalar mu        = FluidSystem::viscosity(Cface);
        const Scalar rho       = FluidSystem::referenceDensity();
        // Flow resistance μ/K weights the diffusive part of F_Ak.
        // K and μ do NOT appear in the buoyancy term — see derivation.
        const Scalar resistance = mu / K;
        // Buoyancy concentration  S = Σᵢ βᵢ C̄ᵢ,  buoyancy flux = ρ₀ S (ê_k × g)
        Scalar S = 0.0;
        for (int i = 0; i < nComp; ++i)
            S += FluidSystem::solutalExpansionCoefficient(i) * Cface[i];
        const Scalar B = rho * S;

        NumEqVector flux(0.0);

        // ----- Vector-potential equations ----------------------------------
        // F_Ak = (μ/K) ∇A_k − ρ₀ S (ê_k × g)
        if constexpr (dimWorld == 2)
        {
            // Single potential (k ≡ z): ê_z × g_2D = (-g_y, g_x)
            GlobalPosition eZcrossG{ -g[1], g[0] };
            auto F = gradA[0];
            F *= resistance;       // weight gradient by flow resistance μ/K
            F.axpy(-B, eZcrossG);  // subtract buoyancy ρ₀ S (ê_z × g)
            flux[Indices::vectorPotentialEqIdx(0)] = (F * n) * area;
        }
        else
        {
            // Three potentials: (ê_k × g)_j = -g_l, (ê_k × g)_l = g_j
            // with j = (k+1)%3,  l = (k+2)%3
            for (int k = 0; k < nPot; ++k)
            {
                const int j = (k + 1) % 3;
                const int l = (k + 2) % 3;
                GlobalPosition ekCrossG(0.0);
                ekCrossG[j] = -g[l];
                ekCrossG[l] =  g[j];
                auto F = gradA[k];
                F *= resistance;       // weight gradient by flow resistance μ/K
                F.axpy(-B, ekCrossG);  // subtract buoyancy ρ₀ S (ê_k × g)
                flux[Indices::vectorPotentialEqIdx(k)] = (F * n) * area;
            }
        }

        // ----- Darcy velocity  u = ∇ × A ----------------------------------
        Scalar volumetricFlux;
        if constexpr (dimWorld == 2)
        {
            // ψ ≡ A[0]:  u = (∂ψ/∂y, -∂ψ/∂x)
            volumetricFlux = (gradA[0][1]*n[0] - gradA[0][0]*n[1]) * area;
        }
        else
        {
            // A[0]=A_x, A[1]=A_y, A[2]=A_z
            const Scalar ux = gradA[2][1] - gradA[1][2];
            const Scalar uy = gradA[0][2] - gradA[2][0];
            const Scalar uz = gradA[1][0] - gradA[0][1];
            volumetricFlux  = (ux*n[0] + uy*n[1] + uz*n[2]) * area;
        }

        // ----- Transport equations  (upwind advection + Fickian diffusion) -
        for (int i = 0; i < nComp; ++i)
        {
            const Scalar insideC  = elemVolVars[scvf.insideScvIdx()].concentration(i);
            const Scalar outsideC = elemVolVars[scvf.outsideScvIdx()].concentration(i);
            const Scalar Cup      = (volumetricFlux >= 0) ? insideC : outsideC;
            const Scalar Di       = FluidSystem::solutalDiffusionCoefficient(i, Cface);
            flux[Indices::transportEqIdx(i)] = Cup * volumetricFlux
                                             - Di * (gradC[i] * n) * area;
        }

        return flux;
    }
};

} // namespace Dumux

#endif
