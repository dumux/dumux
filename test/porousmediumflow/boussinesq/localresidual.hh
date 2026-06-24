// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Local residual for the streamfunction-Boussinesq formulation.
 *
 * Primary variables: (ψ, C)   where ψ is the streamfunction and C the solute
 * mass fraction (both dimensionless).
 *
 * Equations assembled:
 *
 *   Eq 0 (replaceCompEqIdx = 0):
 *       ∇·(∇ψ − C·êₓ) = 0
 *   which is the weak form of the streamfunction Poisson equation
 *       ∇²ψ = ∂C/∂x
 *   derived from Darcy's law + incompressibility + Boussinesq buoyancy.
 *   There is NO time derivative (ψ is instantaneously determined by C).
 *
 *   Eq 1 (solute transport):
 *       φ ∂C/∂t + ∇·(u C − D∇C) = 0,    u = (∂ψ/∂y, −∂ψ/∂x)
 *   This is handled entirely by the parent CompositionalLocalResidual
 *   together with the StreamfunctionAdvection type.
 */
#ifndef DUMUX_BOUSSINESQ_STREAMFUNCTION_LOCAL_RESIDUAL_HH
#define DUMUX_BOUSSINESQ_STREAMFUNCTION_LOCAL_RESIDUAL_HH

#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/porousmediumflow/compositional/localresidual.hh>

namespace Dumux {

/*!
 * \brief Streamfunction-Boussinesq local residual.
 *
 * Inherits the transport equation (eq 1) from CompositionalLocalResidual.
 * Overrides eq 0 storage (zero — no time derivative for ψ) and eq 0 flux
 * (the streamfunction Poisson equation).
 */
template<class TypeTag>
class BoussinesqStreamfunctionLocalResidual
    : public CompositionalLocalResidual<TypeTag>
{
    using ParentType = CompositionalLocalResidual<TypeTag>;

    using Scalar         = GetPropType<TypeTag, Properties::Scalar>;
    using Problem        = GetPropType<TypeTag, Properties::Problem>;
    using FluidSystem    = GetPropType<TypeTag, Properties::FluidSystem>;
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
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using Extrusion   = Extrusion_t<GridGeometry>;

    static constexpr int dimWorld = GridView::dimensionworld;

public:
    using ParentType::ParentType;

    /*!
     * \brief Storage term: only C (eq 1) carries a time derivative.
     *
     * ψ is determined instantaneously from C, so storage[0] = 0.
     */
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        auto storage = ParentType::computeStorage(problem, scv, volVars);
        storage[0] = 0.0;
        return storage;
    }

    /*!
     * \brief Flux term.
     *
     * Eq 1 (transport) is delegated to the parent; eq 0 is replaced by the
     * streamfunction Poisson flux:
     *
     *   F₀ = (∇ψ − C·êₓ) · n · |σ|
     *
     * The face-averaged C is interpolated with CVFE shape values.
     */
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        // Parent handles diffusion + advection for the concentration equation (eq 1).
        auto flux = ParentType::computeFlux(problem, element, fvGeometry,
                                            elemVolVars, scvf, elemFluxVarsCache);

        // Override eq 0: streamfunction Poisson ∇·(∇ψ + C·êₓ) = 0
        const auto& fluxCache   = elemFluxVarsCache[scvf];
        const auto& shapeValues = fluxCache.shapeValues();

        Dune::FieldVector<Scalar, dimWorld> gradPsi(0.0);
        Scalar C_face(0.0);
        for (auto&& scv : scvs(fvGeometry))
        {
            const auto& vv = elemVolVars[scv];
            const Scalar N = shapeValues[scv.indexInElement()][0];
            gradPsi.axpy(vv.pressure(0), fluxCache.gradN(scv.indexInElement()));
            C_face += N * vv.massFraction(0, FluidSystem::soluteIdx);
        }

        // F = ∇ψ - C·êₓ   (weak form of ∇²ψ = ∂C/∂x)
        auto F = gradPsi;
        F[0] -= C_face;
        flux[0] = (F * scvf.unitOuterNormal()) * Extrusion::area(fvGeometry, scvf);

        return flux;
    }
};

} // namespace Dumux

#endif
