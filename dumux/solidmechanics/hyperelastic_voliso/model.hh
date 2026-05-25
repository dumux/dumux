// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup HyperelasticVolIso
 * \brief Hyperelastic model with volumetric-isochoric split
 *
 * \page HyperelasticVolIsoModel Hyperelastic model with volumetric-isochoric split
 *
 * # Physical model
 *
 * This model solves the finite-strain quasi-static momentum balance for a
 * hyperelastic solid using a **mixed displacement-pressure (u–p) formulation**.
 * The split between volumetric and isochoric (deviatoric) response is handled
 * by an independent pressure-like variable \f$ p_s \f$, which prevents
 * volumetric locking for (nearly-)incompressible materials.
 *
 * ## Kinematics
 *
 * The deformation map \f$ \chi : \Omega_0 \to \Omega \f$ maps material points
 * \f$ \mathbf{X} \f$ to spatial points \f$ \mathbf{x} = \chi(\mathbf{X}) \f$.
 * The displacement field is \f$ \mathbf{d} = \mathbf{x} - \mathbf{X} \f$.
 * The deformation gradient and its Jacobian are
 * \f[
 *   \mathbf{F} = \mathbf{I} + \nabla_X \mathbf{d}, \quad J = \det\mathbf{F}.
 * \f]
 *
 * ## Governing equations
 *
 * **Subdomain 0 — Momentum** (displacement \f$ \mathbf{d} \f$, PQ1Bubble or PQ2):
 * \f[
 *   -\nabla_X \cdot \mathbf{P} = \mathbf{0}, \quad
 *   \mathbf{P} = \mu\!\left(\mathbf{F} - \mathbf{F}^{-T}\right) + J p_s \mathbf{F}^{-T},
 * \f]
 * where \f$ \mu \f$ is the shear modulus (read from `spatialParams().shearModulus()`).
 * This stress corresponds to the isochoric part of the neo-Hookean energy
 * \f$ \psi_{\mathrm{iso}} = \tfrac{\mu}{2}(\mathrm{tr}\,\mathbf{C} - 3 - 2\ln J) \f$
 * augmented by the volumetric pressure contribution.
 *
 * **Subdomain 1 — Bulk pressure** (\f$ p_s \f$, Box/P1):
 * \f[
 *   p_s^{\mathrm{eq}}(J) - p_s = 0,
 * \f]
 * where the equilibrium pressure \f$ p_s^{\mathrm{eq}}(J) = \partial W_{\mathrm{vol}}/\partial J \f$
 * is defined by the strain energy function and must be supplied by the pressure problem via
 * `problem.volumetricPressure(J)`.  Two common choices are:
 * - \f$ W_{\mathrm{vol}} = \tfrac{\lambda}{2}(J-1)^2 \f$:
 *   \f$ p_s^{\mathrm{eq}} = \lambda(J-1) \f$
 * - Ogden \f$ \psi_1 \f$:
 *   \f$ p_s^{\mathrm{eq}} = \lambda(J^2-1)/(2J) \f$
 *
 * ## Discretization schemes
 *
 * The two type tags `HyperelasticVolIsoMomentumModel` and
 * `HyperelasticVolIsoPressureModel` must be combined with a suitable
 * discretization tag (e.g. `PQ1BubbleModel`/`BoxModel` for MINI or
 * `PQ2HybridModel`/`BoxModel` for Taylor-Hood) in the application's
 * properties file, and coupled via `HyperelasticVolIsoCouplingManager`.
 */
#ifndef DUMUX_SOLIDMECHANICS_HYPERELASTIC_VOLISO_MODEL_HH
#define DUMUX_SOLIDMECHANICS_HYPERELASTIC_VOLISO_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>

#include "localresidual.hh"
#include "volumevariables.hh"

namespace Dumux {

// ============================================================
// Momentum subdomain traits
// ============================================================

/*!
 * \ingroup HyperelasticVolIso
 * \brief Index set for the momentum subdomain.
 */
template<int dim>
struct HyperelasticVolIsoMomentumIndices
{
    static constexpr int displacement(int i) { return i; }
    static constexpr int momentum(int i) { return i; }
};

/*!
 * \ingroup HyperelasticVolIso
 * \brief Model traits for the momentum subdomain (\f$ \dim \f$ displacement DOFs per node).
 */
template<int dim>
struct HyperelasticVolIsoMomentumModelTraits
{
    using Indices = HyperelasticVolIsoMomentumIndices<dim>;
    static constexpr int numEq() { return dim; }
};

template<class PV, class MT>
struct HyperelasticVolIsoMomentumVVTraits
{ using PrimaryVariables = PV; using ModelTraits = MT; };

// ============================================================
// Pressure subdomain traits
// ============================================================

/*!
 * \ingroup HyperelasticVolIso
 * \brief Index set for the pressure subdomain.
 */
struct HyperelasticVolIsoPressureIndices
{
    static constexpr int pressureIdx   = 0;
    static constexpr int pressureEqIdx = 0;
};

/*!
 * \ingroup HyperelasticVolIso
 * \brief Model traits for the pressure subdomain (1 pressure DOF per node).
 */
struct HyperelasticVolIsoPressureModelTraits
{
    using Indices = HyperelasticVolIsoPressureIndices;
    static constexpr int numEq() { return 1; }
};

template<class PV, class MT>
struct HyperelasticVolIsoPressureVVTraits
{ using PrimaryVariables = PV; using ModelTraits = MT; };

} // end namespace Dumux

// ============================================================
// Property registrations
// ============================================================
namespace Dumux::Properties::TTag {
/*!
 * \ingroup HyperelasticVolIso
 * \brief Type tag for the momentum subdomain of the mixed u–p model.
 */
struct HyperelasticVolIsoMomentumModel { using InheritsFrom = std::tuple<ModelProperties>; };
/*!
 * \ingroup HyperelasticVolIso
 * \brief Type tag for the pressure subdomain of the mixed u–p model.
 */
struct HyperelasticVolIsoPressureModel { using InheritsFrom = std::tuple<ModelProperties>; };
} // end namespace Dumux::Properties::TTag

namespace Dumux::Properties {

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::HyperelasticVolIsoMomentumModel>
{
    using type = HyperelasticVolIsoMomentumModelTraits<
        GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimension>;
};
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::HyperelasticVolIsoMomentumModel>
{ using type = HyperelasticVolIsoMomentumLocalResidual<TypeTag>; };
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::HyperelasticVolIsoMomentumModel>
{
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using type = HyperelasticVolIsoMomentumVolumeVariables<
        HyperelasticVolIsoMomentumVVTraits<PV, MT>>;
};

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::HyperelasticVolIsoPressureModel>
{ using type = HyperelasticVolIsoPressureModelTraits; };
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::HyperelasticVolIsoPressureModel>
{ using type = HyperelasticVolIsoPressureLocalResidual<TypeTag>; };
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::HyperelasticVolIsoPressureModel>
{
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using type = HyperelasticVolIsoPressureVolumeVariables<
        HyperelasticVolIsoPressureVVTraits<PV, MT>>;
};

} // end namespace Dumux::Properties
#endif
