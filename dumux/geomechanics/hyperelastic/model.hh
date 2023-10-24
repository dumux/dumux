// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Hyperelastic
 * \brief Hyperelastic model
 * \brief Deformation of a solid body using the theory of (nonlinear) elasticity (large deformations)
 *
 * This model describes the deformation of a solid body using the finite strain theory.
 * The displacement vector of a material point following the deformation function
 * \f$ \mathbf{x} = \chi(\mathbf{X}) \f$ is given by the difference of current position
 * and position in the reference configuration:
 \f[
 \mathbf{d} = \mathbf{x} - \mathbf{X}.
 \f]
 * An important kinematic quantity characterizing the deformation is the deformation gradient
 \f[
 \mathbf{F} = \frac{\partial \mathbf{x}}{\partial \mathbf{X}} = \mathbf{I} + \nabla \mathbf{d}.
 \f]
 * The equilibrium equation (nonlinear elastostatics) in the reference frame is given by
 \f[
 - \nabla\cdot\mathbf{P} = \mathbf{f}
 \f]
 * where \f$ \mathbf{f} \f$ in \f$ \mathrm{N/m^3} \f$ is the external force acting
 * on the body per unit volume and \f$ \mathbf{P} = \mathbf{P}(\mathbf{F}) \f$ is
 * the first Piola-Kirchhoff stress tensor which relates the stress in the current
 * and the reference configuration. It is related to the Cauchy stress,
 * \f$ \mathbf{P} = J \boldsymbol{\sigma} \mathbf{F}^{-T} \f$,
 * where \f$ J = \operatorname{det}{\mathbf{F}} \f$.
 *
 * A suitable constitutive law for \f$ \mathbf{P} = \mathbf{P}(\mathbf{F}) \f$ completes
 * the model. All hyperelastic materials can be described in terms of a strain energy
 * density function \f$ \psi(\mathbf{F}) \f$ @cite ogden1997
 * and the first Piola-Kirchhoff stress tensor can be computed as
 \f[
 \mathbf{P} = \frac{\partial \psi}{\partial \mathbf{F}}
            = \mathbf{F}2\frac{\partial \psi}{\partial \mathbf{C}},
 \f]
 * where \f$ \mathbf{C} = \mathbf{F}^T\mathbf{F} \f$ is the right Cauchy-Green tensor
 * and the term \f$ \mathbf{S} = 2 (\partial \psi / \partial \mathbf{C}) \f$ is the
 * second Piola-Kirchhoff stress tensor.
 * This model expects the user problem implementation to provide a function
 * `firstPiolaKirchhoffStressTensor(F)` implementing the constitutive law.
 */
#ifndef DUMUX_GEOMECHANICS_HYPERELASTIC_MODEL_HH
#define DUMUX_GEOMECHANICS_HYPERELASTIC_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>

#include "spatialparams.hh"
#include "localresidual.hh"
#include "volumevariables.hh"

namespace Dumux {

template<class PV, class MT>
struct HyperelasticVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using ModelTraits = MT;
};

struct HyperelasticIndices
{
    static constexpr int displacementIdx(int i) { return i; };
    static constexpr int equationIdx(int i) { return i; };
};

/*!
 * \ingroup Hyperelastic
 * \brief HyperelasticModelTraits
 */
template<int dim>
struct HyperelasticModelTraits
{
    using Indices = HyperelasticIndices;
    static constexpr int numEq() { return dim; }
};

namespace Properties {

// Create new type tags
namespace TTag {
struct Hyperelastic { using InheritsFrom = std::tuple<ModelProperties>; };
} // end namespace TTag

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::Hyperelastic>
{ using type = HyperelasticModelTraits<GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimension>; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::Hyperelastic>
{ using type = HyperelasticLocalResidual<TypeTag>; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::Hyperelastic>
{
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using Traits = HyperelasticVolumeVariablesTraits<PV, MT>;
    using type = HyperelasticVolumeVariables<Traits>;
};

// Set the default spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Hyperelastic>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = DefaultHyperelasticSpatialParams<GridGeometry, Scalar>;
};

} // end namespace Properties
} // end namespace Dumux

#endif
