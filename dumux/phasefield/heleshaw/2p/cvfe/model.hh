// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup HeleShawModel
 * \brief The Hele-Shaw (Darcy) Cahn-Hilliard two-phase model, discretized
 *        with the hybrid control-volume/finite-element (CVFE) scheme.
 *
 * Three coupled equations for pressure \f$ p \f$, phase field \f$ \phi \f$,
 * chemical potential \f$ \mu \f$:
 \f[
 \nabla\cdot\left( \lambda(\phi) \left(\nabla p - \mu \nabla \phi - \rho_{mix}(\phi) \mathbf{g}\right) \right) = 0
 \f]
 \f[
 \frac{\partial \phi}{\partial t} + \nabla\cdot(\mathbf{v}\phi) = \nabla\cdot(M \nabla \mu),
   \qquad \mathbf{v} = -\lambda(\phi)\left(\nabla p - \mu \nabla \phi - \rho_{mix} \mathbf{g}\right)
 \f]
 \f[
 \mu = \frac{\partial f}{\partial \phi} - \gamma \nabla^2 \phi
 \f]
 * where \f$ \lambda = K/\eta_{mix}(\phi) \f$ is the Hele-Shaw (Darcy)
 * mobility, \f$ \mu \nabla \phi \f$ is the Korteweg capillary body force,
 * \f$ M \f$ is the Cahn-Hilliard mobility, \f$ \gamma \f$ is the surface
 * tension coefficient, and \f$ f(\phi) \f$ is a double-well free-energy
 * density on the symmetric \f$ [-1,1] \f$ range.
 *
 * This is the hybrid CVFE (PQ2-capable) port of the classic Box
 * implementation at `dumux/phasefield/heleshaw/2p/model.hh`: vertex local
 * dofs get classic control-volume storage/flux/source contributions;
 * local dofs without an associated sub-control volume (e.g. the edge dofs
 * of PQ2/PQ3) get finite-element (weak form) contributions added on top.
 * Unlike Cahn-Hilliard/Allen-Cahn/the conservative level set, this model has
 * an advective term (the Darcy-driven transport of \f$ \phi \f$); the
 * finite-element (edge-dof) treatment of this term is a plain Galerkin
 * (non-upwinded) discretization, since upwinding is inherently a
 * face-based concept that has no clean analog for dofs without an
 * associated sub-control-volume face.
 */
#ifndef DUMUX_PHASEFIELD_HELESHAW_2P_CVFE_MODEL_HH
#define DUMUX_PHASEFIELD_HELESHAW_2P_CVFE_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>

#include "indices.hh"
#include "localresidual.hh"
#include "variables.hh"

namespace Dumux::Properties {

//! Type tag for the hybrid CVFE Hele-Shaw two-phase model
namespace TTag {
struct HeleShawTwoPCVFE { using InheritsFrom = std::tuple<ModelProperties>; };
} // end namespace TTag

//! Use the local residual of the hybrid CVFE Hele-Shaw two-phase model
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::HeleShawTwoPCVFE>
{ using type = HeleShawTwoPCVFELocalResidual<TypeTag>; };

//! The model traits of the hybrid CVFE Hele-Shaw two-phase model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::HeleShawTwoPCVFE>
{
    struct type
    {
        using Indices = HeleShawTwoPCVFEIndices;
        static constexpr int numEq() { return 3; }
    };
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::HeleShawTwoPCVFE>
{
    struct Traits
    {
        using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
        using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    };

    using type = HeleShawTwoPCVFEVariables<Traits>;
};

} // end namespace Dumux::Properties

#endif
