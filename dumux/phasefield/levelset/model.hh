// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PhaseFieldLevelSetModel
 * \brief The phase-field level set equation
 *
 * The conservative phase-field level set represents an interface by a smoothed indicator
 * function \f$ \phi \in [0,1] \f$. Its profile across the interface is kept
 * at a fixed hyperbolic-tangent shape by evolving it in pseudo-time
 * \f$ \tau \f$ according to the compression/reinitialization equation
 \f[
 \frac{\partial \phi}{\partial \tau} + \nabla\cdot\left( \phi (1-\phi)\, \mathbf{n} \right)
   = \nabla\cdot\left( \varepsilon\, (\nabla \phi \cdot \mathbf{n})\, \mathbf{n} \right)
 \f]
 * where \f$ \mathbf{n} = \nabla \phi / |\nabla \phi| \f$ is the (regularized)
 * interface normal and \f$ \varepsilon \f$ controls the interface thickness.
 * Since \f$ (\nabla \phi \cdot \mathbf{n})\,\mathbf{n} = \nabla \phi \f$
 * (the gradient has no component perpendicular to its own direction), the
 * right-hand side term is just standard diffusion \f$ \varepsilon \nabla^2 \phi \f$;
 * the compression term \f$ \phi(1-\phi)\mathbf{n} \f$ is the only genuinely
 * new (nonlinear, direction-dependent) term compared to Allen-Cahn/Cahn-Hilliard.
 * This model implements the reinitialization step only (no advection by an
 * external velocity field), so it is a single second-order equation and can
 * be discretized directly with (piecewise linear or quadratic) CVFE schemes.
 *
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * TODO!!! The normal vector used by this model should come in as a parameter!
 * It should be computed from the phase field of the previous time step, or from a
 * separate phase field (e.g. the Hele-Shaw two-phase model's order parameter).
 *
 * Ideally what we would test here: New model phase field advection (simple advection)
 * This will be used together: Advection step, compute normal, reinitialization step,
 * Then move on and so forth. This is the conservative phase-field level set model, but it is not a complete model yet.
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 */
#ifndef DUMUX_PHASEFIELD_LEVEL_SET_MODEL_HH
#define DUMUX_PHASEFIELD_LEVEL_SET_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>

#include "indices.hh"
#include "localresidual.hh"
#include "variables.hh"

namespace Dumux::Properties {

//! Type tag for the  level-set model
namespace TTag {
struct PhaseFieldLevelSet { using InheritsFrom = std::tuple<ModelProperties>; };
} // end namespace TTag

//! Use the local residual of the  level-set model
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::PhaseFieldLevelSet>
{ using type = LevelSetLocalResidual<TypeTag>; };

//! The model traits of the  level-set model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::PhaseFieldLevelSet>
{
    struct type
    {
        using Indices = LevelSetIndices;
        static constexpr int numEq() { return 1; }
    };
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::PhaseFieldLevelSet>
{
    struct Traits
    {
        using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
        using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    };

    using type = LevelSetVariables<Traits>;
};

} // end namespace Dumux::Properties

#endif
