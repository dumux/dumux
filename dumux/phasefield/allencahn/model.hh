// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup AllenCahnModel
 * \brief The Allen-Cahn equation
 *
 * The Allen-Cahn equation models the motion of a non-conserved phase-field
 * order parameter and is given by
 \f[
 \frac{\partial \phi}{\partial t} = M \left( \gamma \nabla^2 \phi - \frac{\partial f}{\partial \phi} \right)
 \f]
 * where \f$ \phi \f$ is the phase field, \f$ M \f$ is the mobility,
 * \f$ \gamma \f$ is the gradient energy coefficient, and \f$ f(\phi) \f$ is a
 * (typically double-well or double-obstacle) free energy density. Unlike the
 * Cahn-Hilliard equation, this is a single second-order equation, so it can
 * be discretized directly with (piecewise linear or quadratic) CVFE schemes.
 */
#ifndef DUMUX_PHASEFIELD_ALLENCAHN_MODEL_HH
#define DUMUX_PHASEFIELD_ALLENCAHN_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>

#include "indices.hh"
#include "localresidual.hh"
#include "variables.hh"

namespace Dumux::Properties {

//! Type tag for the Allen-Cahn model
namespace TTag {
struct AllenCahn { using InheritsFrom = std::tuple<ModelProperties>; };
} // end namespace TTag

//! Use the local residual of the Allen-Cahn model
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::AllenCahn>
{ using type = AllenCahnLocalResidual<TypeTag>; };

//! The model traits of the Allen-Cahn model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::AllenCahn>
{
    struct type
    {
        using Indices = AllenCahnIndices;
        static constexpr int numEq() { return 1; }
    };
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::AllenCahn>
{
    struct Traits
    {
        using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
        using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    };

    using type = AllenCahnVariables<Traits>;
};

} // end namespace Dumux::Properties

#endif
