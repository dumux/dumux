// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MembranePlate
 * \brief Membrane plate model
 *
 * A membrane is an infinitely thin plate with no bending stiffness.
 * Under a transverse load \f$ F \f$ (force per unit area), the vertical
 * deformation \f$ w \f$ satisfies the Poisson equation
 * \f[
 *   -\nabla\cdot(T\,\nabla w) = F,
 * \f]
 * where \f$ T \f$ is the membrane tension (force per unit length).
 * This is the 2D analogue of the elastic string equation.
 *
 * \par Primary variables
 * - vertical deformation \f$ w \f$
 */
#ifndef DUMUX_MEMBRANE_PLATE_MODEL_HH
#define DUMUX_MEMBRANE_PLATE_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>

#include <dumux/solidmechanics/plate/membrane/spatialparams.hh>

#include "localresidual.hh"
#include "volumevariables.hh"

namespace Dumux {

template<class PV, class MT>
struct MembranePlateVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using ModelTraits = MT;
};

struct MembranePlateIndices
{
    static constexpr int deformationIdx = 0;
    static constexpr int deformationEqIdx = 0;
};

/*!
 * \ingroup MembranePlate
 * \brief Model traits for the membrane plate model
 */
struct MembranePlateModelTraits
{
    using Indices = MembranePlateIndices;
    static constexpr int numEq() { return 1; }
};

} // end namespace Dumux

namespace Dumux::Properties::TTag {
struct MembranePlate { using InheritsFrom = std::tuple<ModelProperties>; };
} // end namespace Dumux::Properties::TTag

namespace Dumux::Properties {

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::MembranePlate>
{ using type = MembranePlateModelTraits; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::MembranePlate>
{ using type = MembranePlateLocalResidual<TypeTag>; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::MembranePlate>
{
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using Traits = MembranePlateVolumeVariablesTraits<PV, MT>;
    using type = MembranePlateVolumeVariables<Traits>;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::MembranePlate>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = MembranePlateDefaultSpatialParams<GridGeometry, Scalar>;
};

} // end namespace Dumux::Properties

#endif
