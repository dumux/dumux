// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 *
 * \brief A single-phase, isothermal Navier-Stokes model in a moving mesh
 */

#ifndef DUMUX_NAVIERSTOKES_CVFE_MOVINGMESH_MODEL_HH
#define DUMUX_NAVIERSTOKES_CVFE_MOVINGMESH_MODEL_HH

#include <dumux/common/properties.hh>

#include <dumux/freeflow/navierstokes/momentum/cvfe/model.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/movingmesh/localresidual.hh>

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
//! The type tag for the single-phase, isothermal Navier-Stokes model
struct NavierStokesMomentumMovingMeshCVFE
{ using InheritsFrom = std::tuple<NavierStokesMomentumCVFE, FreeFlow>; };
} // end namespace TTag

//! The local residual
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::NavierStokesMomentumMovingMeshCVFE>
{ using type = NavierStokesMomentumMovingMeshCVFELocalResidual<TypeTag>; };

namespace TTag {
//! The type tag for the single-phase, isothermal Navier-Stokes model
struct NavierStokesMassMovingMeshCVFE
{ using InheritsFrom = std::tuple<NavierStokesMassOneP, FreeFlow>; };
} // end namespace TTag

//! The local residual
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::NavierStokesMassMovingMeshCVFE>
{ using type = NavierStokesMassMovingMeshCVFELocalResidual<TypeTag>; };

} // end namespace Dumux::Properties

#endif
