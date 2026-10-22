// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Grid variables with the variables cache selected depending on the element discretization
 */
#ifndef DUMUX_DISCRETIZATION_DEFAULT_GRID_VARIABLES_HH
#define DUMUX_DISCRETIZATION_DEFAULT_GRID_VARIABLES_HH

#include <type_traits>

#include <dumux/common/properties.hh>
#include <dumux/discretization/concepts.hh>
#include <dumux/discretization/gridvariables.hh>
#include <dumux/discretization/cvfe/gridvariablescache_.hh>
#include <dumux/discretization/cvfe/hybrid/gridvariablescache.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>
#include <dumux/discretization/cvfe/variablesadapter.hh>
#include <dumux/discretization/fem/gridvariablescache.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup Discretization
 * \brief The grid variables cache selected depending on the element discretization
 *
 * Element discretizations with sub-control volumes only get the control-volume cache,
 * those with additional degrees of freedom without control volume the hybrid cache,
 * and those without sub-control volumes the finite element cache.
 *
 * \tparam GridDiscretization The grid discretization
 * \tparam Problem The problem
 * \tparam Variables The variables per local degree of freedom
 * \tparam enableCache Whether the variables and the interpolation point data are stored for the whole grid
 * \tparam InterpolationPointData The data stored per interpolation point
 */
template<class GridDiscretization, class Problem, class Variables, bool enableCache,
         class InterpolationPointData = Dumux::CVFE::LocalBasisInterpolationPointData<GridDiscretization>>
using DefaultGridVariablesCache = std::conditional_t<
    Concepts::FEElementDiscretization<typename GridDiscretization::LocalView>,
    FEGridVariablesCache<FEDefaultGridVariablesCacheTraits<Problem, Variables, InterpolationPointData>, enableCache>,
    std::conditional_t<
        Concepts::HybridElementDiscretization<typename GridDiscretization::LocalView>,
        CVFE::HybridCVFEGridVariablesCache<CVFE::HybridCVFEDefaultGridVariablesCacheTraits<Problem, Variables, InterpolationPointData>, enableCache>,
        CVFE::CVFEGridVariablesCache<CVFE::CVFEDefaultGridVariablesCacheTraits<Problem, Variables, InterpolationPointData>, enableCache>
    >
>;

/*!
 * \ingroup Discretization
 * \brief The grid variables of a model whose variables are defined per local degree of freedom
 *
 * The grid discretization, the problem and the caching policy (`EnableGridVolumeVariablesCache`)
 * are taken from the type tag.
 *
 * \tparam TypeTag The type tag of the model
 * \tparam Variables The variables per local degree of freedom
 */
template<class TypeTag,
         class Variables = Dumux::Detail::CVFE::VariablesAdapter<GetPropType<TypeTag, Properties::VolumeVariables>>>
using DefaultGridVariables = GridVariables<
    GetPropType<TypeTag, Properties::GridGeometry>,
    DefaultGridVariablesCache<
        GetPropType<TypeTag, Properties::GridGeometry>,
        GetPropType<TypeTag, Properties::Problem>,
        Variables,
        getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>()
    >
>;

} // end namespace Dumux::Experimental

#endif
