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
#include <dumux/discretization/fvgridvariables.hh>
#include <dumux/discretization/cvfe/gridvariablescache_.hh>
#include <dumux/discretization/cvfe/hybrid/gridvariablescache.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>
#include <dumux/discretization/cvfe/variablesadapter.hh>
#include <dumux/discretization/fem/gridvariablescache.hh>

#ifndef DOXYGEN
namespace Dumux::Detail {

template<class TypeTag, bool = Properties::hasDefinedType<TypeTag, Properties::Variables>()>
struct DefaultVariables
{ using type = GetPropType<TypeTag, Properties::Variables>; };

template<class TypeTag>
struct DefaultVariables<TypeTag, false>
{ using type = Dumux::Detail::CVFE::VariablesAdapter<GetPropType<TypeTag, Properties::VolumeVariables>>; };

} // end namespace Dumux::Detail
#endif // DOXYGEN

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
 * \tparam Variables The variables per local degree of freedom, by default the property `Variables`
 *         or, for models that do not define it, the volume variables
 */
template<class TypeTag,
         class Variables = typename Dumux::Detail::DefaultVariables<TypeTag>::type>
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

namespace Dumux::Detail {

/*!
 * \ingroup Discretization
 * \brief Whether a type tag assembles with the variables per local degree of freedom
 *
 * This holds for all element discretizations with degrees of freedom without a sub-control
 * volume, and for the others if `Properties::Transition::UseVariables` is set.
 */
template<class TypeTag>
constexpr bool usesVariables()
{
    using ElementDiscretization = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    if constexpr (Experimental::Concepts::HybridElementDiscretization<ElementDiscretization>
                  || Experimental::Concepts::FEElementDiscretization<ElementDiscretization>)
        return true;
    else if constexpr (Properties::hasDefinedType<TypeTag, Properties::Transition::UseVariables>())
        return getPropValue<TypeTag, Properties::Transition::UseVariables>();
    else
        return false;
}

#ifndef DOXYGEN
template<class TypeTag, bool = usesVariables<TypeTag>()>
struct CVFEGridVariables
{
    using type = FVGridVariables<GetPropType<TypeTag, Properties::GridGeometry>,
                                 GetPropType<TypeTag, Properties::GridVolumeVariables>,
                                 GetPropType<TypeTag, Properties::GridFluxVariablesCache>>;
};

template<class TypeTag>
struct CVFEGridVariables<TypeTag, true>
{ using type = Experimental::DefaultGridVariables<TypeTag>; };
#endif // DOXYGEN

/*!
 * \ingroup Discretization
 * \brief The grid variables of a control-volume finite element type tag:
 *        `Experimental::DefaultGridVariables` if `usesVariables<TypeTag>()`,
 *        otherwise `FVGridVariables` with the volume variables
 */
template<class TypeTag>
using CVFEGridVariables_t = typename CVFEGridVariables<TypeTag>::type;

} // end namespace Dumux::Detail

#endif
