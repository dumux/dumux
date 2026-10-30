// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Typetraits
 * \brief Type traits for classes providing a grid discretization
 */
#ifndef DUMUX_TYPETRAITS_GRID_DISCRETIZATION_HH
#define DUMUX_TYPETRAITS_GRID_DISCRETIZATION_HH

#include <utility>

#include <dune/common/std/type_traits.hh>

#include <dumux/common/typetraits/typetraits.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace Detail {

template<class T>
struct GridDiscretizationType
{
    static_assert(AlwaysFalse<T>::value, "Type exports neither GridGeometry nor GridDiscretization");
};

//! new interfaces, which export GridDiscretization
template<class T>
    requires requires { typename T::GridDiscretization; }
struct GridDiscretizationType<T>
{ using type = typename T::GridDiscretization; };

//! classic types, which only export GridGeometry
template<class T>
    requires (!requires { typename T::GridDiscretization; }
              && requires { typename T::GridGeometry; })
struct GridDiscretizationType<T>
{ using type = typename T::GridGeometry; };

//! Placeholders for element discretizations that do not partition the element into
//! sub-control volumes, such as the finite element schemes
struct NoSubControlVolume {};
struct NoSubControlVolumeFace {};

template<class ED> using SpecifiedSubControlVolume = typename ED::SubControlVolume;
template<class ED> using SpecifiedSubControlVolumeFace = typename ED::SubControlVolumeFace;

} // end namespace Detail
#endif // DOXYGEN

/*!
 * \ingroup Typetraits
 * \brief The grid discretization type of a class exporting it
 * \note Generic code that has to serve both the GridGeometry and the
 *       GridDiscretization interfaces should use this instead of naming
 *       either alias directly.
 */
template<class T>
using GridDiscretization_t = typename Detail::GridDiscretizationType<T>::type;

/*!
 * \ingroup Typetraits
 * \brief The sub-control volume type of an element discretization, or a placeholder
 *        type for discretizations that do not define sub-control volumes
 * \note This allows generic code to declare sub-control-volume based interfaces without
 *       excluding finite element discretizations, for which those interfaces are unusable.
 */
template<class ElementDiscretization>
using SubControlVolume_t = Dune::Std::detected_or_t<
    Detail::NoSubControlVolume, Detail::SpecifiedSubControlVolume, ElementDiscretization
>;

/*!
 * \ingroup Typetraits
 * \brief The sub-control volume face type of an element discretization, or a placeholder
 *        type for discretizations that do not define sub-control volume faces
 * \copydetails SubControlVolume_t
 */
template<class ElementDiscretization>
using SubControlVolumeFace_t = Dune::Std::detected_or_t<
    Detail::NoSubControlVolumeFace, Detail::SpecifiedSubControlVolumeFace, ElementDiscretization
>;

/*!
 * \ingroup Typetraits
 * \brief The grid discretization
 * \note Generic code that has to serve both the gridGeometry() and the
 *       gridDiscretization() interfaces should use this instead of calling
 *       either accessor directly.
 * \note the args make it work for multidomain problems
 */
template<class T, class... Args>
decltype(auto) gridDiscretization(const T& t, Args&&... args)
{
    if constexpr (requires { t.gridDiscretization(std::forward<Args>(args)...); })
        return t.gridDiscretization(std::forward<Args>(args)...);
    else
        return t.gridGeometry(std::forward<Args>(args)...);
}

} // end namespace Dumux

#endif
