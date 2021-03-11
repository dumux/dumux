// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Typetraits
 * \brief Type traits for problem classes
 */
#ifndef DUMUX_TYPETRAITS_PROBLEM_HH
#define DUMUX_TYPETRAITS_PROBLEM_HH

#include <type_traits>

#include <dumux/common/typetraits/isvalid.hh>
#include <dumux/timestepping/timelevel.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

// forward declare
namespace Detail {
template<class Problem, DiscretizationMethod dm>
struct ProblemTraits;

template<class GridGeometry, class Scalar>
struct hasTransientDirichlet
{
private:
    using TimeLevel = Dumux::Experimental::TimeLevel<Scalar>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using BoundaryEntity = std::conditional_t<GridGeometry::discMethod == DiscretizationMethod::box,
                                              typename GridGeometry::SubControlVolume,
                                              typename GridGeometry::SubControlVolumeFace>;
public:
    template<class Problem>
    auto operator()(const Problem& p)
    -> decltype(p.dirichlet(std::declval<Element>(),
                            std::declval<BoundaryEntity>(),
                            std::declval<TimeLevel>()))
    {}
};

template<class Problem, class GridGeometry, class Scalar>
inline constexpr bool hasTransientDirichletInterface
    = decltype(isValid(hasTransientDirichlet<GridGeometry, Scalar>())(std::declval<Problem>()))::value;

} // end namespace Detail

/*!
 * \ingroup Common
 * \brief Type traits for problem classes.
 */
template<class Problem>
struct ProblemTraits
{
    using GridGeometry = std::decay_t<decltype(std::declval<Problem>().gridGeometry())>;
    using BoundaryTypes = typename Detail::template ProblemTraits<Problem, GridGeometry::discMethod>::BoundaryTypes;

    static constexpr bool hasTransientDirichletInterface
        = Detail::hasTransientDirichletInterface<Problem, GridGeometry, typename Problem::Traits::Scalar>;
};

} // end namespace Dumux

#endif
