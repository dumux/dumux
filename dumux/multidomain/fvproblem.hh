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
 * \ingroup MultiDomain
 * \brief Multidomain wrapper for multiple problems
 */
#ifndef DUMUX_MULTIDOMAIN_PROBLEM_HH
#define DUMUX_MULTIDOMAIN_PROBLEM_HH

#include <tuple>
#include <memory>
#include <utility>

#include <dune/common/hybridutilities.hh>
#include <dune/common/indices.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \brief A multidomain wrapper for multiple problems
 * \tparam MDTraits The multidomain traits
 */
template<class MDTraits>
class MultiDomainFVProblem
{
    using SolutionVector = typename MDTraits::SolutionVector;
    static constexpr std::size_t numSubDomains = MDTraits::numSubDomains;

    template<std::size_t i>
    using GridGeometry = typename MDTraits::template SubDomain<i>::GridGeometry;
    using GridGeometries = typename MDTraits::template Tuple<GridGeometry>;

public:
    //! export base types of the stored type
    template<std::size_t i>
    using Type = typename MDTraits::template SubDomain<i>::Problem;

    //! export pointer types the stored type
    template<std::size_t i>
    using PtrType = std::shared_ptr<Type<i>>;

    //! export type of tuple of pointers
    using TupleType = typename MDTraits::template Tuple<PtrType>;

    /*!
     * \brief The default constructor
     */
    MultiDomainFVProblem() = default;

    /*!
     * \brief Contruct the problem
     * \param gridGeometries a tuple of grid geometry shared pointers
     */
    MultiDomainFVProblem(GridGeometries&& gridGeometries)
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<numSubDomains>{}, [&](auto&& id)
        {
            constexpr auto i = std::decay_t<decltype(id)>::value;
            elementAt(problems_, id) = std::make_shared<Type<i>>(std::get<i>(gridGeometries));
        });
    }

    /*!
     * \brief Applies the initial solution for all degrees of freedom of the grid.
     * \param sol the initial solution vector
     */
    void applyInitialSolution(SolutionVector& sol) const
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<numSubDomains>{}, [&](auto&& id)
        {
            elementAt(problems_, id)->applyInitialSolution(sol[id]);
        });
    }

    //! return the problem for domain with index i
    template<std::size_t i>
    const Type<i>& operator[] (Dune::index_constant<i> id) const
    { return *Dune::Hybrid::elementAt(problems_, id); }

    //! return the problem for domain with index i
    template<std::size_t i>
    Type<i>& operator[] (Dune::index_constant<i> id)
    { return *Dune::Hybrid::elementAt(problems_, id); }

    //! return the problem for domain with index i
    template<std::size_t i>
    PtrType<i> get(Dune::index_constant<i> id = Dune::index_constant<i>{})
    { return Dune::Hybrid::elementAt(problems_, id); }

    //! set the pointer for sub domain i
    template<std::size_t i>
    void set(PtrType<i> p, Dune::index_constant<i> id = Dune::index_constant<i>{})
    { Dune::Hybrid::elementAt(problems_, Dune::index_constant<i>{}) = p; }

    /*!
     * \brief return the grid variables tuple we are wrapping
     * \note the copy is not expensive since it is a tuple of shared pointers
     */
    TupleType getTuple()
    { return problems_; }

private:

    //! a tuple of pointes to all grid variables
    TupleType problems_;
};

} // end namespace Dumux

#endif
