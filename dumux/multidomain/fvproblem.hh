// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \brief Multidomain wrapper for multiple problems
 */
#ifndef DUMUX_MULTIDOMAIN_FV_PROBLEM_HH
#define DUMUX_MULTIDOMAIN_FV_PROBLEM_HH

#include <tuple>
#include <memory>
#include <utility>

#include <dune/common/hybridutilities.hh>
#include <dune/common/indices.hh>

#include <dumux/multidomain/fvgridgeometry.hh>

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
     * \brief Construct the problem
     * \param gridGeometries a tuple of grid geometry shared pointers
     */
    MultiDomainFVProblem(MultiDomainFVGridGeometry<MDTraits> gridGeometries)
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<numSubDomains>{}, [&](auto&& id)
        {
            constexpr auto i = std::decay_t<decltype(id)>::value;
            std::get<i>(problems_) = std::make_shared<Type<i>>(gridGeometries.template get<i>());
        });
    }

    /*!
     * \brief Construct wrapper from a tuple of problems
     * \param problemTuple a tuple of shared_ptrs to the problems
     */
    MultiDomainFVProblem(TupleType problemTuple)
    : problems_(std::move(problemTuple))
    {}

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
    { return *std::get<i>(problems_); }

    //! return the problem for domain with index i
    template<std::size_t i>
    Type<i>& operator[] (Dune::index_constant<i> id)
    { return *std::get<i>(problems_); }

    //! access the problem ptr for domain with index i
    template<std::size_t i>
    const PtrType<i>& get(Dune::index_constant<i> id = Dune::index_constant<i>{}) const
    { return std::get<i>(problems_); }

    //! access the problem ptr for domain with index i
    template<std::size_t i>
    PtrType<i>& get(Dune::index_constant<i> id = Dune::index_constant<i>{})
    { return std::get<i>(problems_); }

    /*!
     * \brief Access the underlying tuple representation
     */
    TupleType& asTuple()
    { return problems_; }

    /*!
     * \brief Access the underlying tuple representation
     */
    const TupleType& asTuple() const
    { return problems_; }

private:

    //! a tuple of points to all grid variables
    TupleType problems_;
};

} // end namespace Dumux

#endif
