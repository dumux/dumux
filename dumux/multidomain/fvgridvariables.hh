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
 * \brief Multidomain wrapper for multiple grid variables
 */
#ifndef DUMUX_MULTIDOMAIN_FV_GRID_VARIABLES_HH
#define DUMUX_MULTIDOMAIN_FV_GRID_VARIABLES_HH

#include <tuple>
#include <memory>
#include <utility>

#include <dune/common/hybridutilities.hh>
#include <dune/common/indices.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \brief A multidomain wrapper for multiple grid variables
 * \tparam MDTraits the multidomain traits
 */
template<class MDTraits>
class MultiDomainFVGridVariables
{
    using SolutionVector = typename MDTraits::SolutionVector;
    static constexpr std::size_t numSubDomains = MDTraits::numSubDomains;

    template<std::size_t i>
    using GridGeometry = typename MDTraits::template SubDomain<i>::GridGeometry;
    using GridGeometries = typename MDTraits::template TupleOfSharedPtrConst<GridGeometry>;

    template<std::size_t i>
    using Problem = typename MDTraits::template SubDomain<i>::Problem;
    using Problems = typename MDTraits::template TupleOfSharedPtrConst<Problem>;

public:
    //! export base types of the stored type
    template<std::size_t i>
    using Type = typename MDTraits::template SubDomain<i>::GridVariables;

    //! export pointer types the stored type
    template<std::size_t i>
    using PtrType = std::shared_ptr<Type<i>>;

    //! export type of tuple of pointers
    using TupleType = typename MDTraits::template Tuple<PtrType>;

    /*!
     * \brief The default constructor
     */
    MultiDomainFVGridVariables() = default;

    /*!
     * \brief Contruct the grid variables
     * \param gridGeometries a tuple of grid geometry shared pointers
     * \param problems a tuple of problem shared pointers
     */
    MultiDomainFVGridVariables(GridGeometries&& gridGeometries, Problems&& problems)
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<numSubDomains>{}, [&](auto&& id)
        {
            constexpr auto i = std::decay_t<decltype(id)>::value;
            elementAt(gridVars_, id) = std::make_shared<Type<i>>( std::get<i>(problems), std::get<i>(gridGeometries));
        });
    }

    //! initialize all variables
    void init(const SolutionVector& sol)
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<numSubDomains>{}, [&](auto&& id)
        {
            elementAt(gridVars_, id)->init(sol[id]);
        });
    }

    //! update all variables
    void update(const SolutionVector& sol, bool forceFluxCacheUpdate = false)
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<numSubDomains>{}, [&](auto&& id)
        {
            elementAt(gridVars_, id)->update(sol[id], forceFluxCacheUpdate);
        });
    }

    //! update all variables after grid adaption
    void updateAfterGridAdaption(const SolutionVector& sol)
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<numSubDomains>{}, [&](auto&& id)
        {
            elementAt(gridVars_, id)->updateAfterGridAdaption(sol[id]);
        });
    }

    /*!
     * \brief Sets the current state as the previous for next time step
     * \note this has to be called at the end of each time step
     */
    void advanceTimeStep()
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<numSubDomains>{}, [&](auto&& id)
        {
            elementAt(gridVars_, id)->advanceTimeStep();
        });
    }

    //! resets state to the one before time integration
    void resetTimeStep(const SolutionVector& sol)
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<numSubDomains>{}, [&](auto&& id)
        {
            elementAt(gridVars_, id)->resetTimeStep(sol[id]);
        });
    }

    //! return the grid variables for domain with index i
    template<std::size_t i>
    const Type<i>& operator[] (Dune::index_constant<i> id) const
    { return *Dune::Hybrid::elementAt(gridVars_, id); }

    //! return the grid variables for domain with index i
    template<std::size_t i>
    Type<i>& operator[] (Dune::index_constant<i> id)
    { return *Dune::Hybrid::elementAt(gridVars_, id); }

    //! return the grid variables tuple we are wrapping
    template<std::size_t i>
    PtrType<i> get(Dune::index_constant<i> id = Dune::index_constant<i>{})
    { return Dune::Hybrid::elementAt(gridVars_, id); }

    //! set the pointer for sub domain i
    template<std::size_t i>
    void set(PtrType<i> p, Dune::index_constant<i> id = Dune::index_constant<i>{})
    { Dune::Hybrid::elementAt(gridVars_, id) = p; }

    /*!
     * \brief return the grid variables tuple we are wrapping
     * \note the copy is not expensive since it is a tuple of shared pointers
     */
    TupleType getTuple()
    { return gridVars_; }

private:

    //! a tuple of pointes to all grid variables
    TupleType gridVars_;
};

} // end namespace Dumux

#endif
