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
 * \brief Multidomain wrapper for multiple grid geometries
 */
#ifndef DUMUX_MULTIDOMAIN_FVGRIDGEOMETRY_HH
#define DUMUX_MULTIDOMAIN_FVGRIDGEOMETRY_HH

#include <tuple>
#include <memory>
#include <utility>

#include <dune/common/hybridutilities.hh>
#include <dune/common/indices.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \brief A multidomain wrapper for multiple grid geometries
 * \tparam MDTraits The multidomain traits
 */
template<class MDTraits>
class MultiDomainFVGridGeometry
{
    static constexpr std::size_t numSubDomains = MDTraits::numSubDomains;

public:
    //! export base types of the stored type
    template<std::size_t i>
    using Type = typename MDTraits::template SubDomain<i>::GridGeometry;

    //! export pointer types the stored type
    template<std::size_t i>
    using PtrType = std::shared_ptr<Type<i>>;

    //! export type of tuple of pointers
    using TupleType = typename MDTraits::template Tuple<PtrType>;

    /*!
     * \brief The default constructor
     */
    MultiDomainFVGridGeometry() = default;

    /*!
     * \brief Contruct the problem
     * \param gridViews a tuple of gridViews
     */
    template<class GridViews>
    MultiDomainFVGridGeometry(GridViews&& gridViews)
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<numSubDomains>{}, [&](auto&& id)
        {
            constexpr auto i = std::decay_t<decltype(id)>::value;
            elementAt(gridGeometries_, id) = std::make_shared<Type<i>>(std::get<i>(gridViews));
        });
    }

    /*!
     * \brief Update all grid geometries (do this again after grid adaption)
     */
    void update()
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<numSubDomains>{}, [&](auto&& id)
        {
            elementAt(gridGeometries_, id)->update();
        });
    }

    //! return the grid geometry for domain with index i
    template<std::size_t i>
    const Type<i>& operator[] (Dune::index_constant<i> id) const
    { return *Dune::Hybrid::elementAt(gridGeometries_, id); }

    //! return the grid geometry for domain with index i
    template<std::size_t i>
    Type<i>& operator[] (Dune::index_constant<i> id)
    { return *Dune::Hybrid::elementAt(gridGeometries_, id); }

    ///! return the grid geometry pointer for domain with index i
    template<std::size_t i>
    PtrType<i> get(Dune::index_constant<i> id = Dune::index_constant<i>{})
    { return Dune::Hybrid::elementAt(gridGeometries_, id); }

    //! set the pointer for sub domain i
    template<std::size_t i>
    void set(PtrType<i> p, Dune::index_constant<i> id = Dune::index_constant<i>{})
    { Dune::Hybrid::elementAt(gridGeometries_, Dune::index_constant<i>{}) = p; }

    /*!
     * \brief return the grid variables tuple we are wrapping
     * \note the copy is not expensive since it is a tuple of shared pointers
     */
    TupleType getTuple()
    { return gridGeometries_; }

private:

    //! a tuple of pointes to all grid variables
    TupleType gridGeometries_;
};

} // end namespace Dumux

#endif
