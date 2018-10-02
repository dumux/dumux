// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup Common
 * \brief Boundary flag to store e.g. in sub control volume faces
 */
#ifndef DUMUX_BOUNDARY_FLAG_HH
#define DUMUX_BOUNDARY_FLAG_HH

#include <cstddef>

// ALUGrid specific includes
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

// dune-subgrid specific includes
#if HAVE_DUNE_SUBGRID
#include <dune/subgrid/subgrid.hh>
#endif


namespace Dumux {

/*!
 * \file
 * \ingroup Common
 * \brief Boundary flag to store e.g. in sub control volume faces
 * \tparam Grid the type of the gri
 */
template<class Grid>
class BoundaryFlag
{
public:
    BoundaryFlag() : flag_(-1) {}

    template<class Intersection>
    BoundaryFlag(const Intersection& i) : flag_(-1)
    {
        if (i.boundary())
            flag_ = i.boundarySegmentIndex();
    }

    using value_type = std::size_t;

    value_type get() const { return flag_; }

private:
    value_type flag_;
};


#if HAVE_DUNE_ALUGRID && DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
//! alu uses boundary id
template<int dim, int dimworld, Dune::ALUGridElementType elType, Dune::ALUGridRefinementType refinementType>
class BoundaryFlag<Dune::ALUGrid<dim, dimworld, elType, refinementType>>
{
public:
    BoundaryFlag() : flag_(-1) {}

    template<class Intersection>
    BoundaryFlag(const Intersection& i) : flag_(-1)
    {
        if (i.boundary())
            flag_ = i.impl().boundaryId();
    }

    using value_type = int;

    value_type get() const { return flag_; }

private:
    int flag_;
};
#endif

#if HAVE_DUNE_SUBGRID
//! dune-subgrid doesn't have this implemented
template<int dim, class HostGrid>
class BoundaryFlag<Dune::SubGrid<dim, HostGrid>>
{
public:
    BoundaryFlag() : flag_(-1) {}

    template<class Intersection>
    BoundaryFlag(const Intersection& i) : flag_(-1) {}

    using value_type = int;

    value_type get() const
    { DUNE_THROW(Dune::NotImplemented, "Sub-grid doesn't implement boundary segment indices!"); }

private:
    int flag_;
};
#endif

}  // end namespace Dumux

#endif
