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
 *
 * \brief Finite element map for P0 elements
 */

#ifndef DUMUX_P0FEM_HH
#define DUMUX_P0FEM_HH

#if HAVE_DUNE_PDELAB

#include <dune/geometry/type.hh>

#include <dune/localfunctions/lagrange/p0.hh>
#include <dune/pdelab/finiteelementmap/finiteelementmap.hh>

namespace Dumux {

/*!
 * \brief Finite element map for P0 elements
 */
template<class D, class R, int d, bool isCube = true>
class P0LocalFiniteElementMap
: public Dune::PDELab::SimpleLocalFiniteElementMap< Dune::P0LocalFiniteElement<D,R,d> >
{
    typedef Dune::PDELab::SimpleLocalFiniteElementMap< Dune::P0LocalFiniteElement<D,R,d> > ParentType;
public:
      P0LocalFiniteElementMap ()
      : ParentType (Dune::P0LocalFiniteElement<D,R,d>(isCube ? Dune::GeometryType(Dune::GeometryType::cube, d)
                                                             : Dune::GeometryType(Dune::GeometryType::simplex, d)))
      DUNE_DEPRECATED_MSG("should not be needed anymore")
      {}
};

} // namespace Dumux

#endif // HAVE_DUNE_PDELAB
#endif // DUMUX_P0FEM_HH
