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
 * \ingroup FEMDiscretization
 * \brief Grid geometry local view, which is a wrapper around a
 *        finite element basis local view.
 */
#ifndef DUMUX_DISCRETIZATION_FE_ELEMENT_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_FE_ELEMENT_GEOMETRY_HH

namespace Dumux {

/*!
 * \ingroup FEMDiscretization
 * \brief Grid geometry local view, which is a wrapper around a
 *        finite element basis local view.
 * \tparam The grid geometry type
 * \tparam The FEBasis local view
 */
template<class GridGeometry>
class FEElementGeometry
{
    using GridView = typename GridGeometry::GridView;

    using FEBasis = typename GridGeometry::FEBasis;
    using FEBasisLocalView = typename FEBasis::LocalView;

public:
    //! export type of the element
    using Element = typename GridView::template Codim<0>::Entity;

    //! constructor taking grid geometry
    FEElementGeometry(const GridGeometry& gg)
    : gridGeometry_(gg)
    , feBasisLocalView_(gg.feBasis().localView())
    {}

    //! Prepare element-local data
    void bind(const Element& element)
    {
        feBasisLocalView_.bind(element);
    }

    //! Prepare element-local data
    void bindElement(const Element& element)
    { bind(element); }

    //! Return the finite element basis local view
    const FEBasisLocalView& feBasisLocalView() const
    { return feBasisLocalView_; }

    //! Return reference to the grid geometry
    const GridGeometry& gridGeometry() const
    { return gridGeometry_; }

private:
    const GridGeometry& gridGeometry_;
    FEBasisLocalView feBasisLocalView_;
};

} // end namespace Dumux

#endif
