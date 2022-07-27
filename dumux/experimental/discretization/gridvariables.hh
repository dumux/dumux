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
 * \ingroup Experimental
 * \ingroup Discretization
 * \brief Base class for grid variables
 */
#ifndef DUMUX_DISCRETIZATION_GRID_VARIABLES_HH
#define DUMUX_DISCRETIZATION_GRID_VARIABLES_HH

#include <utility>
#include <memory>

#include <dumux/experimental/common/variables.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup Experimental
 * \ingroup Discretization
 * \brief Base class for grid variables.
 * \tparam GG The grid geometry type
 * \tparam X The type used for solution vectors
 */
template<class GG, class X>
class GridVariables
: public Variables<X>
{
    using ParentType = Variables<X>;

public:
    //! export the grid geometry type
    using GridGeometry = GG;

    /*!
     * \brief Constructor from a grid geometry. The remaining arguments must
     *        be valid arguments for the construction of the Variables class.
     */
    template<class... Args>
    GridVariables(std::shared_ptr<const GridGeometry> gridGeometry,
                  Args&&... args)
    : ParentType(std::forward<Args>(args)...)
    , gridGeometry_(gridGeometry)
    {}

    //! Return a reference to the grid geometry
    const GridGeometry& gridGeometry() const
    { return *gridGeometry_; }

private:
    std::shared_ptr<const GridGeometry> gridGeometry_;
};

} // end namespace Dumux::Experimental

#endif
