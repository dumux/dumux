// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
