// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Markus Wolff                                      *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#ifndef DUMUX_VELOCITYDEFAULT_HH
#define DUMUX_VELOCITYDEFAULT_HH

/**
 * @file
 * @brief  Default implementation of velocity class.
 * @author Markus Wolff
 */

#include <dune/grid/common/gridenums.hh>

#include <dumux/decoupled/common/decoupledproperties.hh>

namespace Dumux
{
//! \ingroup FV2p
//! \brief Default implementation of velocity class.
/*! Provides functions to make transport part compile
 * \tparam TypeTag The Type Tag
 */

template<class TypeTag>
class VelocityDefault
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;
    typedef typename GridView::Intersection Intersection;


public:
    //! Constructs a FVVelocityDefault object
    /*!
     * \param problem a problem class object
     */
    VelocityDefault(Problem& problem)
    {}

    //! Empty functions!
    /*!
     */
    void calculateVelocity(const Intersection& intersection, CellData& cellData)
    {}

    void calculateVelocityOnBoundary(const Intersection& intersection, CellData& cellData)
    {}

    bool calculateVelocityInTransport()
    {
        return false;
    }

    //! \brief Write data files
    /*  \param name file name */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {}


};
}
#endif
