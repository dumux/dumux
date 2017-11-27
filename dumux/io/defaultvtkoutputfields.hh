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
 * \brief Adds vtk output fields specific to a model, this is the default if a
          model doesn't implement this functionality
 */
#ifndef DUMUX_DEFAULT_VTK_OUTPUT_FIELDS_HH
#define DUMUX_DEFAULT_VTK_OUTPUT_FIELDS_HH

#include <dune/common/exceptions.hh>

namespace Dumux
{

/*!
 * \ingroup Common, InputOutput
 * \brief Adds vtk output fields specific to a model
 */
class DefaultVtkOutputFields
{
public:
    template<class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        DUNE_THROW(Dune::NotImplemented, "This model doesn't implement default vtk fields!");
    }
};

} // end namespace Dumux

#endif
