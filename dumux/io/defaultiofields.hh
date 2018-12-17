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
 * \ingroup InputOutput
 * \brief Adds output fields to a given output module, this is the default if a
          model doesn't implement this functionality
 */
#ifndef DUMUX_IO_DEFAULT_IO_FIELDS_HH
#define DUMUX_IO_DEFAULT_IO_FIELDS_HH

#include <dune/common/exceptions.hh>

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief Adds output fields to a given output module
 */
class DefaultIOFields
{
public:
    template<class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        DUNE_THROW(Dune::NotImplemented, "This model doesn't implement default output fields!");
    }

    template <class FluidSystem = void, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx = 0, int state = 0)
    {
        DUNE_THROW(Dune::NotImplemented, "This model doesn't implement primaryVariableName!");
    }
};

} // end namespace Dumux

#endif
