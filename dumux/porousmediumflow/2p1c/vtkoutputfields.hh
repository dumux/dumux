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
 * \ingroup TwoPOneCModel
 * \copydoc Dumux::TwoPOneCVtkOutputFields
 */
#ifndef DUMUX_TWOP_ONEC_VTK_OUTPUT_FIELDS_HH
#define DUMUX_TWOP_ONEC_VTK_OUTPUT_FIELDS_HH

#include <dumux/porousmediumflow/2p/vtkoutputfields.hh>

namespace Dumux {

/*!
 * \ingroup TwoPOneCModel
 * \brief Adds vtk output fields specific to two-phase one-component model.
 */
class TwoPOneCVtkOutputFields
{
public:
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        // use default fields from the 2p model
        TwoPVtkOutputFields::init(vtk);

        // output additional to TwoP output:
        vtk.addVolumeVariable([](const auto& v){ return v.priVars().state(); }, "phase presence");
    }
};

} // end namespace Dumux

#endif
