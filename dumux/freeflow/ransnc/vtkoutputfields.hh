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
 * \ingroup RANSNCModel
 * \copydoc Dumux::RANSNCVtkOutputFields
 */
#ifndef DUMUX_RANS_NC_VTK_OUTPUT_FIELDS_HH
#define DUMUX_RANS_NC_VTK_OUTPUT_FIELDS_HH

#include <dumux/freeflow/rans/vtkoutputfields.hh>
#include <dumux/freeflow/navierstokesnc/vtkoutputfields.hh>

namespace Dumux
{

/*!
 * \ingroup RANSNCModel
 * \brief Adds vtk output fields specific to the RANSNC model
 */
template<class FVGridGeometry, class FluidSystem, int phaseIdx>
class RANSNCVtkOutputFields
{

public:
    //! Initialize the RANSNC specific vtk output fields.
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        RANSVtkOutputFields<FVGridGeometry>::init(vtk);
        NavierStokesNCVtkOutputFields<FVGridGeometry, FluidSystem, phaseIdx>::init(vtk);
        vtk.addVolumeVariable([](const auto& v){ return v.eddyDiffusivity(); }, "D_t");
    }
};

} // end namespace Dumux

#endif
