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
 * \ingroup SweModel
 * \copydoc Dumux::SweVtkOutputFields
 */
#ifndef DUMUX_SWE_VTK_OUTPUT_FIELDS_HH
#define DUMUX_SWE_VTK_OUTPUT_FIELDS_HH

#include <dune/common/fvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{

/*!
 * \ingroup SweModel
 * \brief Adds vtk output fields for the SWEs model
 */
template<class TypeTag>
class SweVtkOutputFields
{
    //using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);

    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimensionworld>;


public:
    //! Initialize the Navier-Stokes specific vtk output fields.
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.getH(); }, "h");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.getU(); }, "u");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.getV(); }, "v");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.getBottom(); }, "bottom");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.getBottom() + v.getH(); }, "theta");
    }

};

} // end namespace Dumux

#endif
