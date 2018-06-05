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
 * \ingroup FreeflowNIModel
 * \copydoc Dumux::FreeflowNonIsothermalVtkOutputFields
 */
#ifndef DUMUX_FREEFLOW_NI_OUTPUT_FIELDS_HH
#define DUMUX_FREEFLOW_NI_OUTPUT_FIELDS_HH

namespace Dumux
{

/*!
 * \ingroup FreeflowNIModel
 * \brief Adds vtk output fields specific to non-isothermal free-flow models
 */
template<class IsothermalVtkOutputFields, class ModelTraits>
class FreeflowNonIsothermalVtkOutputFields
{

public:
    //! Initialize the non-isothermal specific vtk output fields.
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        IsothermalVtkOutputFields::init(vtk);
        add(vtk);
    }

    //! Add the non-isothermal specific vtk output fields.
    template <class VtkOutputModule>
    static void add(VtkOutputModule& vtk)
    {
        vtk.addVolumeVariable([](const auto& v){ return v.temperature(); }, "temperature");
        vtk.addVolumeVariable([](const auto& v){ return v.thermalConductivity(); }, "lambda");
        if (ModelTraits::usesTurbulenceModel())
            vtk.addVolumeVariable([](const auto& v){ return v.effectiveThermalConductivity() - v.thermalConductivity(); }, "lambda_t");
    }
};

} // end namespace Dumux

#endif
