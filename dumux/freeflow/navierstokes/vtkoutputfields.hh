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
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesVtkOutputFields
 */
#ifndef DUMUX_NAVIER_STOKES_VTK_OUTPUT_FIELDS_HH
#define DUMUX_NAVIER_STOKES_VTK_OUTPUT_FIELDS_HH

#include <dune/common/fvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{

/*!
 * \ingroup NavierStokesModel
 * \brief Adds vtk output fields for the Navier-Stokes model
 */
template<class TypeTag>
class NavierStokesVtkOutputFields
{
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using FaceVariables = typename GET_PROP_TYPE(TypeTag, FaceVariables);

    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimensionworld>;

    // Helper type used for tag dispatching (to add discretization-specific fields).
    template<DiscretizationMethod method>
    using MethodType = std::integral_constant<DiscretizationMethod, method>;

public:
    //! Initialize the Navier-Stokes specific vtk output fields.
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.pressure(); }, "p");

        // add discretization-specific fields
        const auto discType = MethodType<GET_PROP_VALUE(TypeTag, DiscretizationMethod)>();
        additionalOutput_(vtk, discType);
    }

private:

    //! Adds discretization-specific fields (nothing by default).
    template <class VtkOutputModule, class AnyMethod>
    static void additionalOutput_(VtkOutputModule& vtk, AnyMethod)
    { }

    //! Adds discretization-specific fields (velocity vectors on the faces for the staggered discretization).
    template <class VtkOutputModule>
    static void additionalOutput_(VtkOutputModule& vtk, MethodType<DiscretizationMethod::staggered>)
    {
        const bool writeFaceVars = getParamFromGroup<bool>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Vtk.WriteFaceData", false);
        if(writeFaceVars)
        {
            auto faceVelocityVector = [](const SubControlVolumeFace& scvf, const FaceVariables& f)
                                      {
                                          GlobalPosition velocity(0.0);
                                          velocity[scvf.directionIndex()] = f.velocitySelf();
                                          return velocity;
                                      };

            vtk.addFaceVariable(faceVelocityVector, "faceVelocity");
        }
    }
};

} // end namespace Dumux

#endif
