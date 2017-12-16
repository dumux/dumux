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
 * \brief This file declares and defines the properties required by
 *        the kinetic modules the M-phase N-component model.
 */
#ifndef DUMUX_MPNC_PROPERTY_DEFAULTS_KINETIC_HH
#define DUMUX_MPNC_PROPERTY_DEFAULTS_KINETIC_HH

#include <dumux/porousmediumflow/mpnc/implicit/properties.hh>

#include <dumux/porousmediumflow/mpnc/implicit/modelkinetic.hh>

// interfacial area volume variables
#include "volumevariablesiakinetic.hh"

// kinetic mass module
#include "mass/indiceskinetic.hh"
#include "mass/localresidualkinetic.hh"
#include "mass/volumevariableskinetic.hh"
#include "mass/vtkwriterkinetic.hh"

// kinetic energy module
#include "energy/indiceskinetic.hh"
#include "energy/localresidualkinetic.hh"
#include "energy/fluxvariableskinetic.hh"
#include "energy/volumevariableskinetic.hh"
#include "energy/vtkwriterkinetic.hh"

namespace Dumux
{
namespace Properties
{

/*!
 * \brief Set the property for the model.
 */
SET_TYPE_PROP(BoxMPNCKinetic, Model, MPNCModelKinetic<TypeTag>);

/*!
 * \brief Set the property for the awn surface parameters by extracting
 *        it from awn surface.
 */
SET_PROP(BoxMPNCKinetic, AwnSurfaceParams)
{
private:
    using AwnSurface = typename GET_PROP_TYPE(TypeTag, AwnSurface);

public:
    using type = typename AwnSurface::Params;
};

/*!
 * \brief Set the property for the aws surface parameters by extracting
 *        it from aws surface.
 */
SET_PROP(BoxMPNCKinetic, AwsSurfaceParams)
{
private:
    using AwsSurface = typename GET_PROP_TYPE(TypeTag, AwsSurface);

public:
    using type = typename AwsSurface::Params;
};

/*!
 * \brief Set the property for the ans surface parameters by extracting
 *        it from the surface.
 */
SET_PROP(BoxMPNCKinetic, AnsSurfaceParams)
{
private:
    using AnsSurface = typename GET_PROP_TYPE(TypeTag, AnsSurface);

public:
    using type = typename AnsSurface::Params;
};

SET_BOOL_PROP(BoxMPNCKinetic, VelocityAveragingInModel, true);

/*!
 * \brief Set the default formulation for the nusselt correlation
 *        Other possible parametrizations can be found in the dimensionlessnumbers
 */
SET_PROP(BoxMPNCKinetic, NusseltFormulation ){
    private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using DimLessNum = DimensionlessNumbers<Scalar>;
    public:
    static constexpr int value = DimLessNum::NusseltFormulation::WakaoKaguei;};

/*!
 * \brief Set the default formulation for the sherwood correlation
 *        Other possible parametrizations can be found in the dimensionlessnumbers
 */
SET_PROP(BoxMPNCKinetic, SherwoodFormulation ){
    private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using DimLessNum = DimensionlessNumbers<Scalar>;
    public:
    static constexpr int value = DimLessNum::SherwoodFormulation::WakaoKaguei;};

} // end namespace Properties
} // end namespace Dumux

#endif
