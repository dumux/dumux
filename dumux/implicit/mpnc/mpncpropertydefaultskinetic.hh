/*****************************************************************************
 *   Copyright (C) 2010-2011 by Philipp Nuske                                *
 *   Copyright (C) 2011 by Andreas Lauser                                    *
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
/*!
 * \brief This file declares and defines the properties required by
 *        the kinetic modules the M-phase N-component model.
 */
#ifndef DUMUX_MPNC_PROPERTY_DEFAULTS_KINETIC_HH
#define DUMUX_MPNC_PROPERTY_DEFAULTS_KINETIC_HH

#include <dumux/implicit/mpnc/mpncproperties.hh>

#include <dumux/implicit/mpnc/mpncmodelkinetic.hh>

// interfacial area volume variables
#include "mpncvolumevariablesiakinetic.hh"

// kinetic mass module
#include "mass/mpncindicesmasskinetic.hh"
#include "mass/mpnclocalresidualmasskinetic.hh"
#include "mass/mpncvolumevariablesmasskinetic.hh"
#include "mass/mpncvtkwritermasskinetic.hh"

// kinetic energy module
#include "energy/mpncindicesenergykinetic.hh"
#include "energy/mpnclocalresidualenergykinetic.hh"
#include "energy/mpncfluxvariablesenergykinetic.hh"
#include "energy/mpncvolumevariablesenergykinetic.hh"
#include "energy/mpncvtkwriterenergykinetic.hh"


namespace Dumux
{
namespace Properties
{

/*!
 * \brief Set the property for the model.
 */
SET_TYPE_PROP(BoxMPNCKinetic, Model, MPNCModelKinetic<TypeTag>);

/*!
 * \brief Set the property for the awn surface by retrieving it from
 *        the spatial parameters.
 */
SET_PROP(BoxMPNCKinetic, AwnSurface)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;

public:
    typedef typename SpatialParams::AwnSurface type;
};

/*!
 * \brief Set the property for the awn surface by retrieving it from
 *        the spatial parameters.
 */
SET_PROP(BoxMPNCKinetic, AnsSurface)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;

public:
    typedef typename SpatialParams::AnsSurface type;
};

/*!
 * \brief Set the property for the awn surface by retrieving it from
 *        the spatial parameters.
 */
SET_PROP(BoxMPNCKinetic, AwsSurface)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;

public:
    typedef typename SpatialParams::AwsSurface type;
};

/*!
 * \brief Set the property for the awn surface parameters by extracting
 *        it from awn surface??.
 */
SET_PROP(BoxMPNCKinetic, AwnSurfaceParams)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, AwnSurface) AwnSurface;

public:
    typedef typename AwnSurface::Params type;
};

/*!
 * \brief Set the property for the awn surface parameters by extracting
 *        it from awn surface??.
 */
SET_PROP(BoxMPNCKinetic, AwsSurfaceParams)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, AwsSurface) AwsSurface;

public:
    typedef typename AwsSurface::Params type;
};

/*!
 * \brief Set the property for the awn surface parameters by extracting
 *        it from awn surface??.
 */
SET_PROP(BoxMPNCKinetic, AnsSurfaceParams)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, AnsSurface) AnsSurface;

public:
    typedef typename AnsSurface::Params type;
};

SET_BOOL_PROP(BoxMPNCKinetic, VelocityAveragingInModel, true);
}
}

#endif
