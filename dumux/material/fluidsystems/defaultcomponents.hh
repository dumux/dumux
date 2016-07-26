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
 * \ingroup Fluidsystems
 *
 * \brief Provides defaults for the members of a fluidsystem.
 */
#ifndef DUMUX_DEFAULT_COMPONENTS_HH
#define DUMUX_DEFAULT_COMPONENTS_HH

#include <dumux/common/propertysystem.hh>

#include <dumux/common/basicproperties.hh>

#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/n2.hh>
#include <dumux/material/components/o2.hh>
#include <dumux/material/components/h2.hh>
#include <dumux/material/components/ch4.hh>
#include <dumux/material/components/simpleco2.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/brine.hh>
#include <dumux/material/components/tabulatedcomponent.hh>

#include <dune/common/stdstreams.hh>

namespace Dumux
{
namespace Properties
{
//! Defines the components which are being used by the fluid system by
//! default and how they are initialized
NEW_PROP_TAG(DefaultComponents);

//! Defines, if a detailed description for members of the fluidsystem is used
NEW_PROP_TAG(EnableComplicatedFluidSystem);

//! Defines the components which are actually being used by the fluidsystem
NEW_PROP_TAG(Components);

//! Specifies default component names and initializes the H2O fluid properties
SET_PROP(NumericModel, DefaultComponents)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dumux::H2O<Scalar> H2O_IAPWS;

public:
    typedef TabulatedComponent<Scalar, H2O_IAPWS> H2O;
    typedef Dumux::N2<Scalar> N2;
    typedef Dumux::O2<Scalar> O2;
    typedef Dumux::H2<Scalar> H2;
    typedef Dumux::CH4<Scalar> CH4;
    typedef Dumux::SimpleCO2<Scalar> SimpleCO2;
    typedef Dumux::SimpleH2O<Scalar> SimpleH2O;
    typedef Dumux::Brine<Scalar, Dumux::H2O<Scalar> > BrineRawComponent;
    typedef TabulatedComponent<Scalar, BrineRawComponent > Brine;

    static void init()
    {
        int nT = 100;
        int nP = 200;
        Dune::dinfo << "Initializing tables for the H2O fluid properties ("
                    << nT*nP
                    << " entries).\n";
        H2O::init(273.15, 623.15, nT, -10, 20e6, nP);
    }
};

//! Initialize the components with default behavior
SET_PROP(NumericModel, Components) : public GET_PROP(TypeTag, DefaultComponents) {};

/*!
 * \brief Enables a detailed description of the fluidsystem
 *
 * Complicated but detailed members of fluidsystems (e.g. phase viscosity,
 * phase density) can be simplified for efficiency reasons with this property.
 * Typically, such high demands on accuracy are not needed, so this property
 * is set to "false" as the default.
 */
SET_BOOL_PROP(NumericModel, EnableComplicatedFluidSystem, false);

} // namespace Properties
} // namespace Dumux

#endif
