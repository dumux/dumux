// -**- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
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
 * \ingroup Properties
 * \ingroup ImplicitProperties
 * \ingroup TwoPNCMinModel
 * \file
 *
 * \brief Defines default values for most properties required by the
 *        two-phase n-component mineralization fully implicit model.
 */
#ifndef DUMUX_2PNCMIN_PROPERTY_DEFAULTS_HH
#define DUMUX_2PNCMIN_PROPERTY_DEFAULTS_HH

#include "2pncminindices.hh"
#include "2pncminmodel.hh"
#include "2pncminindices.hh"
#include "2pncminfluxvariables.hh"
#include "2pncminvolumevariables.hh"
#include "2pncminproperties.hh"

#include <dumux/implicit/2pnc/2pncnewtoncontroller.hh>
#include <dumux/implicit/common/implicitdarcyfluxvariables.hh>
#include <dumux/material/spatialparams/implicitspatialparams.hh>

namespace Dumux
{

namespace Properties {
//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////

/*!
 * \brief Set the property for the number of secondary components.
 * Secondary components are components calculated from
 * primary components by equilibrium relations and
 * do not have mass balance equation on their own.
 * These components are important in the context of bio-mineralization applications.
 * We just forward the number from the fluid system
 *
 */
SET_PROP(TwoPNCMin, NumSecComponents)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static const int value = FluidSystem::numSecComponents;

};
/*!
 * \brief Set the property for the number of solid phases, excluding the non-reactive matrix.
 *
 * We just forward the number from the fluid system
 *
 */
SET_PROP(TwoPNCMin, NumSPhases)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static const int value = FluidSystem::numSPhases;
};

/*!
 * \brief Set the property for the number of equations.
 * For each component and each precipitated mineral/solid phase one equation has to
 * be solved.
 */
SET_PROP(TwoPNCMin, NumEq)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static const int value = FluidSystem::numComponents + FluidSystem::numSPhases;
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(TwoPNCMin, FluidState){
    private:
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    public:
        typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> type;
};

//! Use the 2pncmin local residual operator
SET_TYPE_PROP(TwoPNCMin,
              LocalResidual,
              TwoPNCMinLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(TwoPNCMin, Model, TwoPNCMinModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(TwoPNCMin, VolumeVariables, TwoPNCMinVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(TwoPNCMin, FluxVariables, TwoPNCMinFluxVariables<TypeTag>);

//! The indices required by the isothermal 2pNcMin model
SET_TYPE_PROP(TwoPNCMin, Indices, TwoPNCMinIndices <TypeTag, /*PVOffset=*/0>);

//! disable useSalinity for the calculation of osmotic pressure by default
SET_BOOL_PROP(TwoPNCMin, useSalinity, false);


//! default value for the forchheimer coefficient
// Source: Ward, J.C. 1964 Turbulent flow in porous media. ASCE J. Hydraul. Div 90.
//        Actually the Forchheimer coefficient is also a function of the dimensions of the
//        porous medium. Taking it as a constant is only a first approximation
//        (Nield, Bejan, Convection in porous media, 2006, p. 10)
SET_SCALAR_PROP(TwoPNCMin, SpatialParamsForchCoeff, 0.55);

}

}

#endif
