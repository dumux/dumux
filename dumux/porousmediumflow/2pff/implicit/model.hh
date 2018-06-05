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
 * \ingroup TwoPFractionalFlowModel
 * \brief Adaption of the fully implicit scheme to the two-phase flow model in fractional flow formulation
 */

#ifndef DUMUX_TWOP_FF_MODEL_HH
#define DUMUX_TWOP_FF_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/2p/model.hh>

#include "localresidual.hh"
#include "darcyslaw.hh"
#include "upwindscheme.hh"

namespace Dumux {

template<TwoPFormulation formulation>
struct TwoPFractionalFlowModelTraits
{
    using Indices = TwoPIndices;

    static constexpr TwoPFormulation priVarFormulation()
    { return formulation; }

    static constexpr int numEq() { return 2; } // vt as variable
    static constexpr int numPhases() { return 2; }
    static constexpr int numComponents() { return 2; }

    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return false; }
    static constexpr bool enableEnergyBalance() { return false; }
};

namespace Properties {

//////////////////////////////////////////////////////////////////

//! The type tag for the isothermal two-phase model
NEW_TYPE_TAG(TwoPFractionalFlow, INHERITS_FROM(TwoP));

//! We have a special local resdiual
SET_TYPE_PROP(TwoPFractionalFlow, LocalResidual, TwoPFractionalFlowLocalResidual<TypeTag>);

//! The model traits class
SET_TYPE_PROP(TwoPFractionalFlow, ModelTraits, TwoPFractionalFlowModelTraits<GET_PROP_VALUE(TypeTag, Formulation)>);

//! We have a special advection type
SET_TYPE_PROP(TwoPFractionalFlow, AdvectionType, FractionalFlowCCTpfaDarcysLaw<typename GET_PROP_TYPE(TypeTag, Scalar),
                                                                               typename GET_PROP_TYPE(TypeTag, FVGridGeometry)>);

//! We have a special upwind scheme
SET_TYPE_PROP(TwoPFractionalFlow, FluxVariables, PorousMediumFluxVariables<TypeTag, TwoPFractionalFlowUpwindScheme>);
} // end namespace Properties
} // end namespace Dumux

#endif
