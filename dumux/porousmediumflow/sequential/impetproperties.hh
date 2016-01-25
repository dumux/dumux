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
#ifndef DUMUX_IMPET_PROPERTIES_HH
#define DUMUX_IMPET_PROPERTIES_HH

#include <dumux/porousmediumflow/sequential/properties.hh>
#include <dumux/porousmediumflow/sequential/pressureproperties.hh>
#include <dumux/porousmediumflow/sequential/transportproperties.hh>

/*!
 * \ingroup IMPET
 * \ingroup IMPETProperties
 */
/*!
 * \file
 * \brief Base file for properties related to sequential IMPET algorithms
 */
namespace Dumux
{

template<class TypeTag>
class IMPET;

namespace Properties
{
/*!
 *
 * \brief General properties for sequential IMPET algorithms
 *
 * This class holds properties necessary for the sequential IMPET solution.
 */

//////////////////////////////////////////////////////////////////
// Type tags tags
//////////////////////////////////////////////////////////////////

//! The type tag for models based on the diffusion-scheme
NEW_TYPE_TAG(IMPET, INHERITS_FROM(DecoupledModel));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(ImpetCFLFactor);         //!< Scalar factor for additional scaling of the time step
NEW_PROP_TAG( ImpetSubCFLFactor );//!< Scalar factor for scaling of local sub-time-step
NEW_PROP_TAG(ImpetIterationFlag); //!< Flag to switch the iteration type of the IMPET scheme
NEW_PROP_TAG(ImpetIterationNumber); //!< Number of iterations if IMPET iterations are enabled by the IterationFlags
NEW_PROP_TAG(ImpetMaximumDefect); //!< Maximum Defect if IMPET iterations are enabled by the IterationFlags
NEW_PROP_TAG(ImpetRelaxationFactor); //!< Used for IMPET iterations
NEW_PROP_TAG(ImpetDtVariationRestrictionFactor);
NEW_PROP_TAG(ImpetPorosityThreshold);

//forward declaration!
NEW_PROP_TAG( Model );//! The model of the specific problem
}
}

#include <dumux/porousmediumflow/sequential/impet.hh>

namespace Dumux
{
namespace Properties
{
//set impet model
SET_TYPE_PROP(IMPET, Model, IMPET<TypeTag>);

//Set defaults
SET_SCALAR_PROP(IMPET, ImpetSubCFLFactor, 1.0);
SET_SCALAR_PROP(IMPET, ImpetCFLFactor, 1.0);
//! 0 = no iterations, 1 = iterate IterationNumber iterations, 2 = iterate until converged or IterationNumber is reached
SET_INT_PROP(IMPET, ImpetIterationFlag, 0);
SET_INT_PROP(IMPET, ImpetIterationNumber, 2);
SET_SCALAR_PROP(IMPET, ImpetMaximumDefect, 1e-5);
//! 1 = new solution is new solution, 0 = old solution is new solution
SET_SCALAR_PROP(IMPET, ImpetRelaxationFactor, 1.0);
SET_SCALAR_PROP(IMPET, ImpetDtVariationRestrictionFactor, std::numeric_limits<double>::max());
SET_SCALAR_PROP(IMPET, ImpetPorosityThreshold, 1e-6);
}
}

#endif
