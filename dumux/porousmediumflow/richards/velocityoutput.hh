// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup RichardsModel
 * \brief Velocity output for the Richards model.
 */

#ifndef DUMUX_POROUSMEDIUMFLOW_RICHARDS_VELOCITYOUTPUT_HH
#define DUMUX_POROUSMEDIUMFLOW_RICHARDS_VELOCITYOUTPUT_HH

#include <dumux/porousmediumflow/velocityoutput.hh>

namespace Dumux {

/*!
 * \ingroup RichardsModel
 * \brief Velocity output policy for the Richards model.
 */
template<class GridVariables, class FluxVariables>
class RichardsVelocityOutput : public PorousMediumFlowVelocityOutput<GridVariables, FluxVariables>
{
    using ParentType = PorousMediumFlowVelocityOutput<GridVariables, FluxVariables>;

public:
    using ParentType::ParentType;

    //! Returns the number of phases.
    int numFluidPhases() const override { return 1; }

};

} // end namespace Dumux

#endif
