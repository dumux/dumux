// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
