// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup GeomechanicsModels
 * \brief Velocity output for geomechanical models
 */
#ifndef DUMUX_GEOMECHANICS_VELOCITYOUTPUT_HH
#define DUMUX_GEOMECHANICS_VELOCITYOUTPUT_HH


#include <dune/common/exceptions.hh>
#include <dumux/io/velocityoutput.hh>

namespace Dumux {

/*!
 * \ingroup GeomechanicsModels
 * \brief Velocity output for geomechanical models.
 *        This class could be used to compute the temporal derivative
 *        of the displacement. Currently this is not implemented and
 *        we simply define this here in order to be able to reuse the
 *        VtkOutputModule which expects a VelocityOutput class.
 */
template<class GridVariables>
class GeomechanicsVelocityOutput
: public VelocityOutput<GridVariables>
{
public:
    //! The constructor
    template< typename... Args >
    GeomechanicsVelocityOutput(Args&&... args)
    { DUNE_THROW(Dune::NotImplemented, "Velocity output for geomechanical models."); }

    //! Output is currently disabled (not implemented)
    bool enableOutput() const override { return false; }
};

} // end namespace Dumux

#endif
