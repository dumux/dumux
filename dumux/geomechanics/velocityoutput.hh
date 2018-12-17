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
 * \ingroup Geomechanics
 * \brief Velocity output for geomechanical models
 */
#ifndef DUMUX_GEOMECHANICS_VELOCITYOUTPUT_HH
#define DUMUX_GEOMECHANICS_VELOCITYOUTPUT_HH


#include <dune/common/exceptions.hh>
#include <dumux/io/velocityoutput.hh>

namespace Dumux {

/*!
 * \ingroup Geomechanics
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
