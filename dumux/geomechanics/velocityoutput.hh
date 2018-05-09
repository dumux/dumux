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
 * \brief Velocity output for geomechanical models
 */
#ifndef DUMUX_GEOMECHANICS_VELOCITYOUTPUT_HH
#define DUMUX_GEOMECHANICS_VELOCITYOUTPUT_HH

namespace Dumux {

/*!
 * \brief Velocity output for geomechanical models.
 *        This class could be used to compute the temporal derivative
 *        of the displacement. Currently this is not implemented and
 *        we simply define this here in order to be able to reuse the
 *        VtkOutputModule which expects a VelocityOutput class.
 */
class GeomechanicsVelocityOutput
{
public:
    //! The constructor
    template< typename... Args >
    GeomechanicsVelocityOutput(Args&&... args) {}

    //! Output is currently disabled (not implemented)
    static constexpr bool enableOutput() { return false; }

    //! There is always only one solid phase
    static constexpr int numPhaseVelocities() { return 1; }

    //! Returns the name of phase for which velocity is computed
    static std::string phaseName(int phaseIdx)
    { DUNE_THROW(Dune::NotImplemented, "Velocity output for geomechanical models."); }

    //! Calculate the velocities for the scvs in the element
    template< typename... Args >
    void calculateVelocity(Args... args)
    { DUNE_THROW(Dune::NotImplemented, "Velocity output for geomechanical models."); }
};

} // end namespace Dumux

#endif
