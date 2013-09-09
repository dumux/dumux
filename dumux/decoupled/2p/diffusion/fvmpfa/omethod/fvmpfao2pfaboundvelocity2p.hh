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

#ifndef DUMUX_MPFAO2PFABOUNDVELOCITIES2P_HH
#define DUMUX_MPFAO2PFABOUNDVELOCITIES2P_HH

#warning This file is deprecated. Use fvmpfao2dpressurevelocity2p.hh instead.
#include "fvmpfao2dpressurevelocity2p.hh"

namespace Dumux
{
template<class TypeTag>
class FVMPFAO2PFABoundPressure2P: public FvMpfaO2dPressureVelocity2p<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    DUNE_DEPRECATED_MSG("use FvMpfaO2dPressureVelocity2p(problem) instead")
    FVMPFAO2PFABoundPressure2P(Problem& problem):
        FvMpfaO2dPressureVelocity2p<TypeTag>(problem)
        {}
};
}
#endif
