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
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the regularized version of the van Genuchten's
 *        capillary pressure / relative permeability  <-> saturation relation.
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_TWOP_DEFAULT_REGULARIZATION_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_TWOP_DEFAULT_REGULARIZATION_HH

namespace Dumux {
namespace FluidMatrix {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief The default regularization policy
 */
template<class Scalar>
class TwoPDefaultRegularization
{
public:
    //! Init does nothing
    template<typename... Args> void init(Args&&...) {}
    //! No parameters
    template<class S> struct Params {};
};

} // end namespace FluidMatrix
} // end namespace Dumux

#endif
