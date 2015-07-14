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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/

/*!
 * \file
 * \brief Defines the indices required for the finite volume in the
 *        two phase discrete fracture-matrix model.
 */
#ifndef DUMUX_MODELS_2PDFM_INDICES_HH
#define DUMUX_MODELS_2PDFM_INDICES_HH

#include<dumux/implicit/2p/2pindices.hh>

namespace Dumux
{
// \{

/*!
 * \ingroup TwoPDFMModel
 * \ingroup ImplicitIndices
 * \brief The common indices for the \f$p_w-S_n\f$ formulation of the
 *        isothermal two-phase discrete fracture-matrix model.
 *
 * \tparam TypeTag The problem type tag
 * \tparam formulation The formulation, either pwsn or pnsw
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag,
          int formulation = TwoPFormulation::pwsn,
          int PVOffset = 0>
struct TwoPDFMIndices : public TwoPIndices <TypeTag, formulation, PVOffset>
{
    // Primary variable indices
    static const int pressureIdx = PVOffset + 0; //!< Index for wetting/non-wetting phase pressure (depending on formulation) in a solution vector
    static const int saturationIdx = PVOffset + 1; //!< Index of the saturation of the non-wetting/wetting phase

    // indices of the primary variables
    static const int pwIdx = PVOffset + 0; //!< index of the wetting phase pressure
    static const int snIdx = PVOffset + 1; //!< index of the nonwetting phase saturation

    // indices of the equations
    static const int contiWEqIdx = PVOffset + 0; //!< Index of the continuity equation of the wetting phase
    static const int contiNEqIdx = PVOffset + 1; //!< Index of the continuity equation of the non-wetting phase
};

// \}
} // namespace Dumux

#endif // DUMUX_MODELS_2PDFM_INDICES_HH
