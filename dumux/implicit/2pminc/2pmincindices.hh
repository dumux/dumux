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
 * \brief Defines the indices required for the two-phase minc fully implicit model.
 */
#ifndef DUMUX_BOX_2PMINC_INDICES_HH
#define DUMUX_BOX_2PMINC_INDICES_HH

#include "2pmincproperties.hh"

namespace Dumux
{
// \{


/*!
 * \ingroup TwoPMincModel
 * \ingroup ImplicitIndices
 * \brief Enumerates the formulations which the two-phase model accepts.
 */
struct TwoPMincFormulation
{
    static const int pwsn = 0; //!< pw and sn as primary variables
    static const int pnsw = 1; //!< pn and sw as primary variables
};

/*!
 * \ingroup TwoPMincModel
 * \ingroup ImplicitIndices
 * \brief Defines the indices required for the two-phase fully implicit model.
 *
 * \tparam TypeTag The problem type tag
 */
template <class TypeTag>
struct TwoPMincCommonIndices
{
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    // Phase indices
    static const int wPhaseIdx = FluidSystem::wPhaseIdx; //!< Index of the wetting phase
    static const int nPhaseIdx = FluidSystem::nPhaseIdx; //!< Index of the non-wetting phase
};

/*!
 * \ingroup TwoPMincBoxModel
 * \ingroup ImplicitIndices
 * \brief The indices for the \f$p_w-S_n\f$ formulation of the
 *        isothermal two-phase model.
 *
 * \tparam TypeTag The problem type tag
 * \tparam formulation The formulation, either pwsn or pnsw
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag, 
          int formulation = TwoPMincFormulation::pwsn, 
          int PVOffset = 0 >
struct TwoPMincIndices 
: public TwoPMincCommonIndices<TypeTag>, TwoPMincFormulation
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
    
    static const int pIdxc(int numC) {return pwIdx + 2 *numC;} //!< index of the wetting phase pressure for continuum numC
    static const int sIdxc(int numC) {return snIdx + 2 *numC;} //!< index of the non-wetting phase saturation for continuum numC
    
    static const int contiWEqIdxc(int numC) {return contiWEqIdx + 2 *numC;} //!< Index of the continuity equation of the wetting phase for continuum numC
    static const int contiNEqIdxc(int numC) {return contiNEqIdx + 2 *numC;} //!< Index of the continuity equation of the non-wetting phase for continuum numC
};

/*!
 * \ingroup TwoPMincModel
 * \ingroup ImplicitIndices
 * \brief The indices for the \f$p_n-S_w\f$ formulation of the
 *        isothermal two-phase model.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag, int PVOffset>
struct TwoPMincIndices<TypeTag, TwoPMincFormulation::pnsw, PVOffset>
: public TwoPMincCommonIndices<TypeTag>, TwoPMincFormulation
{
    // Primary variable indices
    static const int pressureIdx = PVOffset + 0; //!< Index for wetting/non-wetting phase pressure (depending on formulation) in a solution vector
    static const int saturationIdx = PVOffset + 1; //!< Index of the saturation of the non-wetting/wetting phase

    // indices of the primary variables
    static const int pnIdx = PVOffset + 0; //!< index of the nonwetting phase pressure
    static const int swIdx = PVOffset + 1; //!< index of the wetting phase saturation

    // indices of the equations
    static const int contiNEqIdx = PVOffset + 0; //!< Index of the continuity equation of the non-wetting phase
    static const int contiWEqIdx = PVOffset + 1; //!< Index of the continuity equation of the wetting phase

    static const int pIdxc(int numC) {return pnIdx + 2 *numC;} //!< index of the nonwetting phase pressure for continuum numC
    static const int sIdxc(int numC) {return swIdx + 2 *numC;} //!< index of the wetting phase saturation for continuum numC
     
    static const int contiWEqIdxc(int numC) {return contiWEqIdx + 2 *numC;} //!< Index of the continuity equation of the wetting phase for continuum numC
    static const int contiNEqIdxc(int numC) {return contiNEqIdx + 2 *numC;} //!< Index of the continuity equation of the non-wetting phase for continuum numC
};

// \}
} // namespace Dumux


#endif
