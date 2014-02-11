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
#ifndef DUMUX_FVTRANSPORT2P2C_ADAPTIVE_HH
#define DUMUX_FVTRANSPORT2P2C_ADAPTIVE_HH

#include <dune/common/deprecated.hh>
#warning This class is deprecated, please use FV2dTransport2P2CAdaptive instead
#include <dumux/decoupled/2p2c/fv2dtransport2p2cadaptive.hh>

/**
 * @file
 * @brief  Finite Volume discretization of the component transport equation
 * @author Markus Wolff, Jochen Fritz, Benjamin Faigle
 */

namespace Dumux
{
//! Old naming of the new class FV2dTransport2P2CAdaptive
/*! \ingroup Adaptive2p2c
 *  Please use the proper class name that distinguishes between the
 *  adaptive implementation in 2d and in 3d !!
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class FVTransport2P2CAdaptive : public FV2dTransport2P2CAdaptive<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    //! Deprecated !! build a FV2dPressure2P2CAdaptive object instead
    /**
     * \param problem a problem class object
     */
public:
    FVTransport2P2CAdaptive(Problem& problem) : FV2dTransport2P2CAdaptive<TypeTag>(problem)
    {    }

    //! This method only exists to ensure a deprecation message is printed
    DUNE_DEPRECATED_MSG("This CLASS is deprecated! Use 2d-specific implementation FV2dTransport2P2CAdaptive instead")
    void initialize()
    {
        FV2dTransport2P2CAdaptive<TypeTag>::initialize();
    }

};
}
#endif
