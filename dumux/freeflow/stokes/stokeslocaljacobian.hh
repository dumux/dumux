// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief Caculates the jacobian of models based on the box scheme element-wise.
 */
#ifndef DUMUX_STOKES_LOCAL_JACOBIAN_HH
#define DUMUX_STOKES_LOCAL_JACOBIAN_HH

#include "stokesproperties.hh"
#include <dumux/boxmodels/common/boxlocaljacobian.hh>
#include <dune/istl/matrix.hh>
#include <dumux/boxmodels/common/boxlocaljacobian.hh>

namespace Dumux
{
/*!
 * \ingroup BoxStokesModel
 * \brief Element-wise calculation of the jacobian matrix for models
 *        based on the box scheme.
 *
 * \todo Please doc me more!
 */
template<class TypeTag>
class StokesLocalJacobian : public BoxLocalJacobian<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    //! \copydoc BoxLocalJacobian::numericEpsilon()
    Scalar numericEpsilon(int scvIdx,
                          int pvIdx) const
    {
        Scalar pv = this->curVolVars_[scvIdx].primaryVars()[pvIdx];
        if (pvIdx == 0 || pvIdx == 1){
            return 1e-7*(std::abs(pv) + 1);
        }
        return 1e-9*(std::abs(pv) + 1);
    }
};
}

#endif
