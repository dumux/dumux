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
 * \brief Element-wise calculation of the Jacobian matrix for the Stokes models.
 */
#ifndef DUMUX_STOKES_LOCAL_JACOBIAN_HH
#define DUMUX_STOKES_LOCAL_JACOBIAN_HH

#include <dune/istl/matrix.hh>
#include "stokesproperties.hh"
#include <dumux/implicit/localjacobian.hh>

namespace Dumux
{
/*!
 * \ingroup BoxStokesModel
 * \ingroup ImplicitLocalJacobian
 * \brief Element-wise calculation of the Jacobian matrix for the Stokes models.
 *
 * This overloads the numericEpsilon method of the ImplicitLocalJacobian.
 * The momentum balance equation uses larger epsilons than the rest.
 */
template<class TypeTag>
class StokesLocalJacobian : public ImplicitLocalJacobian<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

public:
    //! \copydoc ImplicitLocalJacobian::numericEpsilon()
    Scalar numericEpsilon(const int scvIdx,
                          const int pvIdx) const
    {
        Scalar pv = this->curVolVars_[scvIdx].priVars()[pvIdx];
        if (pvIdx < GridView::dimension){
            return 1e-7*(std::abs(pv) + 1);
        }
        return 1e-9*(std::abs(pv) + 1);
    }
};
}

#endif
