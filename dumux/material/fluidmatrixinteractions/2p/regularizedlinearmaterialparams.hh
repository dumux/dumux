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
 * \brief   Parameters that are necessary for the \em regularization of
 *          the linear constitutive relations.
 */
#ifndef REGULARIZED_LINEAR_PARAMS_HH
#define REGULARIZED_LINEAR_PARAMS_HH

#warning "This header is deprecated. Use SmoothedLinearLaw instead. Removal after 3.3"

#include "linearmaterialparams.hh"

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief   Parameters that are necessary for the \em regularization of
 *          the linear constitutive relations.
 */
template<class ScalarT>
class [[deprecated("Use new material laws! Removal after 3.3")]] RegularizedLinearMaterialParams : public LinearMaterialParams<ScalarT>
{
public:
    using Scalar = ScalarT;

    RegularizedLinearMaterialParams()
    {
        setKrLowS(0.05);
        setKrHighS(0.95);
    }

    /*!
     * \brief Set the threshold saturation respective phase below
     *        which the relative permeability gets regularized.
     */
    void setKrLowS(Scalar krLowS)
    {
        krLowS_ = krLowS;
    }

    /*!
     * \brief Return the threshold saturation respective phase below
     *        which the relative permeability gets regularized.
     */
    Scalar krLowS() const
    {
        return krLowS_;
    }

    /*!
     * \brief Set the threshold saturation of the respective phase
     *        above which the relative permeability gets regularized.
     */
    void setKrHighS(Scalar krHighS)
    {
        krHighS_ = krHighS;
    }

    /*!
     * \brief Return the threshold saturation of the respective phase
     *        above which the relative permeability gets regularized.
     */
    Scalar krHighS() const
    {
        return krHighS_;
    }

private:
    Scalar krLowS_;
    Scalar krHighS_;
};
} // namespace Dumux

#endif
