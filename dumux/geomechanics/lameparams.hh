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
 * \brief \copydoc Dumux::LameParams
 */
#ifndef DUMUX_GEOMECHANICS_LAME_PARAMS_HH
#define DUMUX_GEOMECHANICS_LAME_PARAMS_HH

namespace Dumux {

/*!
 * \ingroup Geomechanics
 * \brief Structure encapsulating the lame parameters
 */
template<class Scalar>
struct LameParams
{
    //! Default constructor
    LameParams() = default;

    //! Constructor taking lambda and mu directly
    LameParams(Scalar lambda, Scalar mu)
    : lambda_(lambda) , mu_(mu)
    {}

    //! Return the first lame parameter
    Scalar lambda() const
    { return lambda_; }

    //! Return the second lame parameter
    Scalar mu() const
    { return mu_; }

    //! set the first lame parameter
    void setLambda(Scalar lambda)
    { lambda_ = lambda; }

    //! set the second lame parameter
    void setMu(Scalar mu)
    { mu_ = mu; }

private:
    Scalar lambda_;
    Scalar mu_;
};
} // end namespace Dumux
#endif
