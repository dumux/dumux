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
 * \ingroup EOS
 * \brief Base class for Peng-Robinson parameters of a
 *        single-component fluid or a mixture
 *
 * See:
 *
 * R. Reid, et al. (1987, pp. 43-44) \cite reid1987
 */
#ifndef DUMUX_PENG_ROBINSON_PARAMS_HH
#define DUMUX_PENG_ROBINSON_PARAMS_HH

namespace Dumux {

/*!
 * \ingroup EOS
 * \brief Stores and provides access to the Peng-Robinson parameters
 *
 * See:
 *
 * R. Reid, et al. (1987, pp. 43-44) \cite reid1987
 */
template <class Scalar>
class PengRobinsonParams
{
public:
    /*!
     * \brief Returns the attractive parameter 'a' of the
     *        Peng-Robinson fluid.
     */
    Scalar a() const
    { return a_; }


    /*!
     * \brief Returns the repulsive parameter 'b' of the Peng-Robinson
     *        fluid.
     */
    Scalar b() const
    { return b_; }

    [[deprecated("Will be removed after Dumux 3.3")]]
    void checkDefined() const
    {}

    /*!
     * \brief Set the attractive parameter 'a' of the Peng-Robinson
     *        fluid.
     * \param value value of the attractive parameter
     */
    void setA(Scalar value)
    { a_ = value; }

    /*!
     * \brief Set the repulsive parameter 'b' of the Peng-Robinson
     *        fluid.
     * \param value value of the repulsive parameter
     */
    void setB(Scalar value)
    { b_ = value; }

protected:
    Scalar a_;
    Scalar b_;
};

} // end namespace

#endif
