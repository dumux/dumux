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
 * \brief Specification of the material params for the heat pipe's capillary
 *        pressure model.
 */
#ifndef DUMUX_HEATPIPELAWPARAMS_HH
#define DUMUX_HEATPIPELAWPARAMS_HH

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Reference implementation of a params for the heat pipe's
 *        material law
 */
// "Use new material laws! Removal after 3.3")
template<class ScalarT>
class HeatPipeLawParams
{
public:
    using Scalar = ScalarT;

    HeatPipeLawParams()
    {}

    HeatPipeLawParams(Scalar p0, Scalar gamma)
    {
        setP0(p0);
        setGamma(gamma);
    }

    /*!
     * \brief Return the \f$\mathrm{[\gamma]}\f$ shape parameter.
     */
    Scalar gamma() const
    { return gamma_; }

    /*!
     * \brief Set the \f$\\mathrm{[gamma]}\f$ shape parameter.
     */
    void setGamma(Scalar v)
    { gamma_ = v; }

    /*!
     * \brief Return the entry pressure \f$\mathrm{[p_0]}\f$ in \f$\mathrm{[Pa]}\f$.
     */
    Scalar p0() const
    { return p0_; }

    /*!
     * \brief Return the entry pressure \f$\mathrm{[p_0]}\f$ in \f$\mathrm{[Pa]}\f$.
     */
    void setP0(Scalar v)
    { p0_ = v; }

private:
    Scalar gamma_;
    Scalar p0_;
};
} // namespace Dumux

#endif
