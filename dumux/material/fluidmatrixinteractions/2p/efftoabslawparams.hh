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
 * \brief A default implementation of the parameters for the adapter
 *        class to convert material laws from effective to absolute
 *        saturations.
 */
#ifndef DUMUX_EFF_TO_ABS_LAW_PARAMS_HH
#define DUMUX_EFF_TO_ABS_LAW_PARAMS_HH

#warning "This header is deprecated. Removal after 3.3. Use new material laws."

#include <dune/common/float_cmp.hh>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief A default implementation of the parameters for the adapter
 *        class to convert material laws from effective to absolute
 *        saturations.
 */
template <class EffLawParamsT>
class [[deprecated("Use new material laws! Removal after 3.3")]] EffToAbsLawParams : public EffLawParamsT
{
    using EffLawParams = EffLawParamsT;
public:
    using Scalar = typename EffLawParams::Scalar;

    EffToAbsLawParams()
        : EffLawParams()
    { swr_ = snr_ = 0; }

    /*!
     * \brief Equality comparison with another set of params
     */
    template<class OtherParams>
    bool operator== (const OtherParams& otherParams) const
    {
        return Dune::FloatCmp::eq(swr_, otherParams.swr(), /*eps*/1e-6*swr_)
               && Dune::FloatCmp::eq(snr_, otherParams.snr(), /*eps*/1e-6*snr_)
               && EffLawParamsT::operator==(otherParams);
    }

    /*!
     * \brief Return the residual wetting saturation.
     */
    Scalar swr() const
    { return swr_; }

    /*!
     * \brief Set the residual wetting saturation.
     */
    void setSwr(Scalar v)
    { swr_ = v; }

    /*!
     * \brief Return the residual nonwetting saturation.
     */
    Scalar snr() const
    { return snr_; }

    /*!
     * \brief Set the residual nonwetting saturation.
     */
    void setSnr(Scalar v)
    { snr_ = v; }

private:
    Scalar swr_;
    Scalar snr_;
};

} // end namespace Dumux

#endif
