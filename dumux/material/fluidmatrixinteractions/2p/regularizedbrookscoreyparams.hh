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
 *          the Brooks-Corey capillary pressure model.
 */
#ifndef DUMUX_REGULARIZED_BROOKS_COREY_PARAMS_HH
#define DUMUX_REGULARIZED_BROOKS_COREY_PARAMS_HH

#warning "This header is deprecated. Removal after 3.3. Use new material laws."

#include <dune/common/float_cmp.hh>

#include "brookscoreyparams.hh"

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief   Parameters that are necessary for the \em regularization of
 *          the Brooks-Corey capillary pressure model.
 */
template <class ScalarT>
class [[deprecated("Use new material laws! Removal after 3.3")]] RegularizedBrooksCoreyParams : public BrooksCoreyParams<ScalarT>
{
    using BrooksCoreyParams = Dumux::BrooksCoreyParams<ScalarT>;

public:
    using Scalar = ScalarT;

    RegularizedBrooksCoreyParams()
        : BrooksCoreyParams()
    {
        setThresholdSw(0.01);
    }

    RegularizedBrooksCoreyParams(Scalar pe, Scalar lambda)
        : BrooksCoreyParams(pe, lambda)
    {
        setThresholdSw(0.01);
    }

    /*!
     * \brief Equality comparison with another set of params
     */
    template<class OtherParams>
    bool operator== (const OtherParams& otherParams) const
    {
        return Dune::FloatCmp::eq(thresholdSw_, otherParams.thresholdSw(), /*eps*/1e-6*thresholdSw_)
               && BrooksCoreyParams::operator==(otherParams);
    }

    /*!
     * \brief Set the threshold saturation below which the capillary pressure
     *        is regularized.
     *
     * Most problems are very sensitive to this value (e.g. making it smaller
     * might result in negative pressures).
     */
    void setThresholdSw(Scalar thresholdSw)
    {
        thresholdSw_ = thresholdSw;
    }

    /*!
     * \brief Threshold saturation below which the capillary pressure
     *        is regularized.
     */
    Scalar thresholdSw() const
    {
        return thresholdSw_;
    }

private:
    Scalar thresholdSw_;
};
} // namespace Dumux

#endif
