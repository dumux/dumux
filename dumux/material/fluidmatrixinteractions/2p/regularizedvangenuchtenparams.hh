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
 *          VanGenuchten "material law".
 */
#ifndef REGULARIZED_VAN_GENUCHTEN_PARAMS_HH
#define REGULARIZED_VAN_GENUCHTEN_PARAMS_HH

#warning "This header is deprecated. Removal after 3.3. Use new material laws."

#include <dune/common/float_cmp.hh>

#include "vangenuchtenparams.hh"

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief   Parameters that are necessary for the \em regularization of
 *          VanGenuchten "material law".
 */
template<class ScalarT>
class [[deprecated("Use new material laws! Removal after 3.3")]] RegularizedVanGenuchtenParams : public VanGenuchtenParams<ScalarT>
{
public:
    using Scalar = ScalarT;
    using Parent = VanGenuchtenParams<Scalar>;

    RegularizedVanGenuchtenParams()
    {
        initialize();
    }

    RegularizedVanGenuchtenParams(Scalar vgAlpha, Scalar vgN)
        : Parent(vgAlpha, vgN)
    {
        initialize();
    }

    /*!
     * \brief Equality comparison with another set of params
     */
    template<class OtherParams>
    bool operator== (const OtherParams& otherParams) const
    {
        return Dune::FloatCmp::eq(pcLowSw_, otherParams.pcLowSw(), /*eps*/1e-6*pcLowSw_)
               && Dune::FloatCmp::eq(pcHighSw_, otherParams.pcHighSw(), /*eps*/1e-6*pcHighSw_)
               && Dune::FloatCmp::eq(krnLowSw_, otherParams.krnLowSw(), /*eps*/1e-6*krnLowSw_)
               && Dune::FloatCmp::eq(krwHighSw_, otherParams.krwHighSw(), /*eps*/1e-6*krwHighSw_)
               && Parent::operator==(otherParams);
    }

    /*!
     * \brief Sets some default regularization thresholds
     */
    void initialize()
    {
        setPcLowSw(0.01);
        setPcHighSw(0.99);
        setKrnLowSw(0.1);
        setKrwHighSw(0.9);
    }

    /*!
     * \brief Set the threshold saturation below which the capillary pressure is regularized.
     *
     * Most problems are very sensitive to this value (e.g. making it smaller might
     * result in very high capillary pressures)
     */
    void setPcLowSw(Scalar pcLowSw)
    {
        pcLowSw_ = pcLowSw;
    }

    /*!
     * \brief Threshold saturation below which the capillary pressure is regularized.
     */
    Scalar pcLowSw() const
    {
        return pcLowSw_;
    }

    /*!
     * \brief Set the threshold saturation above which the capillary pressure is regularized.
     */
    void setPcHighSw(Scalar pcHighSw)
    {
        pcHighSw_ = pcHighSw;
    }

    /*!
     * \brief Threshold saturation above which the capillary pressure is regularized.
     *
     * Most problems are very sensitive to this value (e.g. making it smaller might
     * result in negative capillary pressures).
     */
    Scalar pcHighSw() const
    {
        return pcHighSw_;
    }

    /*!
     * \brief Set the threshold saturation below which the relative
     *        permeability of the nonwetting phase gets regularized.
     */
    void setKrnLowSw(Scalar krnLowSw)
    {
        krnLowSw_ = krnLowSw;
    }

    /*!
     * \brief Threshold saturation below which the relative
     *        permeability of the nonwetting phase gets regularized.
     */
    Scalar krnLowSw() const
    {
        return krnLowSw_;
    }

    /*!
     * \brief Set the threshold saturation above which the relative
     *        permeability of the wetting phase gets regularized.
     */
    void setKrwHighSw(Scalar krwHighSw)
    {
        krwHighSw_ = krwHighSw;
    }

    /*!
     * \brief Threshold saturation above which the relative
     *        permeability of the wetting phase gets regularized.
     */
    Scalar krwHighSw() const
    {
        return krwHighSw_;
    }

private:
    Scalar pcLowSw_;
    Scalar pcHighSw_;
    Scalar krnLowSw_;
    Scalar krwHighSw_;
};
} // namespace Dumux

#endif
