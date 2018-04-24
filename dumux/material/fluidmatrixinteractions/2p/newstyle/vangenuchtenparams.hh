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
 *
 * \brief Specification of the material parameters
 *       for the van Genuchten constitutive relations.
 */
#ifndef VAN_GENUCHTEN_PARAMS_HH
#define VAN_GENUCHTEN_PARAMS_HH

namespace Dumux
{

/*!
 *
 * \brief Specification of the material parameters
 *       for the van Genuchten constitutive relations.
 *
 *       In this implementation setting either the \f$\mathrm{n}\f$ or \f$\mathrm{m}\f$ shape parameter
 *       automatically calculates the other. I.e. they cannot be set independently.
 *
 * \ingroup fluidmatrixinteractionsparams
 */
template<class Scalar, class VanGenuchtenType>
class VanGenuchtenParams : public VanGenuchtenType::RegularizationPolicy::Params,
                           public VanGenuchtenType::EffToAbsPolicy::Params
{
public:
    using MaterialLaw = VanGenuchtenType;
    using RegularizationPolicy = typename VanGenuchtenType::RegularizationPolicy;
    using EffToAbsPolicy = typename VanGenuchtenType::EffToAbsPolicy;

    //! Constructor
    //! Initializes some parameters with default values
    VanGenuchtenParams()
    : RegularizationPolicy::Params(),
      EffToAbsPolicy::Params()
    {
        //! Initialize / precompute parameters for the regularization
        initRegularizationParams(*this);
    }

    /*!
     * \brief Return the \f$\mathrm{\alpha}\f$ shape parameter \f$\mathrm{[1/Pa]}\f$ of van Genuchten's
     *        curve.
     */
    Scalar vgAlpha() const
    { return vgAlpha_; }

    /*!
     * \brief Set the \f$\mathrm{\alpha}\f$ shape parameter \f$\mathrm{[1/Pa]}\f$ of van Genuchten's
     *        curve.
     */
    void setVgAlpha(const Scalar v)
    { vgAlpha_ = v; }

    /*!
     * \brief Return the \f$\mathrm{m}\f$ shape parameter \f$\mathrm{[-]}\f$ of van Genuchten's
     *        curve.
     */
    Scalar vgm() const
    { return vgm_; }

    /*!
     * \brief Set the \f$\mathrm{m}\f$ shape parameter \f$\mathrm{[-]}\f$ of van Genuchten's
     *        curve.
     *
     * The \f$\mathrm{n}\f$ shape parameter is set to \f$\mathrm{n = \frac{1}{1 - m}}\f$
     */
    void setVgm(const Scalar m)
    { vgm_ = m; vgn_ = 1/(1 - vgm_); }

    /*!
     * \brief Return the \f$\mathrm{n}\f$ shape parameter \f$\mathrm{[-]}\f$ of van Genuchten's
     *        curve.
     */
    Scalar vgn() const
    { return vgn_; }

    /*!
     * \brief Set the \f$\mathrm{n}\f$ shape parameter \f$\mathrm{[-]}\f$ of van Genuchten's
     *        curve.
     *
     * The \f$\mathrm{n}\f$ shape parameter is set to \f$\mathrm{m = 1 - \frac{1}{n}}\f$
     */
    void setVgn(const Scalar n)
    { vgn_ = n; vgm_ = 1 - 1/vgn_; }

private:
    Scalar vgAlpha_;
    Scalar vgm_;
    Scalar vgn_;
};

} // end namespace Dumux

#endif
