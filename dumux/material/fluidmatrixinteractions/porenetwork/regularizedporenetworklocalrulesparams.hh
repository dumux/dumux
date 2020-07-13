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
 *
 * \brief   Parameters that are necessary for the \em regularization of
 *          the PNMLocalRules.
 */
#ifndef DUMUX_REGULARIZED_PNM_LOCAL_RULES_PARAMS_HH
#define DUMUX_REGULARIZED_PNM_LOCAL_RULES_PARAMS_HH

#include "porenetworklocalrulesparams.hh"

namespace Dumux
{
/*!
 * \brief   Parameters that are necessary for the \em regularization of
 *          the RegularizedPNMLocalRules.
 *
 *        \ingroup fluidmatrixinteractionsparams
 */
template <class Scalar>
class RegularizedPNMLocalRulesParams : public Dumux::PNMLocalRulesParams<Scalar>
{
    using ParentType = Dumux::PNMLocalRulesParams<Scalar>;

public:

    // use the un-regularized class's constructors
    using ParentType::ParentType;

    /*!
     * \brief Threshold saturation below which the capillary pressure
     *        is regularized.
     *
     * This is just 1%. If you need a different value, overload this
     * class.
     */
    Scalar lowSw() const
    {
        // Most problems are very sensitive to this value
        // (e.g. making it smaller might result in negative
        // pressures)
        //
        // If you want to use a different regularization threshold,
        // overload this class and supply the new class as second
        // template parameter for the RegularizedPNMLocalRules!
        static const Scalar value = getParam<Scalar>("Regularization.LowSw", 1e-2);
        return value;
    }

    /*!
     * \brief Threshold saturation below which the capillary pressure
     *        is regularized.
     *
     * This is just 95%. If you need a different value, overload this
     * class.
     */
    Scalar highSw() const
    {
        // Most problems are very sensitive to this value
        // (e.g. making it smaller might result in negative
        // pressures)
        //
        // If you want to use a different regularization threshold,
        // overload this class and supply the new class as second
        // template parameter for the RegularizedPNMLocalRules!
        static const Scalar value = getParam<Scalar>("Regularization.HighSw", 0.95);
        return value;
    }

    /*!
     * \brief Slope for the pc-S curve for Sw > 1.0
     *
     * We use a very high negative to trigger very low values for pn, thus
     * avoiding negative non-wetting phase saturations.
     */
    Scalar slopeHighSw() const
    {
        static const Scalar slope = getParam<Scalar>("Regularization.SlopeHighSw", -1e9);
        return slope;
    }


};
} // namespace Dumux

#endif
