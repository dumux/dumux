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
 * \brief Base classes for standard pore-local pc-Sw curves.
 */
#ifndef DUMUX_PNM_2P_BASE_LOCAL_RULES_HH
#define DUMUX_PNM_2P_BASE_LOCAL_RULES_HH

#include <dumux/common/parameters.hh>
#include <dumux/porenetworkflow/common/poreproperties.hh>

namespace Dumux
{

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Base class for all standard pore-local pc-Sw curves.
 */
struct TwoPLocalRulesBase
{
    /*!
     * \brief The parameter type used for all standard pore-local pc-Sw curves.
     *
     * \note For sake of compatibility, we need to have one unique set of parameters for different types of curves,
     *       even if this means that some parameters might be unused for certain laws, pore geometries, etc.
     *
     * \tparam Scalar The scalar type
     */
    template<class Scalar>
    struct Params
    {
        Scalar poreRadius, contactAngle, surfaceTension;
        Pore::Shape shape;
    };

    /*!
     * \brief Convenience function to create parameters for standard pore-local pc-Sw curves.
     *
     * \tparam Scalar The scalar type
     */
    template<class Scalar>
    static Params<Scalar> makeParams(const Scalar poreRadius, const Scalar contactAngle,
                                     const Scalar surfaceTension, const Pore::Shape shape)
    {
        return Params<Scalar>{poreRadius, contactAngle, surfaceTension, shape};
    }

    //! This is just for compatibility with the REV-scale models.
    //! Could be removed if the pore-network models' volume variables
    //! do not inherit from the REV-scale volume variables.
    template<class... Args>
    static double krw(Args&&...)
    { return 1.0; }

    //! This is just for compatibility with the REV-scale models.
    //! Could be removed if the pore-network models' volume variables
    //! do not inherit from the REV-scale volume variables.
    template<class... Args>
    static double krn(Args&&...)
    { return 1.0; }

    //! This is just for compatibility with the REV-scale models.
    //! Could be removed if the pore-network models' volume variables
    //! do not inherit from the REV-scale volume variables.
    template<class... Args>
    static double dkrw_dsw(Args&&...)
    { return 0.0; }

    //! This is just for compatibility with the REV-scale models.
    //! Could be removed if the pore-network models' volume variables
    //! do not inherit from the REV-scale volume variables.
    template<class... Args>
    static double dkrn_dsw(Args&&...)
    { return 0.0; }

};

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Base class for all regularized standard pore-local pc-Sw curves.
 */
struct RegularizedTwoPLocalRulesBase : public TwoPLocalRulesBase
{
    /*!
     * \brief The parameter type used for all regularized standard pore-local pc-Sw curves.
     *
     * \note For sake of compatibility, we need to have one unique set of parameters for different types of curves,
     *       even if this means that some parameters might be unused for certain laws, pore geometries, etc.
     *
     * \tparam Scalar The scalar type
     */
    template<class Scalar>
    struct Params : public TwoPLocalRulesBase::Params<Scalar>
    {
        Scalar lowSw, highSw, slopeHighSw;
    };

    /*!
     * \brief Convenience function to create parameters for regularized standard pore-local pc-Sw curves.
     *
     * \tparam Scalar The scalar type
     */
    template<class Scalar>
    static Params<Scalar> makeParams(const Scalar poreRadius, const Scalar contactAngle,
                                     const Scalar surfaceTension, const Pore::Shape shape)
    {
        static const Scalar lowSw = getParam<Scalar>("Regularization.LowSw", 1e-2);
        static const Scalar highSw = getParam<Scalar>("Regularization.HighSw", 0.95);
        static const Scalar slopeHighSw = getParam<Scalar>("Regularization.SlopeHighSw", -1e9);
        return Params<Scalar>{{poreRadius, contactAngle, surfaceTension, shape}, lowSw, highSw, slopeHighSw};
    }
};

}

#endif
