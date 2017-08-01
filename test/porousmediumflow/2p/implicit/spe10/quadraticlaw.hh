// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Markus Wolff                                      *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Implementation of quadratic relative-permeability-saturation relations.
 * @author Markus Wolff
 */
#ifndef DUMUX_QUADRATICLAW_HH
#define DUMUX_QUADRATICLAW_HH

#include <algorithm>

namespace Dumux
{
/*!
 * \ingroup fluidmatrixinteractionslaws
 *
 * \brief Implementation of quadratic relative-permeability-saturation relations.
 *
 * This class provides quadratic relative-permeability-saturation relations:
 *
 *  \f[
 *       k_{r_w} = \overline{S}_w^2
 *  \f]
 *
 *  \f[
 *       k_{r_n} = \overline{S}_n^2 = (1 - \overline{S}_n)^2
 *  \f]
 *
 * A class providing functions for capillary-pressure-saturation relations has to be given as template argument.
 *
 * @tparam ScalarT The type of a scalar
 * @tparam PCLawT The class providing the capillary-pressure-saturation relations
 *
 */
template <class ScalarT, class PCLawT>
class QuadraticLaw
{
public:
    typedef PCLawT PCLaw;
    typedef typename PCLaw::Params Params;
    typedef ScalarT Scalar;

    /*!
     * \brief The capillary pressure-saturation curve.
     *
     * \param Swe       Effective saturation of the wetting phase \f$\overline{S}_w\f$
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Capillary pressure.
     */
    static Scalar pc(const Params &params, Scalar Swe)
    {
        return PCLaw::pc(params, Swe);
    }

    DUNE_DEPRECATED_MSG("use pc() (uncapitalized 'c') instead")
    static Scalar pC(const Params &params, Scalar Swe)
    {
        return PCLaw::pc(params, Swe);
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     *
     *
     * \param pC        Capillary pressure \f$p_C\f$
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Effective wetting phase saturation.
     */
    static Scalar sw(const Params &params, Scalar pC)
    {
        return PCLaw::sw(params, pC);
    }

    DUNE_DEPRECATED_MSG("use sw() (uncapitalized 's') instead")
    static Scalar Sw(const Params &params, Scalar pC)
    {
        return PCLaw::sw(params, pC);
    }

    /*!
     * \brief The partial derivative of the capillary
     *        pressure w.r.t. the effective saturation.
     *
     * \param Swe       Effective saturation of the wetting phase \f$\overline{S}_w\f$
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Partial derivative of \f$p_c\f$ w.r.t. effective saturation.
    */
    static Scalar dpc_dsw(const Params &params, Scalar Swe)
    {
        return PCLaw::dpc_dsw(params, Swe);
    }

    DUNE_DEPRECATED_MSG("use dpc_dsw() (uncapitalized 'c', 's') instead")
    static Scalar dpC_dSw(const Params &params, Scalar Swe)
    {
        return PCLaw::dpc_dsw(params, Swe);
    }

    /*!
     * \brief The partial derivative of the effective
     *        saturation w.r.t. the capillary pressure.
     *
     * \param pC        Capillary pressure \f$p_C\f$
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Partial derivative of effective saturation w.r.t. \f$p_c\f$.
     */
    static Scalar dsw_dpc(const Params &params, Scalar pC)
    {
        return PCLaw::dsw_dpc(params, pC);
    }

    DUNE_DEPRECATED_MSG("use dsw_dpc() (uncapitalized 's', 'c') instead")
    static Scalar dSw_dpC(const Params &params, Scalar pC)
    {
        return PCLaw::dsw_dpc(params, pC);
    }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium.
     *
     * \param Swe       The mobile saturation of the wetting phase.
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Relative permeability of the wetting phase.
     */
    static Scalar krw(const Params &params, Scalar Swe)
    {
        if (Swe <= 0.0)
            return 0.0;
        else if (Swe >= 1.0)
            return 1.0;

        return Swe*Swe;
    };

    /*!
     * \brief The derivative of the relative permeability for the
     *        wetting phase with regard to the wetting saturation of the
     *        medium.
     *
     * \param Swe       The mobile saturation of the wetting phase.
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Derivative of the relative permeability of the wetting phase w.r.t. effective wetting phase saturation calculated as implied by Brooks & Corey.
     */
    static Scalar dkrw_dsw(const Params &params, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);

        return 2*Swe;
    };

    DUNE_DEPRECATED_MSG("use dkrw_dsw() (uncapitalized 's') instead")
    static Scalar dkrw_dSw(const Params &params, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);

        return 2*Swe;
    };

    /*!
     * \brief The relative permeability for the non-wetting phase of
     *        the medium.
     *
     * \param Swe       The mobile saturation of the wetting phase.
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Relative permeability of the non-wetting phase.
     */
    static Scalar krn(const Params &params, Scalar Swe)
    {
        if (Swe >= 1.0)
            return 0.0;
        else if (Swe <= 0.0)
            return 1.0;

        return (1-Swe)*(1-Swe);
    }

    /*!
     * \brief The derivative of the relative permeability for the
     *        non-wetting phase in regard to the wetting saturation of
     *        the medium.
     *
     * \param Swe       The mobile saturation of the wetting phase.
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Derivative of the relative permeability of the non-wetting phase w.r.t. effective wetting phase saturation calculated as implied by Brooks & Corey.
     */
    static Scalar dkrn_dsw(const Params &params, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);

        return -2*(1-Swe);
    };

    DUNE_DEPRECATED_MSG("use dkrn_dsw() (uncapitalized 's') instead")
    static Scalar dkrn_dSw(const Params &params, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);

        return -2*(1-Swe);
    };


};
}

#endif // QUADRATICLAW_HH
