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
 * \brief   Implementation of the capillary pressure and
 *          relative permeability <-> saturation relations with a simple linear law
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_LINEAR_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_LINEAR_HH

#include <cmath>
#include <algorithm>

#include <dumux/common/parameters.hh>

namespace Dumux {
namespace FluidMatrix {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the capillary pressure <->
 *        saturation relation, and relative permeability.
 */
class Linear
{

public:
    /*!
     * \brief The parameter type
     * \tparam Scalar The scalar type
     * \note The van Genuchten laws are parameterized with three parameters: \f$\mathrm{n, m, \alpha}\f$.
     *       For the respective formulas check out the description of the free function.
     */
    template<class Scalar>
    struct Params
    {
        Scalar entryPc, maxPc;
    };

    /*!
     * \brief Construct from a subgroup from the global parameter tree
     * \note This will give you nice error messages if a mandatory parameter is missing
     */
    template<class Scalar>
    static Params<Scalar> makeParams(const std::string& paramGroup)
    {
        const auto entryPc = getParam<Scalar>(paramGroup + ".EntryPc");
        const auto maxPc = getParam<Scalar>(paramGroup + ".MaxPc");
        return {entryPc, maxPc};
    }

    /*!
     * \brief The capillary pressure-saturation curve according to van Genuchten.
     * \param swe Effective saturation of the wetting phase \f$\mathrm{\overline{S}_w}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     */
    template<class Scalar>
    static Scalar pc(Scalar swe, const Params<Scalar>& params)
    {
        return (1.0 - swe)*(params.maxPc - params.entryPc) + params.entryPc;
    }

    /*!
     * \brief The saturation-capillary pressure curve according to van Genuchten.
     * \param pc Capillary pressure \f$\mathrm{p_C}\f$ in \f$\mathrm{[Pa]}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     */
    template<class Scalar>
    static Scalar sw(Scalar pc, const Params<Scalar>& params)
    {
        return 1.0 - (pc - params.entryPc)/(params.maxPc - params.entryPc);
    }

    /*!
     * \brief The capillary pressure at Swe = 1.0 also called end point capillary pressure
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     */
    template<class Scalar>
    static Scalar endPointPc(const Params<Scalar>& params)
    { return params.entryPc; }

    /*!
     * \brief The partial derivative of the capillary
     *        pressure w.r.t. the effective saturation
     * \param swe Effective saturation of the wetting phase \f$\mathrm{\overline{S}_w}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     */
    template<class Scalar>
    static Scalar dpc_dswe(Scalar swe, const Params<Scalar>& params)
    {
        return params.entryPc-params.maxPc;
    }

    /*!
     * \brief The partial derivative of the effective
     *        saturation to the capillary pressure
     * \param pc Capillary pressure \f$\mathrm{p_C}\f$ in \f$\mathrm{[Pa]}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     */
    template<class Scalar>
    static Scalar dswe_dpc(Scalar pc, const Params<Scalar>& params)
    {
        return 1.0/(params.entryPc-params.maxPc);
    }

    /*!
     * \brief The relative permeability for the wetting phase of the medium
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     */
    template<class Scalar>
    static Scalar krw(Scalar swe, const Params<Scalar>& params)
    {
        using std::max;
        using std::min;
        return max(min(swe, 1.0), 0.0);
    }

    /*!
     * \brief The derivative of the relative permeability for the
     *        wetting phase in regard to the wetting saturation of the medium
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     */
    template<class Scalar>
    static Scalar dkrw_dswe(Scalar swe, const Params<Scalar>& params)
    {
        return 1.0;
    }

    /*!
     * \brief The relative permeability for the non-wetting phase of the medium
     * \param swe The mobile saturation of the wetting phase.
     */
    template<class Scalar>
    static Scalar krn(Scalar swe, const Params<Scalar>& params)
    {
        using std::max;
        using std::min;
        return max(min(1.0-swe, 1.0), 0.0);
    }

    /*!
     * \brief The derivative of the relative permeability for the non-wetting phase in regard to the wetting saturation
     * \param swe The mobile saturation of the wetting phase.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     */
    template<class Scalar>
    static Scalar dkrn_dswe(Scalar swe, const Params<Scalar>& params)
    {
        return -1.0;
    }
};

} // end namespace FluidMatrix
} // end namespace Dumux

#endif
