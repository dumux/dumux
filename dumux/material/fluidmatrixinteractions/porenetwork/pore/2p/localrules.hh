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
 * \brief Implementation of the capillary pressure and
 * relative permeability <-> saturation relations according to Joekar-Niasar et al., 2010.
 *
 */
#ifndef DUMUX_PNM_2P_LOCAL_RULES_HH
#define DUMUX_PNM_2P_LOCAL_RULES_HH

#include <dumux/porenetworkflow/common/poreproperties.hh>
#include "baselocalrules.hh"
#include "localrules/localrulesforcube.hh"
#include "localrules/localrulesfortetrahedron.hh"
#include "localrules/localrulesforoctahedron.hh"
#include "localrules/localrulesforicosahedron.hh"
#include "localrules/localrulesfordodecahedron.hh"


namespace Dumux
{

template<class ScalarT>
struct LocalRulesTraits
{
    using Cube = TwoPLocalRulesCubeJoekarNiasar<ScalarT>;
    using Tetrahedron = TwoPLocalRulesTetrahedron<ScalarT>;
    using Octahedron = TwoPLocalRulesOctahedron<ScalarT>;
    using Icosahedron = TwoPLocalRulesIcosahedron<ScalarT>;
    using Dodecahedron = TwoPLocalRulesDodecahedron<ScalarT>;
};

template<class ScalarT, class LocalRulesT = LocalRulesTraits<ScalarT>>
struct TwoPLocalRules : public TwoPLocalRulesBase
{
    using Scalar = ScalarT;
    using LocalRules = LocalRulesT;
    using Params = TwoPLocalRulesBase::Params<Scalar>;

    static constexpr bool supportsMultipleGeometries()
    { return true; }

    /*!
     * \brief The capillary pressure-saturation curve
     *
     * \param sw Saturation of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     */
    static Scalar pc(const Params& params, const Scalar sw)
    {
        switch (params.shape)
        {
            case Pore::Shape::cube:
                return LocalRules::Cube::pc(params, sw);
            case Pore::Shape::tetrahedron:
                return LocalRules::Tetrahedron::pc(params, sw);
            case Pore::Shape::octahedron:
                return LocalRules::Octahedron::pc(params, sw);
            case Pore::Shape::icosahedron:
                return LocalRules::Icosahedron::pc(params, sw);
            case Pore::Shape::dodecahedron:
                return LocalRules::Dodecahedron::pc(params, sw);
            default:
                DUNE_THROW(Dune::NotImplemented, "Invalid shape");
        }
    }

    /*!
     * \brief The wetting-phase saturation of a pore body
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen,
     *                  and then the params container is constructed accordingly. Afterwards the values are set there, too.
     */
    static Scalar sw(const Params& params, const Scalar pc)
    {
        switch (params.shape)
        {
            case Pore::Shape::cube:
                return LocalRules::Cube::sw(params, pc);
            case Pore::Shape::tetrahedron:
                return LocalRules::Tetrahedron::sw(params, pc);
            case Pore::Shape::octahedron:
                return LocalRules::Octahedron::sw(params, pc);
            case Pore::Shape::icosahedron:
                return LocalRules::Icosahedron::sw(params, pc);
            case Pore::Shape::dodecahedron:
                return LocalRules::Dodecahedron::sw(params, pc);
            default:
                DUNE_THROW(Dune::NotImplemented, "Invalid shape");
        }
    }

    /*!
     * \brief The partial derivative of the capillary
     *        pressure w.r.t. the wetting phase saturation.
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     * \param sw Saturation of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$
     */
    static Scalar dpc_dsw(const Params& params, const Scalar sw)
    {
        switch (params.shape)
        {
            case Pore::Shape::cube:
                return LocalRules::Cube::dpc_dsw(params, sw);
            case Pore::Shape::tetrahedron:
                return LocalRules::Tetrahedron::dpc_dsw(params, sw);
            case Pore::Shape::octahedron:
                return LocalRules::Octahedron::dpc_dsw(params, sw);
            case Pore::Shape::icosahedron:
                return LocalRules::Icosahedron::dpc_dsw(params, sw);
            case Pore::Shape::dodecahedron:
                return LocalRules::Dodecahedron::dpc_dsw(params, sw);
            default:
                DUNE_THROW(Dune::NotImplemented, "Invalid shape");
        }
    }

    /*!
     * \brief DOCU
     *
     *
     * \param sw Saturation of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     */
    static Scalar dsw_dpc(const Params& params, const Scalar pc)
    {
        switch (params.shape)
        {
            case Pore::Shape::cube:
                return LocalRules::Cube::dsw_dpc(params, pc);
            case Pore::Shape::tetrahedron:
                return LocalRules::Tetrahedron::dsw_dpc(params, pc);
            case Pore::Shape::octahedron:
                return LocalRules::Octahedron::dsw_dpc(params, pc);
            case Pore::Shape::icosahedron:
                return LocalRules::Icosahedron::dsw_dpc(params, pc);
            case Pore::Shape::dodecahedron:
                return LocalRules::Dodecahedron::dsw_dpc(params, pc);
            default:
                DUNE_THROW(Dune::NotImplemented, "Invalid shape");
        }
    }

};

}

#endif // DUMUX_PNM_LOCAL_RULES_HH
