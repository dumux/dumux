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
#include "localrulesforcube.hh"

namespace Dumux
{

template<class ScalarT, class LocalRulesForCube = TwoPLocalRulesCubeJoekarNiasarDefault<ScalarT>>
class MultiShapeTwoPLocalRules
{
    using Scalar = ScalarT;
    using Params = TwoPLocalRulesBase::Params<Scalar>;

    static constexpr bool supportsMultipleGeometries()
    { return true; }

    /*!
     * \brief Return the number of fluid phases
     */
    static constexpr int numFluidPhases()
    { return 2; }

    /*!
     * \brief The capillary pressure-saturation curve
     */
    template<bool enableRegularization = isRegularized()>
    Scalar pc(const Scalar sw) const
    {
        switch (params.shape)
        {
            case Pore::Shape::cube:
                return RegularizedLocalRulesForCube::pc(params, sw);
            default:
                DUNE_THROW(Dune::NotImplemented, "Invalid shape");
        }
    }




private:

    LocalRulesForCube localRulesForCube_;



};

}

#endif // DUMUX_PNM_LOCAL_RULES_HH
