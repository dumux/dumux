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
 * \brief Implementation of a NAPL adsorption model.
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_INTERACTIONS_3P_NAPL_ADSORPTION
#define DUMUX_MATERIAL_FLUIDMATRIX_INTERACTIONS_3P_NAPL_ADSORPTION

#include <algorithm>
#include <dune/common/float_cmp.hh>
#include <dumux/common/parameters.hh>
#include <dumux/material/fluidmatrixinteractions/fluidmatrixinteraction.hh>

namespace Dumux::FluidMatrix {

template<class Scalar>
class ThreePNAPLAdsorption : Adapter<ThreePNAPLAdsorption<Scalar>, Adsorption>
{
public:

    struct Params
    {
        Params(const Scalar rhoBulk, const Scalar kdNAPL)
        : rhoBulk_(rhoBulk), kdNAPL_(kdNAPL) {}

        Scalar rhoBulk() const { return rhoBulk_; }
        void setRhoBulk(const Scalar rhoBulk) { rhoBulk_ = rhoBulk; }

        Scalar kdNAPL() const { return kdNAPL_; }
        void setKdNAPL(const Scalar kdNAPL) { kdNAPL_ = kdNAPL; }

        bool operator== (const Params& p) const
        {
            return Dune::FloatCmp::eq(kdNAPL(), p.kdNAPL(), 1e-6)
                   && Dune::FloatCmp::eq(rhoBulk(), p.rhoBulk(), 1e-6);
        }

    private:
        Scalar rhoBulk_, kdNAPL_;
    };

    /*!
     * \brief Construct from a subgroup from the global parameter tree
     * \note This will give you nice error messages if a mandatory parameter is missing
     */
    static Params makeParams(const std::string& paramGroup)
    {
        const auto rhoBulk = getParamFromGroup<Scalar>(paramGroup, "ThreePNAPLAdsorptionRhoBulk");
        const auto kdNAPL = getParamFromGroup<Scalar>(paramGroup, "ThreePNAPLAdsorptionKdNAPL");
        return {rhoBulk, kdNAPL};
    }

    ThreePNAPLAdsorption(const Params& params): params_(params) {}

    ThreePNAPLAdsorption(const std::string& paramGroup)
    : params_(makeParams(paramGroup)) {}

   /*!
    * \brief the basis for calculating adsorbed NAPL in storage term
    * \param params Array of parameters
    */
   Scalar bulkDensTimesAdsorpCoeff () const
   {
      return params_.rhoBulk() * params_.kdNAPL();
   }

private:
    Params params_;
};
} // end namespace Dumux::FluidMatrix

#endif
