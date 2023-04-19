// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of a NAPL adsorption model
 */
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
