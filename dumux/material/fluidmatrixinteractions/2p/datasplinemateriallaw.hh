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
 * \brief Pc- and Kr-sw curves based on monotone splines through given data points
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_DATA_SPLINE_MATERIAL_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_DATA_SPLINE_MATERIAL_HH

#include <algorithm>
#include <dumux/common/parameters.hh>
#include <dumux/common/monotonecubicspline.hh>
#include <dumux/material/fluidmatrixinteractions/fluidmatrixinteraction.hh>
#include <dumux/material/fluidmatrixinteractions/2p/materiallaw.hh>

namespace Dumux::FluidMatrix {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Pc- and Kr-sw curves based on monotone splines through given data points
 * \tparam S the type for scalar numbers
 * \tparam approximatePcSwInverse if this is set true, the
 *          spline approximates sw(pc) and evaluating pc(sw) needs spline inversion.
 *          if this is false, the spline approximates pc(sw) and evaluating
 *          sw(pc) needs spline inversion. Spline inversion is rather expensive
 *          since it has to be done numerically.
 */
template <class S, bool approximatePcSwInverse = false>
class DataSplineTwoPMaterialLaw
: public Adapter<DataSplineTwoPMaterialLaw<S, approximatePcSwInverse>, PcKrSw>
{
public:

    using Scalar = S;

    /*!
     * \brief Return the number of fluid phases
     */
    static constexpr int numFluidPhases()
    { return 2; }

    static constexpr bool isRegularized()
    { return false; }

    /*!
     * \brief Deleted default constructor (so we are never in an undefined state)
     * \note store owning pointers to laws instead if you need default-constructible objects
     */
    DataSplineTwoPMaterialLaw() = delete;

    /*!
     * \brief Construct from a subgroup from the global parameter tree
     */
    explicit DataSplineTwoPMaterialLaw(const std::string& paramGroup)
    {
        using V = std::vector<Scalar>;
        const std::string swPcGroup = paramGroup.empty() ? "Pc" : paramGroup + ".Pc";
        const auto swPc = getParamFromGroup<V>(swPcGroup, "SwData");
        const std::string krwPcGroup = paramGroup.empty() ? "Krw" : paramGroup + ".Krw";
        const auto swKrw = getParamFromGroup<V>(krwPcGroup, "SwData");
        const std::string krnPcGroup = paramGroup.empty() ? "Krn" : paramGroup + ".Krn";
        const auto swKrn = getParamFromGroup<V>(krnPcGroup, "SwData");

        const auto pc = getParamFromGroup<V>(paramGroup, "PcData");
        const auto krw = getParamFromGroup<V>(paramGroup, "KrwData");
        const auto krn = getParamFromGroup<V>(paramGroup, "KrnData");

        updateData_(swPc, pc, swKrw, krw, swKrn, krn);
    }

    /*!
     * \brief Construct from user data
     */
    DataSplineTwoPMaterialLaw(const std::vector<Scalar>& swPc,
                              const std::vector<Scalar>& pc,
                              const std::vector<Scalar>& swKrw,
                              const std::vector<Scalar>& krw,
                              const std::vector<Scalar>& swKrn,
                              const std::vector<Scalar>& krn)
    {
        updateData_(swPc, pc, swKrw, krw, swKrn, krn);
    }

    /*!
     * \brief The capillary pressure-saturation curve
     */
    Scalar pc(const Scalar sw) const
    {
        if constexpr (approximatePcSwInverse)
            return pcSpline_.evalInverse(sw);
        else
            return pcSpline_.eval(sw);
    }

    /*!
     * \brief The partial derivative of the capillary pressure w.r.t. the saturation
     */
    Scalar dpc_dsw(const Scalar sw) const
    {
        if constexpr (approximatePcSwInverse)
            return 1.0/pcSpline_.evalDerivative(pcSpline_.evalInverse(sw));
        else
            return pcSpline_.evalDerivative(sw);
    }

    /*!
     * \brief The saturation-capillary pressure curve
     */
    Scalar sw(const Scalar pc) const
    {
        if constexpr (approximatePcSwInverse)
            return pcSpline_.eval(pc);
        else
            return pcSpline_.evalInverse(pc);
    }

    /*!
     * \brief The partial derivative of the saturation to the capillary pressure
     */
    Scalar dsw_dpc(const Scalar pc) const
    {
        if constexpr (approximatePcSwInverse)
            return pcSpline_.evalDerivative(pc);
        else
            return 1.0/pcSpline_.evalDerivative(pcSpline_.evalInverse(pc));
    }

    /*!
     * \brief The relative permeability for the wetting phase
     */
    Scalar krw(const Scalar sw) const
    {
        return krwSpline_.eval(sw);
    }

    /*!
     * \brief The derivative of the relative permeability for the wetting phase w.r.t. saturation
     */
    Scalar dkrw_dsw(const Scalar sw) const
    {
        return krwSpline_.evalDerivative(sw);
    }

    /*!
     * \brief The relative permeability for the non-wetting phase
     */
    Scalar krn(const Scalar sw) const
    {
        return krnSpline_.eval(sw);
    }

    /*!
     * \brief The derivative of the relative permeability for the non-wetting phase w.r.t. saturation
     */
    Scalar dkrn_dsw(const Scalar sw) const
    {
        return krnSpline_.evalDerivative(sw);
    }

private:
    void updateData_(const std::vector<Scalar>& swPc,
                     const std::vector<Scalar>& pc,
                     const std::vector<Scalar>& swKrw,
                     const std::vector<Scalar>& krw,
                     const std::vector<Scalar>& swKrn,
                     const std::vector<Scalar>& krn)
    {
        if constexpr (approximatePcSwInverse)
            pcSpline_.updatePoints(pc, swPc);
        else
            pcSpline_.updatePoints(swPc, pc);

        krwSpline_.updatePoints(swKrw, krw);
        krnSpline_.updatePoints(swKrn, krn);
    }

    MonotoneCubicSpline<Scalar> pcSpline_;
    MonotoneCubicSpline<Scalar> krwSpline_;
    MonotoneCubicSpline<Scalar> krnSpline_;
};

} // end namespace Dumux::FluidMatrix

#endif
