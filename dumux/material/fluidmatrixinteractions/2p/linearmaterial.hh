// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Fluidmatrixinteractions
 * \brief   Linear capillary pressure and
 *          relative permeability <-> saturation relations
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_LINEAR_MATERIAL_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_LINEAR_MATERIAL_HH

#include <algorithm>
#include <dumux/common/parameters.hh>
#include <dumux/material/fluidmatrixinteractions/2p/materiallaw.hh>
#include <dumux/material/fluidmatrixinteractions/2p/noregularization.hh>

namespace Dumux::FluidMatrix {

/*!
 * \ingroup Fluidmatrixinteractions
 *
 * \brief Linear capillary pressure and
 * relative permeability <-> saturation relations
 *
 * The entry pressure is reached at \f$\mathrm{ \overline{S}_w = 1}\f$, the maximum
 * capillary pressure is observed at \f$\mathrm{ \overline{S}_w = 0}\f$.
 *
 */
class LinearMaterial
{
public:
    /*!
     * \brief The parameter type
     * \tparam Scalar The scalar type
     */
    template<class Scalar>
    struct Params
    {
        Params(Scalar pcEntry, Scalar pcMax)
        : pcEntry_(pcEntry), pcMax_(pcMax)
        {}

        Scalar pcEntry() const { return pcEntry_; }
        void setPcEntry(Scalar pcEntry) { pcEntry_ = pcEntry; }

        Scalar pcMax() const { return pcMax_; }
        void setPcMax(Scalar pcMax) { pcMax_ = pcMax; }

        bool operator== (const Params& p) const
        {
            return Dune::FloatCmp::eq(pcEntry(), p.pcEntry(), 1e-6)
                   && Dune::FloatCmp::eq(pcMax(), p.pcMax(), 1e-6);
        }

    private:
        Scalar pcEntry_, pcMax_;
    };

    /*!
     * \brief Construct from a subgroup from the global parameter tree
     * \note This will give you nice error messages if a mandatory parameter is missing
     */
    template<class Scalar = double>
    static Params<Scalar> makeParams(const std::string& paramGroup)
    {
        const auto pcEntry = getParamFromGroup<Scalar>(paramGroup, "LinearPcEntry");
        const auto pcMax = getParamFromGroup<Scalar>(paramGroup, "LinearPcMax");
        return {pcEntry, pcMax};
    }

    /*!
     * \brief The capillary pressure-saturation curve
     */
    template<class Scalar>
    static Scalar pc(Scalar swe, const Params<Scalar>& params)
    {
        return (1.0 - swe)*(params.pcMax() - params.pcEntry()) + params.pcEntry();
    }

    /*!
     * \brief The inverse saturation-capillary pressure curve
     */
    template<class Scalar>
    static Scalar swe(Scalar pc, const Params<Scalar>& params)
    {
        return 1.0 - (pc - params.pcEntry())/(params.pcMax() - params.pcEntry());
    }

    /*!
     * \brief The capillary pressure at Swe = 1.0 also called end point capillary pressure
     */
    template<class Scalar>
    static Scalar endPointPc(const Params<Scalar>& params)
    { return params.pcEntry(); }

    /*!
     * \brief The partial derivative of the capillary pressure w.r.t. the effective saturation
     */
    template<class Scalar>
    static Scalar dpc_dswe(Scalar swe, const Params<Scalar>& params)
    {
        return params.pcEntry() - params.pcMax();
    }

    /*!
     * \brief The partial derivative of the effective saturation w.r.t. the capillary pressure
     */
    template<class Scalar>
    static Scalar dswe_dpc(Scalar pc, const Params<Scalar>& params)
    {
        return 1.0/(params.pcEntry() - params.pcMax());
    }

    /*!
     * \brief The relative permeability for the wetting phase
     */
    template<class Scalar>
    static Scalar krw(Scalar swe, const Params<Scalar>& params)
    {
        using std::clamp;
        return clamp(swe, 0.0, 1.0);
    }

    /*!
     * \brief The derivative of the relative permeability
     */
    template<class Scalar>
    static Scalar dkrw_dswe(Scalar swe, const Params<Scalar>& params)
    {
        return 1.0;
    }

    /*!
     * \brief The relative permeability for the non-wetting phase
     */
    template<class Scalar>
    static Scalar krn(Scalar swe, const Params<Scalar>& params)
    {
        using std::clamp;
        return clamp(1.0-swe, 0.0, 1.0);
    }

    /*!
     * \brief The derivative of the relative permeability
     */
    template<class Scalar>
    static Scalar dkrn_dswe(Scalar swe, const Params<Scalar>& params)
    {
        return -1.0;
    }
};

template<typename Scalar = double>
using LinearMaterialDefault = TwoPMaterialLaw<Scalar, LinearMaterial, NoRegularization, TwoPEffToAbsDefaultPolicy>;

} // end namespace Dumux::FluidMatrix

#endif
