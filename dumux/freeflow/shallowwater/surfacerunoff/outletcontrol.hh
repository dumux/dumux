// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup SurfaceRunoff
 * \brief Stage-discharge relations for where a channel network leaves the domain.
 */
#ifndef DUMUX_FREEFLOW_SHALLOWWATER_SURFACERUNOFF_OUTLETCONTROL_HH
#define DUMUX_FREEFLOW_SHALLOWWATER_SURFACERUNOFF_OUTLETCONTROL_HH

#include <cmath>
#include <string>

#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>

namespace Dumux::SurfaceRunoff {

/*!
 * \ingroup SurfaceRunoff
 * \brief What sets the discharge where a channel network leaves the modelled domain.
 *
 * The long-wave approximations neglect inertia and cannot represent a critical-flow control,
 * so the stage-discharge relation of a weir, sill or gate has to be imposed as a boundary
 * condition.
 */
enum class OutletControl
{
    closed, //!< nothing leaves; the domain is a sealed basin
    freeOutfall, //!< water leaves at the rate the local bed slope conveys it
    weir, //!< critical flow over a crest, discharge set by the head on it
    prescribedDischarge, //!< a gate schedule or a measured release
    fixedStage //!< the level is held, e.g. by a reservoir standing downstream
};

inline OutletControl outletControlFromName(const std::string& name)
{
    if (name == "closed") return OutletControl::closed;
    if (name == "freeoutfall") return OutletControl::freeOutfall;
    if (name == "weir") return OutletControl::weir;
    if (name == "prescribeddischarge") return OutletControl::prescribedDischarge;
    if (name == "fixedstage") return OutletControl::fixedStage;
    DUNE_THROW(ParameterException, "Unknown outlet control '" << name << "', expected "
               << "'closed', 'freeoutfall', 'weir', 'prescribeddischarge' or 'fixedstage'");
}

inline std::string outletControlName(const OutletControl control)
{
    switch (control)
    {
        case OutletControl::closed: return "closed";
        case OutletControl::freeOutfall: return "freeoutfall";
        case OutletControl::weir: return "weir";
        case OutletControl::prescribedDischarge: return "prescribeddischarge";
        case OutletControl::fixedStage: return "fixedstage";
    }
    DUNE_THROW(Dune::InvalidStateException, "Unknown outlet control");
}

/*!
 * \ingroup SurfaceRunoff
 * \brief Free discharge over a weir crest, per unit length of crest.
 *
 * `q = C * head^(3/2)` with the head of the water level above the crest, from critical flow
 * over the crest. `C` is about 1.7 m^(1/2)/s for a broad-crested weir and 1.8-2.0 m^(1/2)/s for
 * a sharp-crested one.
 *
 * Zero below the crest and C1 across it. Valid for free (modular) flow only, i.e. with the
 * tailwater below the crest.
 */
template<class Scalar>
Scalar weirDischargePerLength(const Scalar waterLevel, const Scalar crestLevel,
                              const Scalar dischargeCoefficient)
{
    using std::max, std::sqrt;
    const auto head = max(Scalar(0.0), waterLevel - crestLevel);
    return dischargeCoefficient*head*sqrt(head);
}

/*!
 * \ingroup SurfaceRunoff
 * \brief Discharge over a weir that the tailwater can drown, per unit length of crest.
 *
 * The free-flow discharge times Villemonte's (1947) submergence factor
 * `(1 - (hd/hu)^{3/2})^{0.385}`, with the heads `hu` and `hd` of the upstream and downstream
 * levels above the crest. Free flow is recovered whenever the tailwater is below the crest.
 * The factor has an infinite slope at equal levels, so above the head ratio `maxSubmergence`
 * the discharge decreases linearly to zero at equal levels instead.
 */
template<class Scalar>
Scalar weirDischargePerLength(const Scalar upstreamLevel, const Scalar downstreamLevel,
                              const Scalar crestLevel, const Scalar dischargeCoefficient,
                              const Scalar maxSubmergence = 0.999)
{
    using std::max, std::min, std::sqrt, std::pow;
    const auto hu = max(Scalar(0.0), upstreamLevel - crestLevel);
    if (hu <= 0.0)
        return 0.0;

    const auto free = dischargeCoefficient*hu*sqrt(hu);
    const auto hd = max(Scalar(0.0), downstreamLevel - crestLevel);
    if (hd <= 0.0)
        return free;

    const auto ratio = hd/hu;
    if (ratio >= 1.0)
        return 0.0;

    const auto capped = min(ratio, maxSubmergence);
    const auto factor = pow(max(Scalar(0.0), 1.0 - capped*sqrt(capped)), Scalar(0.385));
    if (ratio <= maxSubmergence)
        return free*factor;

    return free*factor*(1.0 - ratio)/(1.0 - maxSubmergence);
}

/*!
 * \ingroup SurfaceRunoff
 * \brief The head a weir needs in order to pass a given discharge per unit length.
 *
 * The inverse of the free-flow `weirDischargePerLength`.
 */
template<class Scalar>
Scalar weirHeadForDischargePerLength(const Scalar dischargePerLength,
                                     const Scalar dischargeCoefficient)
{
    using std::max, std::cbrt;
    const auto ratio = max(Scalar(0.0), dischargePerLength)/dischargeCoefficient;
    return cbrt(ratio*ratio);
}

/*!
 * \ingroup SurfaceRunoff
 * \brief Settings of a single outlet, read from the parameter tree.
 */
template<class Scalar>
struct OutletSettings
{
    OutletControl control = OutletControl::closed;
    Scalar crestLevel = 0.0; //!< weir only, absolute elevation
    Scalar dischargeCoefficient = 1.7; //!< weir only, broad-crested by default
    Scalar discharge = 0.0; //!< prescribedDischarge only, m^3/s
    Scalar stage = 0.0; //!< fixedStage only, absolute elevation

    static OutletSettings fromParams(const std::string& paramGroup)
    {
        OutletSettings settings;
        settings.control = outletControlFromName(
            getParamFromGroup<std::string>(paramGroup, "Problem.OutletControl", "closed")
        );

        if (settings.control == OutletControl::weir)
        {
            settings.crestLevel = getParamFromGroup<Scalar>(paramGroup, "Problem.OutletCrestLevel");
            settings.dischargeCoefficient = getParamFromGroup<Scalar>(
                paramGroup, "Problem.OutletDischargeCoefficient", 1.7
            );
        }
        else if (settings.control == OutletControl::prescribedDischarge)
            settings.discharge = getParamFromGroup<Scalar>(paramGroup, "Problem.OutletDischarge");
        else if (settings.control == OutletControl::fixedStage)
            settings.stage = getParamFromGroup<Scalar>(paramGroup, "Problem.OutletStage");

        return settings;
    }
};

} // end namespace Dumux::SurfaceRunoff

#endif
