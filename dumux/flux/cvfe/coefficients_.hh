// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEFlux
 * \brief Coefficients of the flux laws of control-volume finite element schemes at interpolation points
 */
#ifndef DUMUX_FLUX_CVFE_COEFFICIENTS__HH
#define DUMUX_FLUX_CVFE_COEFFICIENTS__HH

#include <concepts>
#include <type_traits>
#include <utility>

#include <dumux/common/concepts/ipdata_.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/flux/facetensoraverage.hh>

namespace Dumux::Experimental::CVFE {

/*!
 * \ingroup CVFEFlux
 * \brief A coefficient at the points of a sub-control volume face, averaged over the two
 *        sub-control volumes sharing the face (see Dumux::faceTensorAverage)
 *
 * \param coefficient the coefficient per local dof, `coefficient(vars, args...)`
 * \return the coefficient per interpolation point, `(elemDisc, elemVars, ipData, args...)`,
 *         defined for points that know their sub-control volume face
 */
template<class Coefficient>
auto faceTensorAverageAtScvf(Coefficient&& coefficient)
{
    return [c = std::forward<Coefficient>(coefficient)]
           (const auto& elemDisc, const auto& elemVars, const Concept::ScvfIpData auto& ipData, const auto&... args)
    {
        const auto& scvf = elemDisc.scvf(ipData.scvfIndex());
        const auto& insideVars = elemVars[elemDisc.scv(scvf.insideScvIdx())];
        const auto& outsideVars = elemVars[elemDisc.scv(scvf.outsideScvIdx())];
        return Dumux::faceTensorAverage(c(insideVars, args...), c(outsideVars, args...), ipData.unitOuterNormal());
    };
}

} // end namespace Dumux::Experimental::CVFE

namespace Dumux::Detail::CVFE {

//! A coefficient given per local dof, `coefficient(vars, args...)`
template<class Coefficient, class ElementVariables, class... Args>
concept LocalDofCoefficient = std::invocable<const Coefficient&, const typename ElementVariables::Variables&, const Args&...>;

/*!
 * \brief The value of a coefficient at an interpolation point
 * \note A coefficient per local dof is interpolated with the shape values, any other
 *       coefficient is evaluated at the point, `coefficient(elemDisc, elemVars, ipData, args...)`.
 */
template<class ElementDiscretization, class ElementVariables, class IpData, class IpCache, class Coefficient, class... Args>
auto coefficientAtIp(const ElementDiscretization& elemDisc,
                     const ElementVariables& elemVars,
                     const IpData& ipData,
                     const IpCache& ipCache,
                     const Coefficient& coefficient,
                     const Args&... args)
{
    if constexpr (LocalDofCoefficient<Coefficient, ElementVariables, Args...>)
    {
        const auto& shapeValues = ipCache.shapeValues();
        std::decay_t<std::invoke_result_t<const Coefficient&, const typename ElementVariables::Variables&, const Args&...>> value(0.0);
        for (const auto& localDof : localDofs(elemDisc))
            value += shapeValues[localDof.index()][0]*coefficient(elemVars[localDof], args...);
        return value;
    }
    else
        return coefficient(elemDisc, elemVars, ipData, args...);
}

} // end namespace Dumux::Detail::CVFE

#endif
