// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup HyperelasticVolIso
 * \brief Volume variables for the volumetric-isochoric split hyperelastic mixed u–p model.
 */
#ifndef DUMUX_SOLIDMECHANICS_HYPERELASTIC_VOLISO_VOLUMEVARIABLES_HH
#define DUMUX_SOLIDMECHANICS_HYPERELASTIC_VOLISO_VOLUMEVARIABLES_HH

#include <dumux/common/volumevariables.hh>

namespace Dumux {

/*!
 * \ingroup HyperelasticVolIso
 * \brief Volume variables for the momentum subdomain (displacement, PQ1Bubble or PQ2).
 *
 * Provides `displacement(i)` for direction `i`.  Supports the new-style
 * `update(elemSol, problem, fvGeometry, ipData)` used by
 * `Experimental::CVFE::CVFEGridVariablesCache` and its hybrid variant.
 * For non-SCV DOFs (PQ2 edge midpoints) the index is remapped via a proxy.
 */
template<class Traits>
class HyperelasticVolIsoMomentumVolumeVariables : public BasicVolumeVariables<Traits>
{
    using Scalar = typename Traits::PrimaryVariables::value_type;
public:
    using PrimaryVariables = typename Traits::PrimaryVariables;
    using Indices = typename Traits::ModelTraits::Indices;

    Scalar displacement(int i) const { return this->priVar(Indices::displacement(i)); }

    template<class ES, class P, class FVG, class IpD>
    void update(const ES& elemSol, const P& problem, const FVG& fvGeometry, const IpD& ipData)
    {
        const auto idx = ipData.localDofIndex();
        if (idx < fvGeometry.numScv())
        {
            BasicVolumeVariables<Traits>::update(elemSol, problem, fvGeometry.element(), fvGeometry.scv(idx));
        }
        else
        {
            // PQ2 edge-midpoint DOF: remap to scv(0) with index substituted
            const auto& scv0 = fvGeometry.scv(0);
            struct Proxy {
                const ES& s; std::size_t from, to;
                auto operator[](std::size_t i) const { return (i == from) ? s[to] : s[i]; }
            } proxy{elemSol, scv0.localDofIndex(), idx};
            BasicVolumeVariables<Traits>::update(proxy, problem, fvGeometry.element(), scv0);
        }
    }
};

/*!
 * \ingroup HyperelasticVolIso
 * \brief Volume variables for the pressure subdomain (bulk pressure, Box/P1).
 *
 * Provides `pressure()`.  All DOFs are vertex SCVs so no index remapping is needed.
 */
template<class Traits>
class HyperelasticVolIsoPressureVolumeVariables : public BasicVolumeVariables<Traits>
{
    using Scalar = typename Traits::PrimaryVariables::value_type;
public:
    using PrimaryVariables = typename Traits::PrimaryVariables;

    Scalar pressure() const { return this->priVar(0); }

    template<class ES, class P, class FVG, class IpD>
    void update(const ES& elemSol, const P& problem, const FVG& fvGeometry, const IpD& ipData)
    {
        BasicVolumeVariables<Traits>::update(elemSol, problem, fvGeometry.element(),
                                             fvGeometry.scv(ipData.localDofIndex()));
    }
};

} // end namespace Dumux
#endif
