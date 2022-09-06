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
 * \ingroup Assembly
 * \brief Volume variables deflection helper.
 */

#ifndef DUMUX_ASSEMBLY_VOLVARS_DEFLECTION_HELPER_HH
#define DUMUX_ASSEMBLY_VOLVARS_DEFLECTION_HELPER_HH

#include <type_traits>
#include <dune/common/reservedvector.hh>

#ifndef DOXYGEN

namespace Dumux::Detail {

template<class VolVarAccessor, class FVElementGeometry>
class VolVarsDeflectionHelper
{
    static constexpr int maxNumscvs = FVElementGeometry::maxNumElementScvs;

    using SubControlVolume = typename FVElementGeometry::GridGeometry::SubControlVolume;
    using VolumeVariablesRef = std::invoke_result_t<VolVarAccessor, const SubControlVolume&>;
    using VolumeVariables = std::decay_t<VolumeVariablesRef>;
    static_assert(std::is_lvalue_reference_v<VolumeVariablesRef>
                  && !std::is_const_v<std::remove_reference_t<VolumeVariablesRef>>);

public:
    VolVarsDeflectionHelper(VolVarAccessor&& accessor,
                            const FVElementGeometry& fvGeometry,
                            bool deflectAllVolVars)
    : deflectAll_(deflectAllVolVars)
    , accessor_(std::move(accessor))
    , fvGeometry_(fvGeometry)
    {
        if (deflectAll_)
            for (const auto& curScv : scvs(fvGeometry))
                origVolVars_.push_back(accessor_(curScv));
    }

    void setCurrent(const SubControlVolume& scv)
    {
        if (!deflectAll_)
        {
            origVolVars_.clear();
            origVolVars_.push_back(accessor_(scv));
        }
    }

    template<class ElementSolution, class Problem>
    void deflect(const ElementSolution& elemSol,
                 const SubControlVolume& scv,
                 const Problem& problem)
    {
        if (deflectAll_)
            for (const auto& curScv : scvs(fvGeometry_))
                accessor_(curScv).update(elemSol, problem, fvGeometry_.element(), curScv);
        else
            accessor_(scv).update(elemSol, problem, fvGeometry_.element(), scv);
    }

    void restore(const SubControlVolume& scv)
    {
        if (!deflectAll_)
            accessor_(scv) = origVolVars_[0];
        else
            for (const auto& curScv : scvs(fvGeometry_))
                accessor_(curScv) = origVolVars_[curScv.localDofIndex()];
    }

private:
    const bool deflectAll_;
    VolVarAccessor accessor_;
    const FVElementGeometry& fvGeometry_;
    Dune::ReservedVector<VolumeVariables, maxNumscvs> origVolVars_;
};

template<class Accessor, class FVElementGeometry>
VolVarsDeflectionHelper(Accessor&&, FVElementGeometry&&) -> VolVarsDeflectionHelper<Accessor, std::decay_t<FVElementGeometry>>;

} // end namespace Dumux::Detail

#endif // DOXYGEN

#endif
