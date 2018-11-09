// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief Class for the evaluation of the porosity subject to precipitation.
 */
#ifndef DUMUX_EFFECTIVE_SOLID_HEATCAPACITY_HH
#define DUMUX_EFFECTIVE_SOLID_HEATCAPACITY_HH

#include <dune/common/deprecated.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/evalsolution.hh>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Calculates the effective solid heat capacity
 */
template<class TypeTag>
class DUNE_DEPRECATED_MSG("Implement SolidSystems instead!") EffectiveSolidHeatCapacity
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    static const int numComponents = GetPropType<TypeTag, Properties::ModelTraits>::numFluidComponents();
    static const int numSolidPhases = GetPropType<TypeTag, Properties::ModelTraits>::numSolidPhases();

    using Element = typename GridView::template Codim<0>:: Entity;

public:
    void init(const SpatialParams& spatialParams)
    {
        spatialParamsPtr_ = &spatialParams;
    }

    // calculates the effective solid heat capacity of multiple solid phases accordin to
    // their volume fractions
    template<class ElementSolution>
    Scalar effectiveSolidHeatCapacity(const Element& element,
                                      const SubControlVolume& scv,
                                      const ElementSolution& elemSol) const
    {
        auto priVars = evalSolution(element, element.geometry(), elemSol, scv.center());

        Scalar sumPrecipitates = 0.0;
        Scalar effCp = 0.0;
        for (unsigned int solidPhaseIdx = 0; solidPhaseIdx < numSolidPhases; ++solidPhaseIdx)
        {
            sumPrecipitates += priVars[numComponents + solidPhaseIdx];
            effCp += priVars[numComponents + solidPhaseIdx]*spatialParams_().solidPhaseHeatCapacity(element, scv, elemSol, solidPhaseIdx);
        }

        return effCp/sumPrecipitates;
    }

private:
    const SpatialParams& spatialParams_() const
    { return *spatialParamsPtr_; }

    const SpatialParams* spatialParamsPtr_;
};

} // namespace Dumux

#endif
