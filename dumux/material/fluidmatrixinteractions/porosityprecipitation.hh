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
#ifndef DUMUX_POROSITY_PRECIPITATION_HH
#define DUMUX_POROSITY_PRECIPITATION_HH

#include <dumux/discretization/evalsolution.hh>

namespace Dumux
{

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Calculates the porosity depending on the volume fractions of precipitated minerals.
 */
template<class TypeTag>
class PorosityPrecipitation
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SpatialParams = typename GET_PROP_TYPE(TypeTag, SpatialParams);
    using ElementSolution = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    static const int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);
    static const int numSolidPhases = GET_PROP_VALUE(TypeTag, NumSPhases);

    using Element = typename GridView::template Codim<0>:: Entity;

public:
    void init(const SpatialParams& spatialParams)
    {
        spatialParamsPtr_ = &spatialParams;
    }

    /*!
     * \brief calculates the porosity in a sub-control volume
     * \param element element
     * \param elemSol the element solution
     * \param scv sub control volume
     */
    Scalar evaluatePorosity(const Element& element,
                            const SubControlVolume& scv,
                            const ElementSolution& elemSol) const
    {
        auto priVars = evalSolution(element, element.geometry(), elemSol, scv.center());

        Scalar sumPrecipitates = 0.0;
        for (unsigned int solidPhaseIdx = 0; solidPhaseIdx < numSolidPhases; ++solidPhaseIdx)
            sumPrecipitates += priVars[numComponents + solidPhaseIdx];

        auto minPoro = spatialParams_().minPorosity(element, scv);

        using std::max;
        return max(minPoro, spatialParams_().initialPorosity(element, scv) - sumPrecipitates);
    }

private:
    const SpatialParams& spatialParams_() const
    { return *spatialParamsPtr_; }

    const SpatialParams* spatialParamsPtr_;
};

} // namespace Dumux

#endif
