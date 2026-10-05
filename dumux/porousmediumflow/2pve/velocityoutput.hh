// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later

/*!
 * \file
 * \ingroup TwoPVEModel
 * \brief Column-averaged velocity output for the vertically integrated two-phase model.
 */

#ifndef DUMUX_TWOPVE_VELOCITY_OUTPUT_HH
#define DUMUX_TWOPVE_VELOCITY_OUTPUT_HH

#include <dumux/porousmediumflow/velocityoutput.hh>

namespace Dumux {

/*!
 * \ingroup TwoPVEModel
 * \brief Outputs the column-averaged Darcy velocity in metres per second
 *
 * The generic velocity output removes the volume-variable extrusion factor.
 * For VE this factor includes 1/H to reduce the full-height geometry. Dividing
 * the resulting integrated velocity by H recovers the column-averaged velocity.
 */
template<class GridVariables, class FluxVariables>
class TwoPVEVelocityOutput : public PorousMediumFlowVelocityOutput<GridVariables, FluxVariables>
{
    using ParentType = PorousMediumFlowVelocityOutput<GridVariables, FluxVariables>;
    using GridGeometry = typename GridVariables::GridGeometry;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVarsCache = typename GridVariables::GridFluxVariablesCache::LocalView;

public:
    using ParentType::ParentType;
    using VelocityVector = typename ParentType::VelocityVector;

    void calculateVelocity(VelocityVector& velocity,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVarsCache& elemFluxVarsCache,
                           int phaseIdx) const override
    {
        if (!this->enableOutput())
            return;

        ParentType::calculateVelocity(velocity, element, fvGeometry, elemVolVars, elemFluxVarsCache, phaseIdx);
        constexpr int dim = GridGeometry::GridView::dimension;
        const auto& gridGeometry = fvGeometry.gridGeometry();
        const auto height = gridGeometry.bBoxMax()[dim-1] - gridGeometry.bBoxMin()[dim-1];
        velocity[gridGeometry.elementMapper().index(element)] /= height;
    }
};

} // namespace Dumux

#endif
