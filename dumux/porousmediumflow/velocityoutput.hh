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
 * \ingroup PorousmediumflowModels
 * \brief Velocity output for porous media models.
 */

#ifndef DUMUX_POROUSMEDIUMFLOW_VELOCITYOUTPUT_HH
#define DUMUX_POROUSMEDIUMFLOW_VELOCITYOUTPUT_HH

#include <memory>
#include <dune/common/float_cmp.hh>

#include <dumux/common/parameters.hh>
#include <dumux/io/velocityoutput.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/porousmediumflow/velocity.hh>

namespace Dumux {

/*!
 * \ingroup PorousmediumflowModels
 * \brief Velocity output policy for implicit (porous media) models.
 */
template<class GridVariables, class FluxVariables>
class PorousMediumFlowVelocityOutput : public VelocityOutput<GridVariables>
{
    using ParentType = VelocityOutput<GridVariables>;
    using GridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementFluxVarsCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using VolumeVariables = typename GridVariables::VolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using FluidSystem = typename VolumeVariables::FluidSystem;
    using Scalar = typename GridVariables::Scalar;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr bool isBox = GridGeometry::discMethod == DiscretizationMethod::box;
    static constexpr int dofCodim = isBox ? dim : 0;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using Problem = typename GridVolumeVariables::Problem;
    using VelocityBackend = PorousMediumFlowVelocity<GridVariables, FluxVariables>;

public:
    using VelocityVector = typename ParentType::VelocityVector;

    /*!
     * \brief Constructor initializes the static data with the initial solution.
     *
     * \param gridVariables The grid variables
     */
    PorousMediumFlowVelocityOutput(const GridVariables& gridVariables)
    {
        // check, if velocity output can be used (works only for cubes so far)
        enableOutput_ = getParamFromGroup<bool>(gridVariables.curGridVolVars().problem().paramGroup(), "Vtk.AddVelocity");
        if (enableOutput_)
            velocityBackend = std::make_unique<VelocityBackend>(gridVariables);
    }

    //! Returns whether or not velocity output is enabled.
    bool enableOutput() const override { return enableOutput_; }

    //! Returns the phase name of a given phase index.
    std::string phaseName(int phaseIdx) const override { return FluidSystem::phaseName(phaseIdx); }

    //! Returns the number of phases.
    int numFluidPhases() const override { return VolumeVariables::numFluidPhases(); }

    //! Calculates the velocities for the scvs in the element.
    //! We assume the local containers to be bound to the complete stencil.
    void calculateVelocity(VelocityVector& velocity,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVarsCache& elemFluxVarsCache,
                           int phaseIdx) const override
    {
        if (enableOutput_)
            velocityBackend->calculateVelocity(velocity, element, fvGeometry, elemVolVars, elemFluxVarsCache, phaseIdx);
    }

private:
    bool enableOutput_;
    std::unique_ptr<VelocityBackend> velocityBackend;
};

} // end namespace Dumux

#endif
