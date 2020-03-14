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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::StaggeredFreeFlowVelocityOutput
 */
#ifndef DUMUX_STAGGERED_FF_VELOCITYOUTPUT_HH
#define DUMUX_STAGGERED_FF_VELOCITYOUTPUT_HH

#include <dumux/io/velocityoutput.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Velocity output for staggered free-flow models
 */
template<class GridVariables, class SolutionVector>
class StaggeredFreeFlowVelocityOutput : public VelocityOutput<GridVariables>
{
    using ParentType = VelocityOutput<GridVariables>;
    using GridGeometry = typename GridVariables::GridGeometry;
    using Scalar = typename GridVariables::Scalar;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using ElementFluxVarsCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using VolumeVariables = typename GridVariables::VolumeVariables;
    using FluidSystem = typename VolumeVariables::FluidSystem;
    using GridView = typename GridGeometry::GridView;
    // TODO: should be possible to get this somehow
    using Problem = typename std::decay_t<decltype(std::declval<GridVolumeVariables>().problem())>;
    using Element = typename GridView::template Codim<0>::Entity;
    using CoordScalar = typename GridView::ctype;

public:
    using VelocityVector = typename ParentType::VelocityVector;

    /*!
     * \brief Constructor initializes the static data with the initial solution.
     *
     * \param gridVariables The gridVariables
     * \param sol The solution vector
     */
    StaggeredFreeFlowVelocityOutput(const GridVariables& gridVariables, const SolutionVector& sol)
    : gridVariables_(gridVariables)
    , sol_(sol)
    {
        // check if velocity vectors shall be written to the VTK file
        // enable per default
        enableOutput_ = getParamFromGroup<bool>(gridVariables.curGridVolVars().problem().paramGroup(), "Vtk.AddVelocity", true);
    }

    //! Returns whether to enable the velocity output or not
    bool enableOutput() const override { return enableOutput_; }

    //! returns the phase name of a given phase index
    std::string phaseName(int phaseIdx) const override { return FluidSystem::phaseName(phaseIdx); }

    //! returns the number of phases
    int numFluidPhases() const override { return VolumeVariables::numFluidPhases(); }

    //! Calculate the velocities for the scvs in the element
    //! We assume the local containers to be bound to the complete stencil
    void calculateVelocity(VelocityVector& velocity,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVarsCache& elemFluxVarsCache,
                           int phaseIdx) const override
    {
        auto elemFaceVars = localView(gridVariables_.curGridFaceVars());
        elemFaceVars.bindElement(element, fvGeometry, sol_);
        for (auto&& scv : scvs(fvGeometry))
        {
            auto dofIdxGlobal = scv.dofIndex();

            for (auto&& scvf : scvfs(fvGeometry))
            {
                auto dirIdx = scvf.directionIndex();
                velocity[dofIdxGlobal][dirIdx] += 0.5*elemFaceVars[scvf].velocitySelf();
            }
        }
    }

private:
    const GridVariables& gridVariables_;
    const SolutionVector& sol_;

    bool enableOutput_ = true;
};

} // end namespace Dumux

#endif
