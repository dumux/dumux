// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
