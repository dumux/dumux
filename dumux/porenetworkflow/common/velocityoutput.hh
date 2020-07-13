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
 *
 * \brief Velocity output for pore-network models
 */
#ifndef DUMUX_PNM_VELOCITYOUTPUT_HH
#define DUMUX_PNM_VELOCITYOUTPUT_HH

#include <string>
#include <dumux/io/velocityoutput.hh>

namespace Dumux {

/*!
 * \brief Velocity output for pore-network models
 */
template<class GridVariables, class FluxVariables>
class PNMVelocityOutput : public VelocityOutput<GridVariables>
{
    using ParentType = VelocityOutput<GridVariables>;
    using Scalar = typename GridVariables::Scalar;
    using GridGeometry = typename GridVariables::GridGeometry;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVariables::VolumeVariables;
    using FluidSystem = typename VolumeVariables::FluidSystem;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;

public:
    using VelocityVector = typename ParentType::VelocityVector;

    //! Returns the phase name of a given phase index.
    std::string phaseName(int phaseIdx) const override { return FluidSystem::phaseName(phaseIdx); }

    //! Returns the number of phases.
    int numFluidPhases() const override { return VolumeVariables::numFluidPhases(); }

    /*!
     * \brief Constructor initializes the static data with the initial solution.
     *
     * \param problem The problem to be solved
     */
    PNMVelocityOutput(const GridVariables& gridVariables)
    : gridVariables_(gridVariables)
    {
        velocityOutput_ = getParamFromGroup<bool>(problem_().paramGroup(), "Vtk.AddVelocity");
    }

    bool enableOutput() const override
    { return velocityOutput_; }

    //! Calculate the velocities for the scvs in the element
    //! We assume the local containers to be bound to the complete stencil
    void calculateVelocity(VelocityVector& velocity,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVariablesCache& elemFluxVarsCache,
                           int phaseIdx) const override
    {
        if (!velocityOutput_)
            return;

        const auto geometry = element.geometry();

        auto tmpVelocity = (geometry.corner(1) - geometry.corner(0));
        tmpVelocity /= tmpVelocity.two_norm();

        const auto eIdxGlobal = fvGeometry.gridGeometry().elementMapper().index(element);
        velocity[eIdxGlobal] = 0.0;

        for (auto&& scvf : scvfs(fvGeometry))
        {
            if (scvf.boundary())
                continue;

            // get the volume flux divided by the area of the
            // subcontrolvolume face in the reference element
            // TODO: Divide by extrusion factor!!?
            const Scalar flux = getFlux_(element, fvGeometry, scvf, elemVolVars, elemFluxVarsCache, phaseIdx);

            tmpVelocity *= flux;
            velocity[eIdxGlobal] = tmpVelocity;
        }
    }

private:

    Scalar getFlux_(const Element& element,
                    const FVElementGeometry& fvGeometry,
                    const SubControlVolumeFace& scvf,
                    const ElementVolumeVariables& elemVolVars,
                    const ElementFluxVariablesCache& elemFluxVarsCache,
                    const int phaseIdx) const
    {
        const Scalar localArea = elemFluxVarsCache[scvf].throatCrossSectionalArea(phaseIdx);

        // make sure the phase is actually present (2p)
        if (localArea > 0.0)
        {
            // instantiate the flux variables
            FluxVariables fluxVars;
            fluxVars.init(problem_(), element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

            // the upwind term to be used for the volume flux evaluation
            auto upwindTerm = [phaseIdx](const auto& volVars) { return volVars.mobility(phaseIdx); };
            return fluxVars.advectiveFlux(phaseIdx, upwindTerm) / localArea;
        }
        else
            return 0.0;
    }

    const auto& problem_() const { return gridVariables_.curGridVolVars().problem(); }

    bool velocityOutput_;

    const GridVariables& gridVariables_;
};

} // end namespace Dumux

#endif
