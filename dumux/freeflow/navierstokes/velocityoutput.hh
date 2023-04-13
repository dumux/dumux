// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::StaggeredFreeFlowVelocityOutput
 */
#ifndef DUMUX_FREEFLOW_NAVIERSTOKES_VELOCITYOUTPUT_HH
#define DUMUX_FREEFLOW_NAVIERSTOKES_VELOCITYOUTPUT_HH

#include <type_traits>
#include <dune/common/exceptions.hh>
#include <dumux/io/velocityoutput.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/method.hh>
#include <dumux/freeflow/navierstokes/momentum/velocityreconstruction.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Velocity output for staggered free-flow models
 */
template<class GridVariables>
class NavierStokesVelocityOutput : public VelocityOutput<GridVariables>
{
    using ParentType = VelocityOutput<GridVariables>;
    using GridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using ElementFluxVarsCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using VolumeVariables = typename GridVariables::VolumeVariables;
    using FluidSystem = typename VolumeVariables::FluidSystem;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FieldType = typename ParentType::FieldType;

public:
    using VelocityVector = typename ParentType::VelocityVector;

    NavierStokesVelocityOutput(const std::string& paramGroup = "")
    {
        enableOutput_ = getParamFromGroup<bool>(paramGroup, "Vtk.AddVelocity", true);
    }

    //! Returns whether to enable the velocity output or not
    bool enableOutput() const override { return enableOutput_; }

    //! returns the phase name of a given phase index
    std::string phaseName(int phaseIdx) const override { return FluidSystem::phaseName(phaseIdx); }

    //! returns the number of phases
    int numFluidPhases() const override { return VolumeVariables::numFluidPhases(); }

    //! returns the field type
    FieldType fieldType() const override { return FieldType::element; }

    //! Calculate the velocities for the scvs in the element
    //! We assume the local containers to be bound to the complete stencil
    void calculateVelocity(VelocityVector& velocity,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVarsCache& elemFluxVarsCache,
                           int phaseIdx) const override
    {
        using CouplingManager = std::decay_t<decltype(elemVolVars.gridVolVars().problem().couplingManager())>;
        using MomGG = std::decay_t<decltype(std::declval<CouplingManager>().problem(CouplingManager::freeFlowMomentumIndex).gridGeometry())>;
        if constexpr (MomGG::discMethod == DiscretizationMethods::fcstaggered)
            calculateVelocityForStaggeredGrid_(velocity, element, fvGeometry, elemVolVars);
        else if constexpr (DiscretizationMethods::isCVFE<typename MomGG::DiscretizationMethod>)
            calculateVelocityForCVFESchemes_(velocity, element, fvGeometry, elemVolVars);
        else
            DUNE_THROW(Dune::NotImplemented, "Navier-Stokes velocity output for scheme " << MomGG::discMethod);
    }

private:
    void calculateVelocityForStaggeredGrid_(VelocityVector& velocity,
                                            const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& elemVolVars) const
    {
        const auto eIdx = fvGeometry.gridGeometry().elementMapper().index(element);
        const auto getFaceVelocity = [&](const FVElementGeometry& fvG, const auto& scvf)
        {
            return elemVolVars.gridVolVars().problem().faceVelocity(element, fvGeometry, scvf);
        };

        velocity[eIdx] = StaggeredVelocityReconstruction::cellCenterVelocity(getFaceVelocity, fvGeometry);
    }

    void calculateVelocityForCVFESchemes_(VelocityVector& velocity,
                                          const Element& element,
                                          const FVElementGeometry& fvGeometry,
                                          const ElementVolumeVariables& elemVolVars) const
    {
        const auto eIdx = fvGeometry.gridGeometry().elementMapper().index(element);
        velocity[eIdx] = elemVolVars.gridVolVars().problem().elementVelocity(fvGeometry);
    }


    bool enableOutput_;
};

} // end namespace Dumux

#endif
