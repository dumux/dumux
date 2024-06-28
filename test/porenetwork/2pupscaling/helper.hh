// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_PNM_UPSCALING_HELPER_HH
#define DUMUX_PNM_UPSCALING_HELPER_HH
/*!
 * \file
 *
 * \brief Test for the pore network model
 */
 #include <config.h>

 #include <ctime>
 #include <iostream>

 #include <dune/common/parallel/mpihelper.hh>
 #include <dune/common/timer.hh>
 #include <dune/grid/io/file/dgfparser/dgfexception.hh>
 #include <dune/grid/io/file/vtk.hh>
 #include <dune/grid/io/file/vtk/vtksequencewriter.hh>

 #include <dumux/common/initialize.hh>
 #include <dumux/common/properties.hh>
 #include <dumux/common/parameters.hh>
 #include <dumux/common/dumuxmessage.hh>

#include <dumux/common/properties/model.hh>
#include <dumux/common/properties/grid.hh>
#include <dumux/discretization/porenetwork/gridgeometry.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/thresholdcapillarypressures.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/pore/2p/localrulesforplatonicbody.hh>
#include <dumux/io/grid/porenetwork/gridmanager.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility1p.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility2p.hh>
#include <dumux/porenetwork/common/throatproperties.hh>

#include <dumux/porenetwork/2p/static/staticdrainge.hh>
#include <dumux/io/gnuplotinterface.hh>

#include <dumux/porenetwork/1p/model.hh>// for `TTag::PNMOneP`

// The class that contains a collection of single-phase flow throat transmissibilities
// among them the transmisibility model to be used can be specified in AdvectionType class
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility1p.hh>

// The class that provides specializations for both creeping and non-creeping advection types.
#include <dumux/flux/porenetwork/advection.hh>

// The local residual for incompressible flow is included.
// The one-phase flow model (included above) uses a default implementation of the
// local residual for single-phase flow. However, in this example we are using an
// incompressible fluid phase. Therefore, we are including the specialized local
// residual which contains functionality to analytically compute the entries of
// the Jacobian matrix. We will use this in the main file.
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

// We will use a single liquid phase consisting of a component with constant fluid properties.
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

namespace Dumux::PoreNetwork{

template<class Scalar>
class SimpleFluxVariablesCache
{
    struct WettingLayerCache
    {
        using CreviceResistanceFactor = WettingLayerTransmissibility::CreviceResistanceFactorZhou;
        WettingLayerCache(const SimpleFluxVariablesCache& fluxVariablesCache)
        :fluxVariablesCache_(fluxVariablesCache)
        {}

        Scalar creviceResistanceFactor(const int cornerIdx) const
        { return CreviceResistanceFactor::beta(fluxVariablesCache_.cornerHalfAngle_, fluxVariablesCache_.contactAngle_); }

    private:
        const SimpleFluxVariablesCache<Scalar>& fluxVariablesCache_{};
    };

    using NumCornerVector = Dune::ReservedVector<Scalar, 4>;

public:
    template< class GridGeometry, class Element>
    void update(const GridGeometry& gridGeometry, const Element& element, std::array<Scalar,2> pc)
    {   const auto eIdx = gridGeometry.elementMapper().index(element);
        throatCrossSectionShape_ = gridGeometry.throatCrossSectionShape(/*eIdx*/0);
        throatShapeFactor_ = gridGeometry.throatShapeFactor(eIdx);
        throatLength_ = gridGeometry.throatLength(eIdx);
        throatInscribedRadius_ = gridGeometry.throatInscribedRadius(eIdx);
        Scalar totalThroatCrossSectionalArea = gridGeometry.throatCrossSectionalArea(eIdx);

        const auto numCorners = Throat::numCorners(throatCrossSectionShape_);
        cornerHalfAngle_ = Throat::cornerHalfAngles<Scalar>(throatCrossSectionShape_)[0];
        contactAngle_ = getParam<Scalar>("Problem.ContactAngle");
        surfaceTension_ = getParam<Scalar>("Problem.SurfaceTension");
        pc_ = *std::max_element(pc.begin(), pc.end());
        for (std::size_t i = 0U; i<numCorners; ++i)
            wettingLayerArea_[i] = Throat::wettingLayerCrossSectionalArea(curvatureRadius(), contactAngle_, cornerHalfAngle_);

        throatCrossSectionalArea_[wPhaseIdx()] = std::min(
            std::accumulate(wettingLayerArea_.begin(), wettingLayerArea_.end(), 0.0),
            totalThroatCrossSectionalArea
        );

        throatCrossSectionalArea_[nPhaseIdx()] = totalThroatCrossSectionalArea - throatCrossSectionalArea_[wPhaseIdx()];
    }


    Scalar throatLength() const
    { return throatLength_; }

    Scalar surfaceTension() const
    { return surfaceTension_; }

    Scalar curvatureRadius() const
    { return surfaceTension_ / pc_;}

    Scalar throatInscribedRadius() const
    { return throatInscribedRadius_; }

    Scalar throatShapeFactor() const
    { return throatShapeFactor_; }

    Scalar pc() const
    { return pc_; }

    std::size_t wPhaseIdx() const
    { return 1 - nPhaseIdx_; }

    std::size_t nPhaseIdx() const
    { return nPhaseIdx_; }

    Scalar wettingLayerCrossSectionalArea( int cornerIdx) const
    { return wettingLayerArea_[cornerIdx]; }

    Scalar throatCrossSectionalArea(const int phaseIdx) const
    { return throatCrossSectionalArea_[phaseIdx]; }

    Scalar throatCrossSectionalArea() const
    { return throatCrossSectionalArea_[0] + throatCrossSectionalArea_[1]; }

    const auto& wettingLayerFlowVariables() const
    { return wettingLayerCache_; }

    /*!
     * \brief Returns the throats's cross-sectional shape.
     */
    Throat::Shape throatCrossSectionShape() const
    { return throatCrossSectionShape_; }

private:
    Scalar throatLength_{};
    Scalar surfaceTension_{};
    Scalar throatInscribedRadius_{};
    Scalar throatShapeFactor_{};
    Throat::Shape throatCrossSectionShape_;
    Scalar pc_{};
    Scalar wettingLayerCrossSectionalArea_{};
    Scalar cornerHalfAngle_{};
    Scalar contactAngle_{};
    std::size_t nPhaseIdx_ = 1;
    std::array<Scalar, 2> throatCrossSectionalArea_{};
    NumCornerVector wettingLayerArea_;
    WettingLayerCache wettingLayerCache_ = WettingLayerCache(*this);

};

template<class Scalar>
class SinglePhaseTransmissibility
{
public:

    using SinglePhaseCache = EmptyCache;

    template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables, class FluxVariablesCache>
    static Scalar singlePhaseTransmissibility(const Problem& problem,
                                              const Element& element,
                                              const FVElementGeometry& fvGeometry,
                                              const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                              const ElementVolumeVariables& elemVolVars,
                                              const FluxVariablesCache& fluxVarsCache,
                                              const int phaseIdx)
    {
        using FluidSystem = typename ElementVolumeVariables::VolumeVariables::FluidSystem;
        constexpr bool isGas = FluidSystem::isGas();
        const auto eIdx = problem.gridGeometry().gridView().indexSet().index(element);
        return throatTransmissibility_[eIdx][isGas];
    }

    static void importTransmissibility(const std::vector<std::array<Scalar, 2>> throatTransmissibility)
    {
        throatTransmissibility_ = throatTransmissibility;
    }

private:

    inline static std::vector<std::array<Scalar, 2>> throatTransmissibility_;

};

}
#endif
