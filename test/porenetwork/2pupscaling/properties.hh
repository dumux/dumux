// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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

#include "problem_static.hh"

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
        const auto& shape = gridGeometry.throatCrossSectionShape(/*eIdx*/0);
        throatShapeFactor_ = gridGeometry.throatShapeFactor(eIdx);
        throatLength_ = gridGeometry.throatLength(eIdx);
        throatInscribedRadius_ = gridGeometry.throatInscribedRadius(eIdx);
        Scalar totalThroatCrossSectionalArea = gridGeometry.throatCrossSectionalArea(eIdx);

        const auto numCorners = Throat::numCorners(shape);
        cornerHalfAngle_ = Throat::cornerHalfAngles<Scalar>(shape)[0];
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

private:
    Scalar throatLength_{};
    Scalar surfaceTension_{};
    Scalar throatInscribedRadius_{};
    Scalar throatShapeFactor_{};
    Scalar pc_{};
    Scalar wettingLayerCrossSectionalArea_{};
    Scalar cornerHalfAngle_{};
    Scalar contactAngle_{};
    std::size_t nPhaseIdx_ = 1;
    std::array<Scalar, 2> throatCrossSectionalArea_{};
    NumCornerVector wettingLayerArea_;
    WettingLayerCache wettingLayerCache_ = WettingLayerCache(*this);

};
}

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct PNMTWOPStatic { using InheritsFrom = std::tuple<GridProperties, ModelProperties>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::PNMTWOPStatic> { using type = Dune::FoamGrid<1, 3>; };

template<class TypeTag>
struct GridGeometry<TypeTag, TTag::PNMTWOPStatic>
{
private:
    static constexpr bool enableCache = false;
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = Dumux::PoreNetwork::GridGeometry<Scalar, GridView, enableCache>;
};

// The flux variables cache
template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::PNMTWOPStatic>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = PoreNetwork::SimpleFluxVariablesCache<Scalar>;
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PNMTWOPStatic> { using type = Dumux::PoreNetwork::DrainageProblemStatic<TypeTag>; };

} // end namespace Dumux::Properties
