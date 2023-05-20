// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief The spatial parameters class for the Hhat problem with multiple solid spheres
 */

#ifndef DUMUX_TEST_MULTIDOMAIN_DUALNETWORK_SPATIALPARAMS_HH
#define DUMUX_TEST_MULTIDOMAIN_DUALNETWORK_SPATIALPARAMS_HH

#include <dumux/porenetwork/solidenergy/spatialparams.hh>
#include <dumux/porenetwork/1p/spatialparams.hh>

namespace Dumux::PoreNetwork {

/*!
 * \brief The spatial parameters class.
 */
template<class GridGeometry, class Scalar>
class FluidSpatialParams
: public PoreNetwork::OnePDefaultSpatialParams<GridGeometry, Scalar>
{
    using ParentType = PoreNetwork::OnePDefaultSpatialParams<GridGeometry, Scalar>;
    using GlobalPosition = typename GridGeometry::GlobalCoordinate;
public:
    // export permeability type
    using PermeabilityType = Scalar;

    template<class GridData>
    FluidSpatialParams(std::shared_ptr<const GridGeometry> fvGridGeometry, const GridData& gridData)
    : ParentType(fvGridGeometry)
    {
        throatCenter_.resize(fvGridGeometry->gridView().size(0));
        poreExtendedRadius_.resize(fvGridGeometry->gridView().size(1));

        for (const auto& element : elements(fvGridGeometry->gridView()))
        {
            const auto eIdx = fvGridGeometry->elementMapper().index(element);
            const auto& params = gridData.parameters(element);
            const auto posX = params[gridData.parameterIndex("ThroatCenterX")];
            const auto posY = params[gridData.parameterIndex("ThroatCenterY")];
            const auto posZ = params[gridData.parameterIndex("ThroatCenterZ")];
            throatCenter_[eIdx] = GlobalPosition{posX, posY, posZ};
        }

        for (const auto& vertex : vertices(fvGridGeometry->gridView()))
        {
            const auto vIdx = fvGridGeometry->vertexMapper().index(vertex);
            poreExtendedRadius_[vIdx] = gridData.getParameter(vertex, "PoreExtendedRadius");
        }
    }

    const GlobalPosition& throatCenter(const std::size_t eIdx) const
    { return throatCenter_[eIdx]; }

    template<class Element, class SubControlVolume, class ElementSolutionVector>
    Scalar poreExtendedRadius(const Element& element,
                              const SubControlVolume& scv,
                              const ElementSolutionVector& elemSol) const
    { return poreExtendedRadius_[scv.dofIndex()]; }

    Scalar poreExtendedRadius(std::size_t dofIdx) const
    { return poreExtendedRadius_[dofIdx]; }

private:
    std::vector<GlobalPosition> throatCenter_;
    std::vector<Scalar> poreExtendedRadius_;
};

/*!
 * \brief The spatial parameters class.
 */
template<class GridGeometry, class Scalar>
class SolidSpatialParams
: public PoreNetwork::SolidEnergySpatialParams<GridGeometry, Scalar>
{
    using ParentType = PoreNetwork::SolidEnergySpatialParams<GridGeometry, Scalar>;
    using GlobalPosition = typename GridGeometry::GlobalCoordinate;
public:
    // export permeability type
    using PermeabilityType = Scalar;

    template<class GridData>
    SolidSpatialParams(std::shared_ptr<const GridGeometry> fvGridGeometry, const GridData& gridData)
    : ParentType(fvGridGeometry)
    {
        throatCenter_.resize(fvGridGeometry->gridView().size(0));
        poreExtendedRadius_.resize(fvGridGeometry->gridView().size(1));

        for (const auto& element : elements(fvGridGeometry->gridView()))
        {
            const auto eIdx = fvGridGeometry->elementMapper().index(element);
            const auto& params = gridData.parameters(element);
            const auto posX = params[gridData.parameterIndex("ThroatCenterX")];
            const auto posY = params[gridData.parameterIndex("ThroatCenterY")];
            const auto posZ = params[gridData.parameterIndex("ThroatCenterZ")];
            throatCenter_[eIdx] = GlobalPosition{posX, posY, posZ};
        }

        for (const auto& vertex : vertices(fvGridGeometry->gridView()))
        {
            const auto vIdx = fvGridGeometry->vertexMapper().index(vertex);
            poreExtendedRadius_[vIdx] = gridData.getParameter(vertex, "PoreExtendedRadius");
        }
    }

    const GlobalPosition& throatCenter(const std::size_t eIdx) const
    { return throatCenter_[eIdx]; }

    template<class Element, class SubControlVolume, class ElementSolutionVector>
    Scalar poreExtendedRadius(const Element& element,
                              const SubControlVolume& scv,
                              const ElementSolutionVector& elemSol) const
    { return poreExtendedRadius_[scv.dofIndex()]; }

    Scalar poreExtendedRadius(std::size_t dofIdx) const
    { return poreExtendedRadius_[dofIdx]; }

private:
    std::vector<GlobalPosition> throatCenter_;
    std::vector<Scalar> poreExtendedRadius_;
};

} // end namespace Dumux

#endif
