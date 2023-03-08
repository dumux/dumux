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
 * \brief The spatial parameters class for the Hhat problem with multiple solid spheres
 */

#ifndef DUMUX_TEST_MULTIDOMAIN_DUALNETWORK_SPATIALPARAMS_HH
#define DUMUX_TEST_MULTIDOMAIN_DUALNETWORK_SPATIALPARAMS_HH

#include <dumux/porenetwork/solidenergy/spatialparams.hh>
#include <dumux/porenetwork/1p/spatialparams.hh>

namespace Dumux {

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
