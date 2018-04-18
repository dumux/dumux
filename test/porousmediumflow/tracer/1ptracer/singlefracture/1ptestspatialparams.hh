// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup TracerTests
 * \brief The spatial params the incompressible test
 */
#ifndef DUMUX_INCOMPRESSIBLE_ONEP_TEST_SPATIAL_PARAMS_HH
#define DUMUX_INCOMPRESSIBLE_ONEP_TEST_SPATIAL_PARAMS_HH

#include <random>
#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux
{
/*!
 * \ingroup TracerTests
 * \brief The spatial params the incompressible test
 */
template<class TypeTag>
class OnePTestSpatialParams : public FVSpatialParamsOneP<TypeTag>
{
    using ParentType = FVSpatialParamsOneP<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using GridCreator = typename GET_PROP_TYPE(TypeTag, GridCreator);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    enum LayerMarkers
    {
        lowerLayer = 1,
        upperLayer = 2,
        fracture = 3
    };

    using PermeabilityType = Scalar;
    OnePTestSpatialParams(const Problem& problem) : ParentType(problem)
    {
        bottomLayerK_ = getParam<Scalar>("SpatialParams.BottomLayerPermeability");
        upperLayerK_ = getParam<Scalar>("SpatialParams.UpperLayerPermeability");
        fractureK_ = getParam<Scalar>("SpatialParams.FracturePermeability");
        bottomLayerPhi_ = getParam<Scalar>("SpatialParams.BottomLayerPorosity");
        upperLayerPhi_ = getParam<Scalar>("SpatialParams.UpperLayerPorosity");
        fracturePhi_ = getParam<Scalar>("SpatialParams.FracturePorosity");
    }

    //! Function for defining the (intrinsic) permeability \f$[m^2]\f$.
    const Scalar permeability(const Element &element,
                              const SubControlVolume &scv,
                              const ElementSolutionVector &elemSol) const
    { return permeability(element); }

    const Scalar permeability(const Element &element) const
    {
        const auto eIdx = this->problem().fvGridGeometry().elementMapper().index(element);
        const auto marker = GridCreator::getElementDomainMarker(eIdx);
        if (marker == LayerMarkers::lowerLayer)
            return bottomLayerK_;
        else if (marker == LayerMarkers::upperLayer)
            return upperLayerK_;
        else if (marker == LayerMarkers::fracture)
            return fractureK_;
        else
            DUNE_THROW(Dune::InvalidStateException, "Unknown element marker");
    }

    //! Function for defining the porosity.
    Scalar porosity(const Element &element,
                    const SubControlVolume &scv,
                    const ElementSolutionVector &elemSol) const
    { return porosity(element); }

    Scalar porosity(const Element& element) const
    {
        const auto eIdx = this->problem().fvGridGeometry().elementMapper().index(element);
        const auto marker = GridCreator::getElementDomainMarker(eIdx);
        if (marker == LayerMarkers::lowerLayer)
            return bottomLayerPhi_;
        else if (marker == LayerMarkers::upperLayer)
            return upperLayerPhi_;
        else if (marker == LayerMarkers::fracture)
            return fracturePhi_;
        else
            DUNE_THROW(Dune::InvalidStateException, "Unknown element marker");
    }

private:
    Scalar bottomLayerK_;
    Scalar upperLayerK_;
    Scalar fractureK_;
    Scalar bottomLayerPhi_;
    Scalar upperLayerPhi_;
    Scalar fracturePhi_;
};
} // end namespace Dumux

#endif
