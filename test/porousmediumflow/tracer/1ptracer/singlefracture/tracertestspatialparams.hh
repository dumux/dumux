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
 * \brief Definition of the spatial parameters for the tracer problem
 */
#ifndef DUMUX_TRACER_TEST_SPATIAL_PARAMS_HH
#define DUMUX_TRACER_TEST_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux
{
/*!
 * \ingroup TracerTests
 * \brief Definition of the spatial parameters for the tracer problem
 */
template<class TypeTag>
class TracerTestSpatialParams : public FVSpatialParamsOneP<TypeTag>
{
    using ParentType = FVSpatialParamsOneP<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using GridCreator = typename GET_PROP_TYPE(TypeTag, GridCreator);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Dune::FieldVector<Scalar, dimWorld>;

public:

    enum LayerMarkers
    {
        lowerLayer = 1,
        upperLayer = 2,
        fracture = 3
    };

    TracerTestSpatialParams(const Problem& problem) : ParentType(problem)
    {
        bottomLayerPhi_ = getParam<Scalar>("SpatialParams.BottomLayerPorosity");
        upperLayerPhi_ = getParam<Scalar>("SpatialParams.UpperLayerPorosity");
        fracturePhi_ = getParam<Scalar>("SpatialParams.FracturePorosity");
    }

    //! Function for defining the porosity.
    Scalar porosity(const Element &element,
                    const SubControlVolume &scv,
                    const ElementSolutionVector &elemSol) const
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

    //! Define the dispersivity
    Scalar dispersivity(const Element &element,
                        const SubControlVolume& scv,
                        const ElementSolutionVector& elemSol) const
    { return 0; }

    //! Fluid properties that are spatial params in the tracer model
    //! They can possible vary with space but are usually constants

    //! fluid density
    static constexpr Scalar fluidDensity() { return 1000.0; }
    Scalar fluidDensity(const Element &element, const SubControlVolume& scv) const { return fluidDensity(); }

    //! fluid molar mass
    static constexpr Scalar fluidMolarMass() { return 18.0e-3; }
    Scalar fluidMolarMass(const Element &element, const SubControlVolume& scv) const { return fluidMolarMass(); }
    Scalar fluidMolarMass(const GlobalPosition &globalPos) const { return fluidMolarMass(); }

    //! velocity field
    Scalar volumeFlux(const Element &element,
                      const FVElementGeometry& fvGeometry,
                      const ElementVolumeVariables& elemVolVars,
                      const SubControlVolumeFace& scvf) const
    { return volumeFlux_[scvf.index()]; }

    Scalar volumeFlux(const SubControlVolumeFace& scvf) const
    { return volumeFlux_[scvf.index()]; }

    void setVolumeFlux(const std::vector<Scalar>& f)
    { volumeFlux_ = f; }

private:
    Scalar bottomLayerPhi_;
    Scalar upperLayerPhi_;
    Scalar fracturePhi_;
    std::vector<Scalar> volumeFlux_;
};

} // end namespace Dumux

#endif
