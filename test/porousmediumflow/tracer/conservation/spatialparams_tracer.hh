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
 * \brief spatial parameters for the tracer conservation test
 */
#ifndef DUMUX_TEST_SPATIAL_PARAMS_TRACER_HH
#define DUMUX_TEST_SPATIAL_PARAMS_TRACER_HH

#include <dumux/material/spatialparams/fv.hh>

namespace Dumux {

/*!
 * \ingroup TracerTests
 * \brief spatial parameters for the tracer conservation test
 */
template<class GridGeometry, class Scalar>
class TracerConservationTestSpatialParams
: public FVSpatialParams<GridGeometry, Scalar,
                         TracerConservationTestSpatialParams<GridGeometry, Scalar>>
{
    using ParentType = FVSpatialParams<GridGeometry, Scalar,
                                       TracerConservationTestSpatialParams<GridGeometry, Scalar>>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:

    TracerConservationTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry) {}

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

    Scalar fluidDensity(const Element &element,
                        const SubControlVolume& scv) const
    { return density_[this->gridGeometry().elementMapper().index(element)];; }

    void setDensity(const std::vector<Scalar>& d)
    { density_ = d; }

    Scalar fluidMolarMass(const Element &element,
                          const SubControlVolume& scv) const
    { return 0.01; }

    Scalar fluidMolarMass(const GlobalPosition &globalPos) const
    { return 0.01; }

    template<class ElementVolumeVariables>
    Scalar volumeFlux(const Element &element,
                      const FVElementGeometry& fvGeometry,
                      const ElementVolumeVariables& elemVolVars,
                      const SubControlVolumeFace& scvf) const
    { return volumeFlux_[scvf.index()]; }

    void setVolumeFlux(const std::vector<Scalar>& f)
    { volumeFlux_ = f; }

    Scalar saturation(const Element &element,
                      const SubControlVolume& scv) const
    { return saturation_[scv.dofIndex()]; }

    void setSaturation(const std::vector<Scalar>& s)
    { saturation_ = s; }

private:
    std::vector<Scalar> volumeFlux_;
    std::vector<Scalar> density_;
    std::vector<Scalar> saturation_;
};

} // end namespace Dumux

#endif
