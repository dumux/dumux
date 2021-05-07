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
 * \ingroup TracerTests
 * \brief Definition of the spatial parameters for the tracer problem.
 */

#ifndef DUMUX_TRACER_TEST_SPATIAL_PARAMS_HH
#define DUMUX_TRACER_TEST_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux {

/*!
 * \ingroup TracerTests
 * \brief Definition of the spatial parameters for the tracer problem.
 */
template<class GridGeometry, class Scalar>
class AdvectionDiffusionTestSpatialParams
: public FVSpatialParamsOneP<GridGeometry, Scalar,
                             AdvectionDiffusionTestSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParamsOneP<GridGeometry, Scalar,
                                           AdvectionDiffusionTestSpatialParams<GridGeometry, Scalar>>;

    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Dune::FieldVector<Scalar, dimWorld>;

public:
    AdvectionDiffusionTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

    //! Fluid properties that are spatial parameters in the tracer model
    //! They can possibly vary with space but are usually constants

    //! Fluid density
    Scalar fluidDensity(const Element &element,
                        const SubControlVolume& scv) const
    {
        static const Scalar density = getParam<Scalar>("Problem.Density");
        return density;
    }

    //! Fluid molar mass
    Scalar fluidMolarMass(const Element &element,
                          const SubControlVolume& scv) const
    {
        const auto globalPos = scv.center();
        return fluidMolarMass(globalPos);
    }

    Scalar fluidMolarMass(const GlobalPosition &globalPos) const
    {
        static const Scalar fluidMolarMass = getParam<Scalar>("Problem.FluidMolarMass");
        return fluidMolarMass;
    }

    //! Set the velocity field
    void setVelocityAtFace(const std::vector<GlobalPosition>& v)
    { velocityAtFace_ = v; }

    //! The velocity field
    GlobalPosition velocity(const SubControlVolumeFace& scvf) const
    { return velocityAtFace_[scvf.index()]; }

    //! The volume flux
    template<class ElementVolumeVariables>
    Scalar volumeFlux(const Element &element,
                      const FVElementGeometry& fvGeometry,
                      const ElementVolumeVariables& elemVolVars,
                      const SubControlVolumeFace& scvf) const
    { return velocity(scvf) * scvf.unitOuterNormal() * scvf.area() * elemVolVars[fvGeometry.scv(scvf.insideScvIdx())].extrusionFactor(); }

private:
    std::vector<GlobalPosition> velocityAtFace_;

};

} // end namespace Dumux

#endif
