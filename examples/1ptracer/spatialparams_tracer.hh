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

#ifndef DUMUX_TRACER_TEST_SPATIAL_PARAMS_HH
#define DUMUX_TRACER_TEST_SPATIAL_PARAMS_HH

// ## Parameter distributions (`spatialparams_tracer.hh`)
//
//
// In this file, we define spatial properties of the porous medium such as the permeability
// and the porosity in various functions for the tracer problem. Furthermore, spatially-dependent
// properties of the tracer fluid system are defined as well as functions related to setting and retrieving
// the volume fluxes calculated from the solution of the 1p problem.
//
// ### Include files
// We use the properties for porous medium flow models, declared in the file `properties.hh`.
#include <dumux/porousmediumflow/properties.hh>
// As in the 1p spatialparams, we inherit from the spatial parameters for single-phase models using finite volumes, which we include here.
#include <dumux/material/spatialparams/fv1p.hh>

// ### The spatial parameters class
//
// In the `TracerTestSpatialParams` class, we define all functions needed to describe
// spatially-dependent parameters for the tracer_problem.
// [[codeblock]]
namespace Dumux {

template<class GridGeometry, class Scalar>
class TracerTestSpatialParams
: public FVSpatialParamsOneP<GridGeometry, Scalar,
                             TracerTestSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParamsOneP<GridGeometry, Scalar,
                                           TracerTestSpatialParams<GridGeometry, Scalar>>;

    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Dune::FieldVector<Scalar, dimWorld>;

public:

    TracerTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry) {}
    // [[/codeblock]]
    //
    // #### Properties of the porous matrix
    // We define the same porosity for the whole domain as in the 1p spatialparams.
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.2; }

    // We do not consider dispersivity for the tracer transport. Thus, we set the
    // dispersivity coefficient to zero.
    template<class ElementSolution>
    Scalar dispersivity(const Element &element,
                        const SubControlVolume& scv,
                        const ElementSolution& elemSol) const
    { return 0; }

    // #### Properties of the fluid system
    // In the following, we define fluid properties that are spatial parameters in the tracer model.
    // They can possible vary in space but are usually constants.
    // Furthermore, spatially constant values of the fluid system are defined in the `TracerFluidSystem`
    // class in `proerties_tracer.hh`. We define the fluid density to a constant value of 1000 $`\frac{kg}{m^3}`$.
    Scalar fluidDensity(const Element &element,
                        const SubControlVolume& scv) const
    { return 1000; }

    // The following functions define the molar mass of the fluid in function of the
    // elements of the computational grid and the position in the domain.
    // [[codeblock]]
    // This interface defines the fluid molar mass within the sub-control volume `scv`.
    Scalar fluidMolarMass(const Element& element,
                          const SubControlVolume& scv) const
    { return fluidMolarMassAtPos(scv.dofPosition()); }

    // This interface defines the fluid molar mass depending on the position in the domain.
    Scalar fluidMolarMassAtPos(const GlobalPosition& globalPos) const
    { return 18.0; }
    // [[/codeblock]]
    //
    // #### The volume fluxes
    // The following function returns the volume flux across the given sub-control volume face `scvf`.
    // This flux is obtained from the vector `volumeFlux_` that contains the fluxes across al sub-control
    // volume faces of the discretization. This vector can be set using the `setVolumeFlux` function.
    template<class ElementVolumeVariables>
    Scalar volumeFlux(const Element &element,
                      const FVElementGeometry& fvGeometry,
                      const ElementVolumeVariables& elemVolVars,
                      const SubControlVolumeFace& scvf) const
    {
        return volumeFlux_[scvf.index()];
    }

    // This function allows setting the volume fluxes for all sub-control volume faces of the discretization.
    // This is used in the main function after these fluxes have been based on the pressure solution obtained
    // with the single-phase model.
    void setVolumeFlux(const std::vector<Scalar>& f)
    { volumeFlux_ = f; }

    // The remainder of the class contains the private data members, which in this case
    // are only the volume fluxes across the sub-control volume faces of the discretization.
    // [[codeblock]]
private:
    std::vector<Scalar> volumeFlux_;
}; // end class definition TracerTestSpatialParams
} // end namespace Dumux
// [[/codeblock]]

#endif
