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

#ifndef DUMUX_TWOP_INCOMPRESSIBLE_EXAMPLE_SPATIAL_PARAMS_HH
#define DUMUX_TWOP_INCOMPRESSIBLE_EXAMPLE_SPATIAL_PARAMS_HH

// ## The file `spatialparams.hh`
// [[content]]
//
// ### Includes
// We include the basic spatial parameters for finite volumes file from which we will inherit
#include <dumux/material/spatialparams/fv.hh>

// We include all laws which are needed to define the interaction between the solid matrix and the fluids, e.g. laws for capillary pressure saturation relationships.
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>

// ### The spatial parameters class
// In the TwoPTestSpatialParams class we define all functions needed to describe the porous matrix, e.g. porosity and permeability. We inherit from the `FVSpatialParams` class, which is the base class for multiphase porous medium flow applications.
// [[codeblock]]
namespace Dumux {

template<class GridGeometry, class Scalar>
class TwoPTestSpatialParams
: public FVSpatialParams<GridGeometry, Scalar, TwoPTestSpatialParams<GridGeometry, Scalar>>
{
    //we introduce using declarations that are derived from the property system which we need in this class
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using ThisType = TwoPTestSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVSpatialParams<GridGeometry, Scalar, ThisType>;

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PcKrSwCurve = FluidMatrix::VanGenuchtenDefault<Scalar>;

public:
    using PermeabilityType = Scalar;
    // [[/codeblock]]

    //Here, we get parameters for the position of the lens and porosity and permeability from the input file. Additionally, we set the parameters for the Van-Genuchten relationship.
    // [[codeblock]]
    TwoPTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , pcKrSwLens_("SpatialParams.Lens") // read params from input file
    , pcKrSwOuter_("SpatialParams.Outer") // read params from input file
    {
        // We get the position of the lens from the params.input file.
        // The lens is defined by the position of the lower left and the upper right corner
        lensLowerLeft_ = getParam<GlobalPosition>("SpatialParams.LensLowerLeft");
        lensUpperRight_ = getParam<GlobalPosition>("SpatialParams.LensUpperRight");

        // Here, we get the permeabilities from the params.input file.
        // In case that no parameter is set, the default parameters (9.05e-12 and 4.6e-10) are used
        lensK_ = getParam<Scalar>("SpatialParams.Lens.Permeability", 9.05e-12);
        outerK_ = getParam<Scalar>("SpatialParams.Outer.Permeability", 4.6e-10);
    }
     // [[/codeblock]]

    // We define the (intrinsic) permeability $`[m^2]`$. In this test, we use element-wise distributed permeabilities.
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        if (isInLens_(element.geometry().center()))
            return lensK_;
        return outerK_;
    }

    // We set the porosity $`[-]`$ depending on the position
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
         if (isInLens_(globalPos))
            return 0.2;
        return 0.4;
    }

    // We set the parameter object for the Van Genuchten material law.
    template<class ElementSolution>
    auto fluidMatrixInteraction(const Element& element,
                                const SubControlVolume& scv,
                                const ElementSolution& elemSol) const
    {
        if (isInLens_(element.geometry().center()))
            return makeFluidMatrixInteraction(pcKrSwLens_);
        return makeFluidMatrixInteraction(pcKrSwOuter_);
    }


     // Here we can define which phase is to be considered as the wetting phase. Here the wetting phase is the the first phase of the fluidsystem. In this case that is water.
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    {  return FluidSystem::phase0Idx; }

    // The remainder of this class contains a convenient function to determine if
    // a position is inside the lens and defines the data members.
    // [[details]] private data members and member functions
    // [[codeblock]]
private:
    // we have a convenience definition of the position of the lens
    bool isInLens_(const GlobalPosition &globalPos) const
    {
        for (int i = 0; i < dimWorld; ++i) {
            if (globalPos[i] < lensLowerLeft_[i] + eps_ || globalPos[i] > lensUpperRight_[i] - eps_)
                return false;
        }
        return true;
    }

    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;

    Scalar lensK_;
    Scalar outerK_;

    const PcKrSwCurve pcKrSwLens_;
    const PcKrSwCurve pcKrSwOuter_;

    static constexpr Scalar eps_ = 1.5e-7;
};

} // end namespace Dumux
// [[/codeblock]]
// [[/details]]
// [[/content]]
#endif
