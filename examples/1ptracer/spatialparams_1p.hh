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

#ifndef DUMUX_INCOMPRESSIBLE_ONEP_TEST_SPATIAL_PARAMS_HH
#define DUMUX_INCOMPRESSIBLE_ONEP_TEST_SPATIAL_PARAMS_HH

// ## Parameter distributions (`spatialparams_1p.hh`)
//
// This file contains the __spatial parameters class__ which defines the
// distributions for the porous medium parameters permeability and porosity
// over the computational grid
//
// [[content]]
//
// ### Include files
// In this example, we use a randomly generated and element-wise distributed
// permeability field. For this, we use the random number generation facilities
// provided by the C++ standard library.
#include <random>

// We include the spatial parameters class for single-phase models discretized
// by finite volume schemes, from which the spatial parameters defined for this
// example will inherit.
#include <dumux/material/spatialparams/fv1p.hh>
//
// ### The spatial parameters class
//
// In the `OnePTestSpatialParams` class, we define all functions needed to describe
// the porous medium, e.g. porosity and permeability, for the single-phase problem.
// We inherit from the `FVSpatialParamsOneP` class here, which is the base class
// for spatial parameters in the context of single-phase porous medium flow
// applications using finite volume discretization schemes.
// [[codeblock]]
namespace Dumux {

template<class GridGeometry, class Scalar>
class OnePTestSpatialParams
: public FVSpatialParamsOneP<GridGeometry, Scalar,
                             OnePTestSpatialParams<GridGeometry, Scalar>>
{
    // The following convenience aliases will be used throughout this class
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParamsOneP<GridGeometry, Scalar,
                                           OnePTestSpatialParams<GridGeometry, Scalar>>;

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename SubControlVolume::GlobalPosition;

public:
    // The spatial parameters must export the type used to define permeabilities.
    // Here, we are using scalar permeabilities, but tensors are also supported.
    using PermeabilityType = Scalar;
    // [[/codeblock]]
    //
    // #### Generation of the random permeability field
    // We generate a random permeability field upon construction of the spatial parameters class
    // using lognormal distributions. The mean values for the permeability inside and outside of a
    // low-permeable lens (given by the coorinates `lensLowerLeft_` and `lensUpperRight_`) are defined
    // in the variables  `permeabilityLens_` and `permeability_`. The respective values are obtained
    // from the input file making use of the free function `getParam`. We use a standard deviarion
    // of 10% here and compute permeabily values for all elements of the computational grid.
    // [[codeblock]]
    OnePTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry), K_(gridGeometry->gridView().size(0), 0.0)
    {
        // The permeability of the domain and the lens are obtained from the `params.input` file.
        permeability_ = getParam<Scalar>("SpatialParams.Permeability");
        permeabilityLens_ = getParam<Scalar>("SpatialParams.PermeabilityLens");

        // Furthermore, the position of the lens, which is defined by the position of the lower left and the upper right corners, are obtained from the input file.
        lensLowerLeft_ = getParam<GlobalPosition>("SpatialParams.LensLowerLeft");
        lensUpperRight_ =getParam<GlobalPosition>("SpatialParams.LensUpperRight");

        // We generate random fields for the permeability using lognormal distributions,
        // with `permeability_` as mean value and 10 % of it as standard deviation.
        // A separate distribution is used for the lens using `permeabilityLens_`.
        // A permeability value is created for each element of the grid and is stored in the vector `K_`.
        std::mt19937 rand(0);
        std::lognormal_distribution<Scalar> K(std::log(permeability_), std::log(permeability_)*0.1);
        std::lognormal_distribution<Scalar> KLens(std::log(permeabilityLens_), std::log(permeabilityLens_)*0.1);

        // loop over all elements and compute a permeability value
        for (const auto& element : elements(gridGeometry->gridView()))
        {
            const auto eIdx = gridGeometry->elementMapper().index(element);
            const auto globalPos = element.geometry().center();
            K_[eIdx] = isInLens_(globalPos) ? KLens(rand) : K(rand);
        }
    }
    // [[/codeblock]]

    // #### Properties of the porous matrix
    // This function returns the permeability $`[m^2]`$ to be used within a sub-control volume (`scv`) inside the element `element`.
    // One can define the permeability as function of the primary variables on the element, which are given in the provided
    // instance of `ElementSolution`. Here, we use element-wise distributed permeabilities that were randomly generated in
    // the constructor and stored in the vector `K_`(see above).
    template<class ElementSolution>
    const PermeabilityType& permeability(const Element& element,
                                         const SubControlVolume& scv,
                                         const ElementSolution& elemSol) const
    {
        return K_[scv.elementIndex()];
    }

    // We set the porosity $`[-]`$ for the whole domain to a value of $`20 \%`$.
    // Note that in case you want to use solution-dependent porosities, you can
    // use the overload
    // `porosity(const Element&, const SubControlVolume&, const ElementSolution&)`
    // that is defined in the base class `FVSpatialParamsOneP`. Per default, this
    // fowards to the `porosityAtPos` function per default, which we overload here.
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.2; }

    // We reference to the permeability field. This is used in the main function to write an output for the permeability field.
    const std::vector<Scalar>& getKField() const
    { return K_; }

    // The remainder of this class contains a convenient function to determine if
    // a position is inside the lens and defines the data members.
    // [[details]] private data members and member functions
    // [[codeblock]]
private:
    // The following function returns true if a given position is inside the lens.
    // We use an epsilon of 1.5e-7 here for floating point comparisons.
    bool isInLens_(const GlobalPosition& globalPos) const
    {
        for (int i = 0; i < dimWorld; ++i)
        {
            if (globalPos[i] < lensLowerLeft_[i] + 1.5e-7
                || globalPos[i] > lensUpperRight_[i] - 1.5e-7)
                return false;
        }

        return true;
    }

    GlobalPosition lensLowerLeft_, lensUpperRight_;
    Scalar permeability_, permeabilityLens_;
    std::vector<Scalar> K_;
}; // end class definition of OnePTestSpatialParams
} // end namespace Dumux
// [[/codeblock]]
// [[/details]]
// [[/content]]
#endif
