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

// the header guard
#ifndef DUMUX_INCOMPRESSIBLE_ONEP_TEST_SPATIAL_PARAMS_HH
#define DUMUX_INCOMPRESSIBLE_ONEP_TEST_SPATIAL_PARAMS_HH

// In this file we generate the random permeability field in the constructor of the OnePTestSpatialParams class. Thereafter spatial properties of the porous medium like permeability and porosity are defined in various functions for the 1p problem.
// We want to generate a random permeability field. For this we use a random number generation of the C++ standard library.
#include <random>
// In the file properties.hh all properties are declared.
#include <dumux/porousmediumflow/properties.hh>

// We include the spatial parameters for single-phase, finite volumes from which we will inherit.
#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux {

// In the OnePTestSpatialParams class, we define all functions needed to describe the porous matrix, e.g. define porosity and permeability for the 1p_problem.

template<class FVGridGeometry, class Scalar>
class OnePTestSpatialParams
: public FVSpatialParamsOneP<FVGridGeometry, Scalar,
                             OnePTestSpatialParams<FVGridGeometry, Scalar>>
{
    // We introduce using declarations that are derived from the property system which we need in this class.
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParamsOneP<FVGridGeometry, Scalar,
                                           OnePTestSpatialParams<FVGridGeometry, Scalar>>;

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename SubControlVolume::GlobalPosition;

public:
    using PermeabilityType = Scalar;
    OnePTestSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry), K_(fvGridGeometry->gridView().size(0), 0.0)
    {
        // ### Generation of the random permeability field
        // We get the permeability of the domain and the lens from the params.input file.
        permeability_ = getParam<Scalar>("SpatialParams.Permeability");
        permeabilityLens_ = getParam<Scalar>("SpatialParams.PermeabilityLens");

        // Further, we get the position of the lens, which is defined by the position of the lower left and the upper right corner.
        lensLowerLeft_ = getParam<GlobalPosition>("SpatialParams.LensLowerLeft");
        lensUpperRight_ =getParam<GlobalPosition>("SpatialParams.LensUpperRight");

        // We generate random fields for the permeability using a lognormal distribution, with the permeability_ as mean value and 10 % of it as standard deviation. A seperate distribution is used for the lens using permeabilityLens_.
        std::mt19937 rand(0);
        std::lognormal_distribution<Scalar> K(std::log(permeability_), std::log(permeability_)*0.1);
        std::lognormal_distribution<Scalar> KLens(std::log(permeabilityLens_), std::log(permeabilityLens_)*0.1);
        for (const auto& element : elements(fvGridGeometry->gridView()))
        {
            const auto eIdx = fvGridGeometry->elementMapper().index(element);
            const auto globalPos = element.geometry().center();
            K_[eIdx] = isInLens_(globalPos) ? KLens(rand) : K(rand);
        }
    }

    // ### Properties of the porous matrix
    // We define the (intrinsic) permeability $`[m^2]`$ using the generated random permeability field. In this test, we use element-wise distributed permeabilities.
    template<class ElementSolution>
    const PermeabilityType& permeability(const Element& element,
                                         const SubControlVolume& scv,
                                         const ElementSolution& elemSol) const
    {
        return K_[scv.dofIndex()];
    }


    // We set the porosity $`[-]`$ for the whole domain.
    Scalar porosityAtPos(const GlobalPosition &globalPos) const
    { return 0.2; }

    // We reference to the permeability field. This is used in the main function to write an output for the permeability field.
    const std::vector<Scalar>& getKField() const
    { return K_; }

private:
    // We have a convenience definition of the position of the lens.
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

    Scalar permeability_, permeabilityLens_;

    std::vector<Scalar> K_;

    static constexpr Scalar eps_ = 1.5e-7;
};

} // end namespace Dumux

#endif
