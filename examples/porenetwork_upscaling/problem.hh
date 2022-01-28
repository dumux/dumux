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

#ifndef DUMUX_PNM_ONEP_PERMEABILITY_UPSCALING_PROBLEM_HH
#define DUMUX_PNM_ONEP_PERMEABILITY_UPSCALING_PROBLEM_HH

// ## The problem class (`problem.hh`)
// This file contains the __problem class__ which defines the initial and boundary
// conditions for the single-phase flow simulation.
// [[content]]
// ### Includes
#include <dumux/common/boundarytypes.hh> // for `BoundaryTypes`
#include <dumux/common/properties.hh> // for `GetPropType`
#include <dumux/common/parameters.hh> // for `getParam`
#include <dumux/porousmediumflow/problem.hh>  // for `PorousMediumFlowProblem`

// ### The problem class
// We enter the problem class where all necessary boundary conditions and initial conditions are set for our simulation.
// As this is a porous medium flow problem, we inherit from the base class `PorousMediumFlowProblem`.
namespace Dumux {

template<class TypeTag>
class UpscalingProblem : public PorousMediumFlowProblem<TypeTag>
{
    // [[details]] convenience aliases
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using BoundaryTypes = Dumux::BoundaryTypes<PrimaryVariables::size()>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    // [[/details]]

    // ### The constructor of our problem.
    // [[codeblock]]
public:
    template<class SpatialParams>
    UpscalingProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<SpatialParams> spatialParams)
    : ParentType(gridGeometry, spatialParams)
    {
        // the applied pressure gradient
        pressureGradient_ = getParam<Scalar>("Problem.PressureGradient");

        // We can either use pore labels (given in the grid file) to identify inlet and outlet pores
        // or use the network's bounding box to find these pores automatically. Using labels is usually much
        // more accurate, so this is the default here.
        useLabels_ = getParam<bool>("Problem.UseLabels", true);

        // an epsilon value for the bounding box approach
        eps_ = getParam<Scalar>("Problem.Epsilon", 1e-7);
    }
    // [[/codeblock]]

    // #### Boundary conditions
    // This function is used to define the __type of boundary conditions__ used depending on the location.
    // Here, we use Dirichlet boundary conditions (fixed pressures) at the inlet and outlet and Neumann
    // boundary conditions at all remaining boundaries. Note that the PNM only supports Neumann no-flow boundaries.
    // The specify a certain mass flux, we would have to use a source term on the boundary pores (which is not done in this example).
    // [[codeblock]]
    BoundaryTypes boundaryTypes(const Element &element, const SubControlVolume& scv) const
    {
        BoundaryTypes bcTypes;

        // fix the pressure at the inlet and outlet pores
        if (isInletPore_(scv)|| isOutletPore_(scv))
            bcTypes.setAllDirichlet();
        else // Neumann (no-flow) for the remaining boundaries
            bcTypes.setAllNeumann();

        return bcTypes;
    }
    // [[/codeblock]]

    // The following function specifies the __values on Dirichlet boundaries__ (pressures).
    // We set 0 Pa at the outlet and a value based on the given pressure gradient
    // and the length of the domain at the inlet.
    // [[codeblock]]
     PrimaryVariables dirichlet(const Element& element,
                                const SubControlVolume& scv) const
     {
        PrimaryVariables values(0.0);

        if (isInletPore_(scv))
            values[Indices::pressureIdx] = pressureGradient_ * length_[direction_];
        else
            values[Indices::pressureIdx] = 0.0;

        return values;
    }
    // [[/codeblock]]

    // #### Upscaling

    // [[details]] auxiliary functions needed for the upscaling process
    // [[codeblock]]

    // Set the current direction (0:x, 1:y, 2:z) in which the pressure gradient is applied
    void setDirection(int directionIdx)
    { direction_ = directionIdx; }

    // Get the current direction in which the pressure gradient is applied.
    int direction() const
    { return direction_; }

    // Set the side lengths to consider for the upscaling process.
    void setSideLengths(const GlobalPosition& sideLengths)
    { length_ = sideLengths; }

    // Return the side lengths to consider for the upscaling process.
    const GlobalPosition& sideLengths() const
    { return length_; }

    // Return the liquid mass density.
    Scalar liquidDensity() const
    {
        static const Scalar liquidDensity = getParam<Scalar>("Component.LiquidDensity");
        return liquidDensity;
    }

    // Return the liquid dynamic viscosity.
    Scalar liquidDynamicViscosity() const
    {
        static const Scalar liquidDynamicViscosity = getParam<Scalar>("Component.LiquidKinematicViscosity") * liquidDensity();
        return liquidDynamicViscosity;
    }

    // Return the applied pressure gradient.
    Scalar pressureGradient() const
    { return pressureGradient_; }


    // Return the label of inlet pores assuming a previously set direction.
    int inletPoreLabel() const
    {
        static constexpr std::array<int, 3> label = {1, 3, 5};
        return label[direction_];
    }

    // Return the label of outlet pores assuming a previously set direction.
    int outletPoreLabel() const
    {
        static constexpr std::array<int, 3> label = {2, 4, 6};
        return label[direction_];
    }

private:

    bool isInletPore_(const SubControlVolume& scv) const
    {
        if (useLabels_)
            return inletPoreLabel() == this->gridGeometry().poreLabel(scv.dofIndex());
        else
            return scv.dofPosition()[direction_] < this->gridGeometry().bBoxMin()[direction_] + eps_;
    }

    bool isOutletPore_(const SubControlVolume& scv) const
    {
        if (useLabels_)
            return outletPoreLabel() == this->gridGeometry().poreLabel(scv.dofIndex());
        else
            return scv.dofPosition()[direction_] > this->gridGeometry().bBoxMax()[direction_] - eps_;
    }

    // private data members
    Scalar eps_;
    Scalar pressureGradient_;
    int direction_;
    GlobalPosition length_;
    bool useLabels_;

    // [[/codeblock]]
    // [[/details]]
};

} // end namespace Dumux
// [[/codeblock]]
// [[/content]]
#endif
