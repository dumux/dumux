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

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

template<class TypeTag>
class UpscalingProblem : public PorousMediumFlowProblem<TypeTag>
{
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

public:
    template<class SpatialParams>
    UpscalingProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<SpatialParams> spatialParams, const std::string& paramGroup)
    : ParentType(gridGeometry, spatialParams, paramGroup)
    {
        // the applied pressure gradient
        pressureGradient_ = getParamFromGroup<Scalar>(paramGroup, "Problem.PressureGradient");

        // We can either use pore labels (given in the grid file) to identify inlet and outlet pores
        // or use the network's bounding box to find these pores automatically. Using labels is usually much
        // more accurate, so this is the default here.
        useLabels_ = getParamFromGroup<bool>(paramGroup, "Problem.UseLabels", true);

        // an epsilon value for the bounding box approach
        eps_ = getParamFromGroup<Scalar>(paramGroup, "Problem.Epsilon", 1e-7);
    }

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
};

} // end namespace Dumux

#endif
