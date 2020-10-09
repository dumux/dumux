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
    static constexpr int dimworld = GridGeometry::GridView::dimensionworld;

public:
    // In the constructor, we obtain a number of parameters, related to fluid
    // properties and boundary conditions, from the input file.
    // [[codeblock]]
    template<class SpatialParams>
    UpscalingProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<SpatialParams> spatialParams)
    : ParentType(gridGeometry, spatialParams)
    {
        // the applied pressure gradient
        pressureGradient_ = getParam<Scalar>("Problem.PressureGradient");

        for (int i = 0; i < dimworld; ++i)
            lenght_[i] = this->gridGeometry().bBoxMax()[i] - this->gridGeometry().bBoxMin()[i];

        eps_ = getParam<Scalar>("Problem.Epsilon", 1e-7);
        useLabels_ = getParam<bool>("Problem.UseLabels", true);
    }
    // [[/codeblock]]

    // We need to specify a constant temperature for our isothermal problem.
    // Fluid properties that depend on temperature will be calculated with this value.
    Scalar temperature() const
    { return 283.15; }

    // #### Specify the types of boundary conditions
    // This function is used to define the type of boundary conditions used depending on the location.
    // Two types of boundary  conditions can be specified: Dirichlet or Neumann boundary condition.
    // On a Dirichlet boundary, the values of the primary variables need to be fixed. On a Neumann
    // boundary condition, values for derivatives need to be fixed. Here, we use Dirichlet boundary
    // conditions on all boundaries.
    BoundaryTypes boundaryTypes(const Element &element, const SubControlVolume& scv) const
    {
        BoundaryTypes bcTypes;

        if (isInletPore_(scv)|| isOutletPore_(scv))
        {
            bcTypes.setAllDirichlet();
        }
        else // neuman for the remaining boundaries
            bcTypes.setAllNeumann();

        return bcTypes;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param vertex The vertex (pore body) for which the condition is evaluated
     *
     */
     PrimaryVariables dirichlet(const Element &element,
                                const SubControlVolume &scv) const
     {
        PrimaryVariables values(0.0);

        if (isInletPore_(scv))
            values[Indices::pressureIdx] = pressureGradient_ * lenght_[direction_];
        else
            values[Indices::pressureIdx] = 0.0;

        return values;
    }

    /*!
     * \brief Sets the current direction in which the pressure gradient is applied
     * \param directionIdx The index of the direction (0:x, 1:y, 2:z)
     */
    void setDirection(int directionIdx)
    { direction_ = directionIdx; }


private:

    bool isInletPore_(const SubControlVolume& scv) const
    {
        const auto poreLabel = this->gridGeometry().poreLabel(scv.dofIndex());

        if (poreLabel < 0)
            return false;

        if (useLabels_)
            return poreLabel == 1 + 2*direction_;
        else
            return scv.dofPosition()[direction_] < this->gridGeometry().bBoxMin()[direction_] + eps_;
    }

    bool isOutletPore_(const SubControlVolume& scv) const
    {
        const auto poreLabel = this->gridGeometry().poreLabel(scv.dofIndex());

        if (poreLabel < 0)
            return false;

        if (useLabels_)
            return poreLabel == 2 + 2*direction_;
        else
            return scv.dofPosition()[direction_] > this->gridGeometry().bBoxMax()[direction_] - eps_;
    }

    // private data members
    Scalar eps_;
    Scalar pressureGradient_;
    int direction_;
    std::array<Scalar, dimworld> lenght_;
    bool useLabels_;
};

} // end namespace Dumux
// [[/codeblock]]
// [[/content]]
#endif
