// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_PNM_ONEP_PERMEABILITY_UPSCALING_PROBLEM_HH
#define DUMUX_PNM_ONEP_PERMEABILITY_UPSCALING_PROBLEM_HH

/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief The problem class defining initial and boundary conditions for the single-phase
 *        Darcy flow solve used to extract the apparent permeability of one fluid at a given
 *        saturation state.
 */
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/porousmediumflow/problem.hh>

#include "transmissibilityhelper.hh"

namespace Dumux {

/*!
 * \ingroup PoreNetworkModels
 * \brief The problem class for the single-phase upscaling solve.
 *
 * A fixed pressure gradient is applied between the inlet and outlet pores of the network
 * (identified either via pore labels or via the network's bounding box); all other boundaries
 * are no-flow. The resulting single-phase mass flux is used by UpscalingHelper to compute the
 * apparent permeability of the currently selected fluid at the current saturation state.
 */
template<class TypeTag>
class UpscalingProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using BoundaryTypes = Dumux::BoundaryTypes<PrimaryVariables::size()>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

public:
    UpscalingProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        // the applied pressure gradient
        pressureGradient_ = getParam<Scalar>("Problem.PressureGradient", 1e5);

        // We can either use pore labels (given in the grid file) to identify inlet and outlet pores
        // or use the network's bounding box to find these pores automatically. Using labels is usually much
        // more accurate, so this is the default here.
        useLabels_ = getParam<bool>("Problem.UseLabels", true);

        // an epsilon value for the floating point comparisons to determine inlet/outlet pores
        eps_ = getParam<Scalar>("Problem.Epsilon", 1e-7);
    }

    // We need to specify a constant temperature for our isothermal problem.
    // Fluid properties that depend on temperature will be calculated with this value.
    Scalar temperature() const
    { return 283.15; }

    // Set the pressure gradient to be applied to the network
    void setPressureGradient(Scalar pressureGradient)
    { pressureGradient_ = pressureGradient; }

    // This function is used to define the type of boundary conditions used depending on the location.
    // Here, we use Dirichlet boundary conditions (fixed pressures) at the inlet and outlet. Note that
    // the PNM does not support Neumann boundaries: to specify a certain mass flux on a boundary, we
    // would have to use a source term on the boundary pores (not done in this example).
    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolume& scv) const
    {
        BoundaryTypes bcTypes;

        // fix the pressure at the inlet and outlet pores
        if (isInletPore_(scv) || isOutletPore_(scv))
            bcTypes.setAllDirichlet();

        return bcTypes;
    }

    // The following function specifies the values on Dirichlet boundaries (pressures).
    // We set 0 Pa at the outlet and a value based on the given pressure gradient
    // and the length of the domain at the inlet.
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

    // Return the fluid mass density.
    Scalar density() const
    {
        return FluidSystem::density(0.0, 0.0); // dummy values for pressure and temperature
    }

    // Return the fluid dynamic viscosity.
    Scalar dynamicViscosity() const
    {
        return FluidSystem::viscosity(0.0, 0.0); // dummy values for pressure and temperature
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

    // Imports the per-throat phase conductances precomputed by the quasi-static two-phase
    // model (see problem_static.hh) for the current saturation step.
    template<class Transmissibility>
    void importTransmissibility(const Transmissibility& transmissibility)
    {
        Dumux::PoreNetwork::SinglePhaseTransmissibility<Scalar>::importTransmissibility(transmissibility);
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

    Scalar eps_;
    Scalar pressureGradient_;
    int direction_;
    GlobalPosition length_;
    bool useLabels_;
};

} // end namespace Dumux

#endif
