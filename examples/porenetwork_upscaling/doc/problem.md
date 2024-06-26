<!-- Important: This file has been automatically generated by generate_example_docs.py. Do not edit this file directly! -->


| [:arrow_left: Back to the main documentation](../README.md) | [:arrow_right: Continue with part 2](main.md) |
|---|---:|

# Part 1: Simulation setup

The code documentation is structured as follows:

[[_TOC_]]


## Compile-time settings (`properties.hh`)
This file defines the `TypeTag` used for the simulation in this example, for
which we specialize a number of compile-time `properties`.

<details open>
<summary><b>Click to hide/show the file documentation</b> (or inspect the [source code](../properties.hh))</summary>

### Includes
<details><summary> Click to show includes</summary>

```cpp

#include <dune/foamgrid/foamgrid.hh> // for `Dune::FoamGrid`
```

The `OneP` type tag specializes most of the `properties` required for single-
phase flow simulations in DuMu<sup>x</sup>. We will use this in the following to inherit the
respective properties, and subsequently specialize those properties for our
type tag, which we want to modify or for which no meaningful default can be set.

```cpp
#include <dumux/porenetwork/1p/model.hh>// for `TTag::PNMOneP`
```

The class that contains a collection of single-phase flow throat transmissibilities
among them the transmisibility model to be used can be specified in AdvectionType class

```cpp
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility1p.hh>
```

The class that provides specializations for both creeping and non-creeping advection types.

```cpp
#include <dumux/flux/porenetwork/advection.hh>
```

The local residual for incompressible flow is included.
The one-phase flow model (included above) uses a default implementation of the
local residual for single-phase flow. However, in this example we are using an
incompressible fluid phase. Therefore, we are including the specialized local
residual which contains functionality to analytically compute the entries of
the Jacobian matrix. We will use this in the main file.

```cpp
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>
```

We will use a single liquid phase consisting of a component with constant fluid properties.

```cpp
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
```

The classes that define the problem and parameters used in this simulation

```cpp
#include "problem.hh"
#include "spatialparams.hh"
```

</details>

### `TypeTag` definition
Two type tags for our simulation are defined, one for creeping flow (`PNMUpscalingCreepingFlow`) and another for non-creeping flow (`PNMUpscalingNonCreepingFlow`),
which inherit properties from the single-phase pore network model (`ONMOneP`). The non-creeping flow inherits
all properties from the creeping flow simulation but overrides the `AdvectionType` property.

```cpp
namespace Dumux::Properties {
namespace TTag {
struct PNMUpscalingCreepingFlow { using InheritsFrom = std::tuple<PNMOneP>; };
struct PNMUpscalingNonCreepingFlow { using InheritsFrom = std::tuple<PNMUpscalingCreepingFlow>; };
}
```

### Property specializations

In the following piece of code, mandatory `properties` for which no meaningful
default can be set, are specialized for our type tag `PNMUpscaling`.

```cpp
// We use `dune-foamgrid`, which is especially tailored for 1D networks.
template<class TypeTag>
struct Grid<TypeTag, TTag::PNMUpscalingCreepingFlow>
{ using type = Dune::FoamGrid<1, 3>; };

// The problem class specifying initial and boundary conditions:
template<class TypeTag>
struct Problem<TypeTag, TTag::PNMUpscalingCreepingFlow>
{ using type = UpscalingProblem<TypeTag>; };

//! The spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PNMUpscalingCreepingFlow>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = PoreNetwork::UpscalingSpatialParams<GridGeometry, Scalar>;
};

//! The advection type for creeping flow
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::PNMUpscalingCreepingFlow>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using TransmissibilityLaw = PoreNetwork::TransmissibilityPatzekSilin<Scalar, true/*considerPoreBodyResistance*/>;
public:
    using type = PoreNetwork::CreepingFlow<Scalar, TransmissibilityLaw>;
};

//! The advection type for non-creeping flow (includes model for inertia effects)
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::PNMUpscalingNonCreepingFlow>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using TransmissibilityLaw = PoreNetwork::TransmissibilityPatzekSilin<Scalar, true/*considerPoreBodyResistance*/>;
public:
    using type = PoreNetwork::NonCreepingFlow<Scalar, TransmissibilityLaw>;
};

// We use a single liquid phase consisting of a component with constant fluid properties.
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PNMUpscalingCreepingFlow>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};
```

Moreover, here we use a local residual specialized for incompressible flow
that contains functionality related to analytic differentiation.

```cpp
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::PNMUpscalingCreepingFlow>
{ using type = OnePIncompressibleLocalResidual<TypeTag>; };

} // end namespace Dumux::Properties
```


</details>



## The problem class (`problem.hh`)
This file contains the __problem class__ which defines the initial and boundary
conditions for the single-phase flow simulation.

<details open>
<summary><b>Click to hide/show the file documentation</b> (or inspect the [source code](../problem.hh))</summary>

### Includes

```cpp
#include <dumux/common/boundarytypes.hh> // for `BoundaryTypes`
#include <dumux/common/properties.hh> // for `GetPropType`
#include <dumux/common/parameters.hh> // for `getParam`
#include <dumux/porousmediumflow/problem.hh>  // for `PorousMediumFlowProblem`
```

### The problem class
We enter the problem class where all necessary boundary conditions and initial conditions are set for our simulation.
As this is a porous medium flow problem, we inherit from the base class `PorousMediumFlowProblem`.

```cpp
namespace Dumux {

template<class TypeTag>
class UpscalingProblem : public PorousMediumFlowProblem<TypeTag>
{
```

<details><summary> Click to show convenience aliases</summary>

```cpp
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using BoundaryTypes = Dumux::BoundaryTypes<PrimaryVariables::size()>;
```

</details>

#### The constructor of our problem.

```cpp
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
```

#### Temperature
We need to specify a constant temperature for our isothermal problem.
Fluid properties that depend on temperature will be calculated with this value.

```cpp
    Scalar temperature() const
    { return 283.15; }
```

Set the pressure gradient to be applied to the network

```cpp
    void setPressureGradient(Scalar pressureGradient)
    { pressureGradient_ = pressureGradient; }
```

#### Boundary conditions
This function is used to define the __type of boundary conditions__ used depending on the location.
Here, we use Dirichlet boundary conditions (fixed pressures) at the inlet and outlet. Note that the PNM does not support Neumann boundaries.
To specify a certain mass flux on a boundary, we would have to use a source term on the boundary pores (which is not done in this example).

```cpp
    BoundaryTypes boundaryTypes(const Element &element, const SubControlVolume& scv) const
    {
        BoundaryTypes bcTypes;

        // fix the pressure at the inlet and outlet pores
        if (isInletPore_(scv)|| isOutletPore_(scv))
            bcTypes.setAllDirichlet();

        return bcTypes;
    }
```

The following function specifies the __values on Dirichlet boundaries__ (pressures).
We set 0 Pa at the outlet and a value based on the given pressure gradient
and the length of the domain at the inlet.

```cpp
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
```

#### Upscaling
<details><summary> Click to show auxiliary functions needed for the upscaling process</summary>

```cpp

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
```

</details>

Return the label of inlet pores assuming a previously set direction.

```cpp
    int inletPoreLabel() const
    {
        static constexpr std::array<int, 3> label = {1, 3, 5};
        return label[direction_];
    }
```

Return the label of outlet pores assuming a previously set direction.

```cpp
    int outletPoreLabel() const
    {
        static constexpr std::array<int, 3> label = {2, 4, 6};
        return label[direction_];
    }
```

</details>

```cpp

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
```

</details>

</details>


| [:arrow_left: Back to the main documentation](../README.md) | [:arrow_right: Continue with part 2](main.md) |
|---|---:|

