// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_EXAMPLES_DIFFUSION_MODEL_HH
#define DUMUX_EXAMPLES_DIFFUSION_MODEL_HH

// In the file `model.hh`, we define the model equations and
// set all default model properties. The setup consist of three steps:
// 1. Create a model type tag (used to specialize properties)
// 2. Define the local residual class implementing the discrete equation
// 3. Specialize important properties of the model such that Dumux knows how to assemble the system matrix
//
// __Table of contents__
//
// [TOC]
//
// We start in `model.hh` with the necessary header includes:
// [[details]] includes
#include <dune/common/fvector.hh>
#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/volumevariables.hh>
#include <dumux/discretization/method.hh>
// [[/details]]
//
// ## 1. Property Tag
//
// The property tag is simply an empty struct with the name `DiffusionModel`
//
// [[content]]
// [[codeblock]]
namespace Dumux::Properties::TTag {
//! The diffusion model tag that we can specialize properties for
struct DiffusionModel {};
} // end namespace Dumux::Properties::TTag
// [[/codeblock]]
// [[/content]]
//
// ## 2. The local (element-wise) residual
//
// The local residual assembles the contribution to the residual for
// all degrees of freedom associated with an element. Here, we use the
// Box method which is based on $P_1$ basis functions (piece-wise linears)
// and the degrees of freedom are on the nodes. Each node is associate with
// exactly one sub control volume (`scv`) per element and several ($2$ in $\mathbb{R}^2$)
// sub control volume faces (`scvf`). In the local residual, we can implement the
// contribution for one `scv` (storage and source terms) or one `scvf` (flux terms).
//
// Let's have a look at the class implementation.
//
// [[content]]
//
// The class `DiffusionModelLocalResidual` inherits from something called `BaseLocalResidual`.
// This base class differs depending on the chosen discretization scheme. For the box method
// (which is a control-volume finite element scheme) used in this example, the property
// `BaseLocalResidual` is specialized to `CVFELocalResidual<TypeTag>`
// in [dumux/discretization/box.hh](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/dumux/discretization/box.hh).
// Since this local residual only works with control-volume finite element schemes due to
// the flux implementation, we could have also chosen to inherit from `public CVFELocalResidual<TypeTag>`.
namespace Dumux {
template<class TypeTag>
class DiffusionModelLocalResidual
: public GetPropType<TypeTag, Properties::BaseLocalResidual>
{
    // [[exclude]]
    // the base local residual is selected depending on the chosen discretization scheme
    using ParentType = GetPropType<TypeTag, Properties::BaseLocalResidual>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using VolumeVariables = typename GridVariables::GridVolumeVariables::VolumeVariables;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    static constexpr int dimWorld = GridView::dimensionworld;

public:
    using ParentType::ParentType;
    // [[/exclude]]
    //
    // **Storage term:** Evaluate the rate of change of all conserved quantities
    // [[codeblock]]
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        NumEqVector storage;
        storage[Indices::massBalanceEqIdx] = volVars.priVar(Indices::concentrationIdx);
        return storage;
    }
    // [[/codeblock]]

    // **Flux term:** Evaluate the fluxes over a face of a sub control volume.
    // Here we evaluate the (integrated) flux
    //
    // ```math
    // F_{K,\sigma} = -D \sum_{B \in \mathcal{B}_K} c_{h,B} \nabla N_B \cdot\boldsymbol{n} \vert \sigma \vert
    // ````
    //
    // [[codeblock]]
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        static_assert(DiscretizationMethods::isCVFE<typename GridGeometry::DiscretizationMethod>,
            "This local residual is hard-coded to control-volume finite element schemes");

        // Compute ∇c at the integration point of this sub control volume face.
        const auto& fluxVarCache = elemFluxVarsCache[scvf];
        Dune::FieldVector<Scalar, dimWorld> gradConcentration(0.0);
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];
            // v.axpy(a, w) means v += a*w
            gradConcentration.axpy(
                volVars.priVar(Indices::concentrationIdx),
                fluxVarCache.gradN(scv.indexInElement())
            );
        }

        NumEqVector flux;

        // Compute the flux with `vtmv` (vector transposed times matrix times vector) or -n^T D ∇c A.
        // The diffusion coefficient comes from the `problem` (see Part 2 of the example).
        flux[Indices::massBalanceEqIdx] = -1.0*vtmv(
            scvf.unitOuterNormal(), problem.diffusionCoefficient(), gradConcentration
        )*scvf.area();

        return flux;
    }
};
} // end namespace Dumux
// [[/codeblock]]
// [[/content]]
//
// ## 3. The model properties
//
// By specializing properties for our type tag `DiffusionModel`,
// every other class that knows about the type tag (this will be
// for example the assembler or the problem), can extract the
// type information that we specify here.
//
// Note that these types can be overwritten for specific problem
// definitions if this is needed (we will show this on the next page).
//
// [[content]]
namespace Dumux::Properties {

// The type of the local residual is the class defined above
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::DiffusionModel>
{ using type = DiffusionModelLocalResidual<TypeTag>; };

// The default scalar type is double
// we compute with double precision floating point numbers
template<class TypeTag>
struct Scalar<TypeTag, TTag::DiffusionModel>
{ using type = double; };

// The model traits specify some information about our equation system.
// Here we have just one equation. We still specify indices so in the
// places where we access primary variables, we can do so with a named variable.
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::DiffusionModel>
{
    struct type
    {
        struct Indices
        {
            static constexpr int concentrationIdx = 0;
            static constexpr int massBalanceEqIdx = 0;
        };

        static constexpr int numEq() { return 1; }
    };
};

// The primary variable vector has entries of type `Scalar` and is
// as large as the number of equations (here 1) but we keep it general.
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::DiffusionModel>
{
    using type = Dune::FieldVector<
        GetPropType<TypeTag, Properties::Scalar>,
        GetPropType<TypeTag, Properties::ModelTraits>::numEq()
    >;
};

// The `BasicVolumeVariables` are the simplest class of volume variables.
// They only store one instance of `PrimaryVariables` for the
// degree of freedom (here: vertex dof) that they are attached to.
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::DiffusionModel>
{
    struct Traits
    {
        using PrimaryVariables
            = GetPropType<TypeTag, Properties::PrimaryVariables>;
    };
    using type = BasicVolumeVariables<Traits>;
};

} // end namespace Dumux::Properties
// [[/content]]
// [[exclude]]
#endif
// [[/exclude]]
