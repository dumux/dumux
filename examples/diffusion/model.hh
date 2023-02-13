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

#ifndef DUMUX_EXAMPLES_DIFFUSION_MODEL_HH
#define DUMUX_EXAMPLES_DIFFUSION_MODEL_HH

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/volumevariables.hh>

#include <dumux/discretization/method.hh>

namespace Dumux {

namespace Properties::TTag {
//! The diffusion model tag that we can specialize properties for
struct DiffusionModel {};
} // end namespace Properties::TTag

/*!
 * \brief The local residual of the diffusion model
 */
template<class TypeTag>
class DiffusionModelLocalResidual
: public GetPropType<TypeTag, Properties::BaseLocalResidual>
{
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

    /*!
     * \brief Evaluate the rate of change of all conserved quantities
     */
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        NumEqVector storage;
        storage[Indices::massBalanceEqIdx] = volVars.priVar(Indices::concentrationIdx);
        return storage;
    }

    /*!
     * \brief Evaluate the fluxes over a face of a sub control volume
     * Here we evaluate the flow rate, F = -D∇c·n A
     */
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        static_assert(DiscretizationMethods::isCVFE<GridGeometry::DiscretizationMethod>
            "This local residual is hard-coded to control-volume finite element schemes");

        // compute ∇c at the intergration point of this sub control volume face
        const auto& fluxVarCache = elemFluxVarsCache[scvf];
        Dune::FieldVector<Scalar, dimWorld> gradConcentration(0.0);
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];
            gradConcentration.axpy(
                volVars.priVar(Indices::concentrationIdx),
                fluxVarCache.gradN(scv.indexInElement())
            );
        }

        NumEqVector flux;

        // compute the flux as (v^T M v) or -n^T D ∇c A
        flux[Indices::massBalanceEqIdx] = -1.0*vtmv(
            scvf.unitOuterNormal(), problem.diffusionCoefficient(), gradConcentration
        )*scvf.area();

        return flux;
    }
};

} // end namespace Dumux

namespace Dumux::Properties {

//! set the local residual to the one defined above
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::DiffusionModel>
{ using type = DiffusionModelLocalResidual<TypeTag>; };

//! Set the default type of scalar values to double
template<class TypeTag>
struct Scalar<TypeTag, TTag:: DiffusionModel >
{ using type = double; };

//! Set the default primary variable vector to a vector of size of number of equations
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag:: DiffusionModel >
{
    using type = Dune::FieldVector<
        GetPropType<TypeTag, Properties::Scalar>,
        GetPropType<TypeTag, Properties::ModelTraits>::numEq()
    >;
};

//! Set the model traits property
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

//! Set the volume variables property
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

#endif
