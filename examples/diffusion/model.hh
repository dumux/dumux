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

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// VolumeVariables //////// (often in a separate file volumevariables.hh) ////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

namespace Dumux {

template <class Traits>
class DiffusionModelVolumeVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;
    static_assert(Traits::PrimaryVariables::dimension == Traits::ModelTraits::numEq());
public:
    //! export the type used for the primary variables
    using PrimaryVariables = typename Traits::PrimaryVariables;
    //! export the indices type
    using Indices = typename Traits::ModelTraits::Indices;

    /*!
     * \brief Update all quantities for a given control volume
     */
    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    void update(const ElementSolution& elemSol,
                const Problem& problem,
                const Element& element,
                const SubControlVolume& scv)
    {
        priVars_ = elemSol[scv.indexInElement()];
    }

    Scalar concentration() const
    { return priVars_[Indices::concentration0Idx]; }

    Scalar priVar(const int pvIdx) const
    { return priVars_[pvIdx]; }

    const PrimaryVariables& priVars() const
    { return priVars_; }

    // for compatibility with more general models
    Scalar extrusionFactor() const
    { return 1.0; }

private:
    PrimaryVariables priVars_;
};

} // end namespace Dumux

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// LocalResidual //////////// (often in a separate file localresidual.hh) ////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>

namespace Dumux {

template<class TypeTag>
class DiffusionModelLocalResidual
: public GetPropType<TypeTag, Properties::BaseLocalResidual>
{
    // the base local residual is selected depending on the chosen discretization scheme
    using ParentType = GetPropType<TypeTag, Properties::BaseLocalResidual>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using VolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::VolumeVariables;

    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;

    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
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
        storage[Indices::concentration0Idx] = volVars.concentration();
        return storage;
    }

    /*!
     * \brief Evaluate the fluxes over a face of a sub control volume
     * Here we evaluate the flow rate, F = -D∇c·n A
     *
     * TODO: why is this called flux, if we expect it to be integrated already?
     * computeFluxIntegral?
     */
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        const auto& fluxVarCache = elemFluxVarsCache[scvf];
        Dune::FieldVector<Scalar, dimWorld> gradConcentration(0.0);
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];
            gradConcentration.axpy(volVars.concentration(), fluxVarCache.gradN(scv.indexInElement()));
        }

        const auto D = problem.diffusionCoefficient();

        NumEqVector flux;
        flux[Indices::massBalanceEq0Idx] = -1.0*vtmv(scvf.unitOuterNormal(), D, gradConcentration)*scvf.area();

        return flux;
    }
};

} // end namespace Dumux
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Model properties/traits ///////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#include <dumux/common/properties.hh>

namespace Dumux::Properties {

namespace TTag {
struct DiffusionModel {};
} // end namespace TTag

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
    struct Traits
    {
        struct Indices
        {
            static constexpr int concentration0Idx = 0;
            static constexpr int massBalanceEq0Idx = 0;
        };

        static constexpr int numEq() { return 1; }
    };

    using type = Traits;
};

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::DiffusionModel>
{ using type = DiffusionModelLocalResidual<TypeTag>; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::DiffusionModel>
{
    struct Traits
    {
        using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
        using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    };

    using type = DiffusionModelVolumeVariables<Traits>;
};

} // end namespace Dumux::Properties

#endif
