// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief A complex-valued Helmholtz model to test the infrastructure for complex-valued unknowns
 *
 * Solves \f$ -\nabla\cdot(\nabla u) - k^2 u = f \f$ with complex \f$ u, k^2, f \f$
 * on control-volume finite element and cell-centered TPFA discretizations.
 */
#ifndef DUMUX_TEST_COMPLEX_HELMHOLTZ_MODEL_HH
#define DUMUX_TEST_COMPLEX_HELMHOLTZ_MODEL_HH

#include <complex>

#include <dune/common/fvector.hh>
#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/volumevariables.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>

namespace Dumux {

template<class TypeTag>
class ComplexHelmholtzModelLocalResidual
: public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using PrimaryVariable = typename PrimaryVariables::value_type;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

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
    using Extrusion = Extrusion_t<GridGeometry>;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    static constexpr int dimWorld = GridView::dimensionworld;

public:
    using ParentType::ParentType;

    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        return NumEqVector(0.0);
    }

    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        NumEqVector flux(0.0);

        if constexpr (DiscretizationMethods::isCVFE<typename GridGeometry::DiscretizationMethod>)
        {
            const auto& fluxVarCache = elemFluxVarsCache[scvf];
            Dune::FieldVector<PrimaryVariable, dimWorld> gradU(0.0);
            for (const auto& localDof : localDofs(fvGeometry))
                gradU.axpy(elemVolVars[localDof.index()].priVar(Indices::uIdx), fluxVarCache.gradN(localDof.index()));

            flux[Indices::balanceEqIdx] = -(gradU*scvf.unitOuterNormal())*Extrusion::area(fvGeometry, scvf);
        }
        else if constexpr (GridGeometry::discMethod == DiscretizationMethods::cctpfa)
        {
            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
            const Scalar ti = computeTpfaTransmissibility(fvGeometry, scvf, insideScv, 1.0, 1.0);
            Scalar tij = Extrusion::area(fvGeometry, scvf)*ti;
            if (!scvf.boundary())
            {
                const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
                const Scalar tj = -1.0*computeTpfaTransmissibility(fvGeometry, scvf, outsideScv, 1.0, 1.0);
                tij = Extrusion::area(fvGeometry, scvf)*(ti*tj)/(ti + tj);
            }

            const auto uInside = elemVolVars[scvf.insideScvIdx()].priVar(Indices::uIdx);
            const auto uOutside = elemVolVars[scvf.outsideScvIdx()].priVar(Indices::uIdx);
            flux[Indices::balanceEqIdx] = tij*(uInside - uOutside);
        }
        else
            static_assert(Dune::AlwaysFalse<TypeTag>::value, "Discretization method not supported");

        return flux;
    }

    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume& scv) const
    {
        NumEqVector source(0.0);

        if (!problem.pointSourceMap().empty())
            source += problem.scvPointSources(element, fvGeometry, elemVolVars, scv);

        source += problem.source(element, fvGeometry, elemVolVars, scv);

        // the time-derivative term appears as a reaction term in the frequency domain
        const auto& volVars = elemVolVars[scv];
        source[Indices::balanceEqIdx] += problem.waveNumberSquared()*volVars.priVar(Indices::uIdx);

        return source;
    }
};

struct ComplexHelmholtzModelTraits
{
    struct Indices
    {
        static constexpr int uIdx = 0;
        static constexpr int balanceEqIdx = 0;
    };

    static constexpr int numEq() { return 1; }
};

} // end namespace Dumux

namespace Dumux::Properties::TTag {

struct ComplexHelmholtzModel { using InheritsFrom = std::tuple<ModelProperties>; };

} // end namespace Dumux::Properties::TTag

namespace Dumux::Properties {

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::ComplexHelmholtzModel>
{ using type = ComplexHelmholtzModelLocalResidual<TypeTag>; };

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::ComplexHelmholtzModel>
{ using type = ComplexHelmholtzModelTraits; };

//! complex-valued unknowns, the scalar type stays real
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::ComplexHelmholtzModel>
{
    using type = Dune::FieldVector<
        std::complex<GetPropType<TypeTag, Properties::Scalar>>,
        GetPropType<TypeTag, Properties::ModelTraits>::numEq()
    >;
};

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::ComplexHelmholtzModel>
{
    struct Traits
    {
        using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    };
    using type = BasicVolumeVariables<Traits>;
};

} // end namespace Dumux::Properties

#endif
