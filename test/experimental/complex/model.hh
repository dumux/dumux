// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_EXAMPLES_COMPLEX_HELMHOLTZ_MODEL_HH
#define DUMUX_EXAMPLES_COMPLEX_HELMHOLTZ_MODEL_HH

#include <dune/common/fvector.hh>
#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/volumevariables.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/defaultlocaloperator.hh>

namespace Dumux::Properties::TTag {

struct ComplexHelmholtzModel {};

} // end namespace Dumux::Properties::TTag

namespace Dumux {
template<class TypeTag>
class ComplexHelmholtzModelLocalResidual
: public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;

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
        static_assert(DiscretizationMethods::isCVFE<typename GridGeometry::DiscretizationMethod>,
            "This local residual is hard-coded to control-volume finite element schemes");

        const auto& fluxVarCache = elemFluxVarsCache[scvf];
        Dune::FieldVector<std::complex<double>, dimWorld> gradPhi(0.0);
        for (const auto& localDof : localDofs(fvGeometry))
        {
            const auto& volVars = elemVolVars[localDof.index()];
            gradPhi.axpy(
                volVars.priVar(Indices::phiIdx),
                fluxVarCache.gradN(localDof.index())
            );
        }

        NumEqVector flux;
        flux[Indices::balanceEqIdx] = (gradPhi*scvf.unitOuterNormal())*scvf.area();
        return flux;
    }

    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume& scv) const
    {
        NumEqVector source(0.0);

        // add contribution from possible point sources
        if (!problem.pointSourceMap().empty())
            source += problem.scvPointSources(element, fvGeometry, elemVolVars, scv);

        source += problem.source(element, fvGeometry, elemVolVars, scv);

        // add contribution from the time-derivative term (source term in the frequency domain)
        const auto& volVars = elemVolVars[scv];
        source -= problem.waveNumberSquared()*volVars.priVar(Indices::phiIdx);

        return source;
    }
};
} // end namespace Dumux


namespace Dumux::Properties {

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::ComplexHelmholtzModel>
{ using type = ComplexHelmholtzModelLocalResidual<TypeTag>; };

template<class TypeTag>
struct Scalar<TypeTag, TTag::ComplexHelmholtzModel>
{ using type = double; };

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::ComplexHelmholtzModel>
{
    struct type
    {
        struct Indices
        {
            static constexpr int phiIdx = 0;
            static constexpr int balanceEqIdx = 0;
        };

        static constexpr int numEq() { return 1; }
    };
};

template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::ComplexHelmholtzModel>
{
    using type = Dune::FieldVector<
        std::complex<double>,
        GetPropType<TypeTag, Properties::ModelTraits>::numEq()
    >;
};

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::ComplexHelmholtzModel>
{
    struct Traits
    {
        using PrimaryVariables
            = GetPropType<TypeTag, Properties::PrimaryVariables>;
    };
    using type = BasicVolumeVariables<Traits>;
};

template<class TypeTag>
struct JacobianMatrix<TypeTag, TTag::ComplexHelmholtzModel>
{
private:
    using Scalar = std::complex<double>;
    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    using MatrixBlock = typename Dune::FieldMatrix<Scalar, numEq, numEq>;
public:
    using type = typename Dune::BCRSMatrix<MatrixBlock>;
};

} // end namespace Dumux::Properties

#endif
