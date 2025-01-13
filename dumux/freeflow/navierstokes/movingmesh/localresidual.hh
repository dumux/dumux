// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesResidualImpl
 */
#ifndef DUMUX_NAVIERSTOKES_CVFE_MOVINGMESH_LOCAL_RESIDUAL_HH
#define DUMUX_NAVIERSTOKES_CVFE_MOVINGMESH_LOCAL_RESIDUAL_HH

#include <dune/common/hybridutilities.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/transpose.hh>

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>

#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/method.hh>
#include <dumux/assembly/cvfelocalresidual.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Element-wise calculation of the Navier-Stokes momentum residual on a moving mesh for models using CVFE discretizations
 */
template<class TypeTag>
class NavierStokesMomentumMovingMeshCVFELocalResidual
: public GetPropType<TypeTag, Properties::BaseLocalResidual>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseLocalResidual>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    using GridFluxVariablesCache = typename GridVariables::GridFluxVariablesCache;
    using ElementFluxVariablesCache = typename GridFluxVariablesCache::LocalView;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Implementation = GetPropType<TypeTag, Properties::LocalResidual>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    using Extrusion = Extrusion_t<GridGeometry>;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    static constexpr auto dim = GridView::dimension;

    using Tensor = Dune::FieldMatrix<Scalar, dim, dim>;

public:
    //! Use the parent type's constructor
    using ParentType::ParentType;

    NumEqVector computeStorage(const Problem& problem,
                               const FVElementGeometry& fvGeometry,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars,
                               const bool isPreviousStorage) const
    {
        const auto J = problem.meshMotion().detF(fvGeometry, scv);
        return problem.density(fvGeometry.element(), fvGeometry, scv, isPreviousStorage) * volVars.velocity() * J;
    }

    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume& scv) const
    {
        NumEqVector source = ParentType::computeSource(problem, element, fvGeometry, elemVolVars, scv);

        // add rho*g (note that gravity might be zero in case it's disabled in the problem)
        source += problem.density(element, fvGeometry, scv) * problem.gravity();

        // integrate in reference coordinates dv = J*dV
        const auto J = problem.meshMotion().detF(fvGeometry, scv);
        source *= J;

        return source;
    }

    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        const auto& fluxVarCache = elemFluxVarsCache[scvf];

        // from moving mesh coupling
        const auto F = problem.meshMotion().F(fvGeometry, scvf);
        auto invF = F; invF.invert();
        const auto J = F.determinant();
        const auto pressure = problem.pressure(element, fvGeometry, scvf);
        const auto referencePressure = problem.referencePressure();

        // interpolate velocity gradient at scvf
        Tensor gradV(0.0);
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];
            for (int dir = 0; dir < dim; ++dir)
                gradV[dir].axpy(volVars.velocity(dir), fluxVarCache.gradN(scv.indexInElement()));
        }

        const auto mu = problem.effectiveViscosity(element, fvGeometry, scvf);

        // compute transformed Cauchy stress σ = -mu (∇v F^-1 + F^-T (∇v)^T) + pI
        auto sigma = -mu*((gradV * invF) + (transpose(invF) * transpose(gradV)));
        for (int i = 0; i < dim; ++i)
            sigma[i][i] += (pressure-referencePressure);

        const NumEqVector flux =
            J*mv(sigma, mv(transpose(invF), scvf.unitOuterNormal())) // σ J F^-T N dA
            * Extrusion::area(fvGeometry, scvf);

        return flux;
    }
};

/*!
 * \ingroup NavierStokesModel
 * \brief Element-wise calculation of the Navier-Stokes mass residual on a moving mesh for models using CVFE discretizations
 */
template<class TypeTag>
class NavierStokesMassMovingMeshCVFELocalResidual
: public GetPropType<TypeTag, Properties::BaseLocalResidual>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseLocalResidual>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    using GridFluxVariablesCache = typename GridVariables::GridFluxVariablesCache;
    using ElementFluxVariablesCache = typename GridFluxVariablesCache::LocalView;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Implementation = GetPropType<TypeTag, Properties::LocalResidual>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    using Extrusion = Extrusion_t<GridGeometry>;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    static constexpr auto dim = GridView::dimension;

public:
    //! Use the parent type's constructor
    using ParentType::ParentType;

    NumEqVector computeStorage(const Problem& problem,
                               const FVElementGeometry& fvGeometry,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars,
                               const bool isPreviousStorage) const
    {
        return NumEqVector(0.0);
    }

    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume& scv) const
    {
        NumEqVector source = ParentType::computeSource(problem, element, fvGeometry, elemVolVars, scv);

        // integrate in reference coordinates dv = J*dV
        const auto J = problem.meshMotion().detF(fvGeometry, scv);
        source *= J;

        return source;
    }

    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        const auto velocity = problem.faceVelocity(element, fvGeometry, scvf);
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto density = insideVolVars.density();

        const auto F = problem.meshMotion().F(fvGeometry, scvf);
        auto invF = F; invF.invert();
        const auto J = F.determinant();

        const Scalar flux = density*(velocity
            * mv(transpose(invF), scvf.unitOuterNormal()))*J // v J F^-T N dA
            * Extrusion::area(fvGeometry, scvf);
        return flux;
    }
};

} // end namespace Dumux

#endif
