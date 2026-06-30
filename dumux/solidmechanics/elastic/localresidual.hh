// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Elastic
 * \brief Element-wise calculation of the local residual for problems
 *        using the elastic model considering linear elasticity.
 */
#ifndef DUMUX_SOLIDMECHANICS_ELASTIC_LOCAL_RESIDUAL_HH
#define DUMUX_SOLIDMECHANICS_ELASTIC_LOCAL_RESIDUAL_HH

#include <dune/common/fmatrix.hh>
#include <dune/grid/common/intersectioniterator.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/common/typetraits/localdofs_.hh>

namespace Dumux {

/*!
 * \ingroup Elastic
 * \brief Element-wise calculation of the local residual for problems
 *        using the elastic model considering linear elasticity.
 */
template<class TypeTag>
class ElasticLocalResidual
: public DiscretizationDefaultLocalOperator<TypeTag>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = typename ElementVolumeVariables::VolumeVariables;

    // class assembling the stress tensor
    using StressType = GetPropType<TypeTag, Properties::StressType>;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using StressTensor = Dune::FieldMatrix<Scalar, dim, dimWorld>;

    using GV = GetPropType<TypeTag, Properties::GridVariables>;
    using GVC = typename GV::GridVariablesCache;
    using ElementVariables = typename GVC::LocalView;

public:
    using ParentType::ParentType;
    using ElementResidualVector = typename ParentType::ElementResidualVector;

    // =========================================================
    // Old interface (FVAssembler + FVGridVariables)
    // =========================================================

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
        return StressType::force(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
    }

    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume &scv) const
    {
        NumEqVector source(0.0);
        source += problem.source(element, fvGeometry, elemVolVars, scv);
        source += problem.scvPointSources(element, fvGeometry, elemVolVars, scv);

        static const bool gravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");
        if (gravity)
        {
            const auto& g = problem.spatialParams().gravity(scv.center());
            for (int dir = 0; dir < GridView::dimensionworld; ++dir)
                source[Indices::momentum(dir)] += elemVolVars[scv].solidDensity()*g[dir];
        }
        return source;
    }

    // =========================================================
    // New interface (Experimental::Assembler + Experimental::GridVariables)
    // =========================================================

    //! Storage is zero for quasi-static elasticity.
    template<class ElemVars>
    NumEqVector storageIntegral(const FVElementGeometry&,
                                const ElemVars&,
                                const SubControlVolume&, bool) const
    { return NumEqVector(0.0); }

    //! Integrated body-force source (gravity) over a sub-control volume.
    template<class ElemVars>
    NumEqVector sourceIntegral(const FVElementGeometry& fvGeometry,
                               const ElemVars& elemVars,
                               const SubControlVolume& scv) const
    {
        NumEqVector source(0.0);

        static const bool gravity = getParamFromGroup<bool>(
            this->asImp().problem().paramGroup(), "Problem.EnableGravity");
        if (gravity)
        {
            const auto& problem = this->asImp().problem();
            for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scv))
            {
                const auto& g = problem.spatialParams().gravity(qpData.ipData().global());
                for (int dir = 0; dir < dimWorld; ++dir)
                    source[Indices::momentum(dir)] +=
                        elemVars[scv].solidDensity() * g[dir] * qpData.weight();
            }
        }
        return source;
    }

    //! Interior flux integral: -sigma·n integrated over the scvf QPs (CV path, vertex DOFs).
    template<class ElemVars>
    NumEqVector fluxIntegral(const FVElementGeometry& fvGeometry,
                             const ElemVars& elemVars,
                             const SubControlVolumeFace& scvf) const
    {
        const auto& problem = this->asImp().problem();
        NumEqVector f(0.0);

        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scvf))
        {
            const auto& ipCache = cache(elemVars, qpData.ipData());
            const auto& lame = problem.spatialParams().lameParamsAtPos(qpData.ipData().global());
            const auto sigma = this->stressTensor_(fvGeometry, elemVars, ipCache, lame);

            // Residual contribution: -sigma·n·weight  (sign: weak form of -div sigma = 0)
            for (int i = 0; i < dim; ++i)
                for (int j = 0; j < dimWorld; ++j)
                    f[i] -= sigma[i][j] * scvf.unitOuterNormal()[j] * qpData.weight();
        }
        return f;
    }

    //! FE contributions for non-CV (edge/face-midpoint) DOFs — only active for hybrid schemes (PQ2).
    //! Volume part: ∫ sigma : ∇φ_i dX.  Boundary part: Neumann traction on boundary intersections.
    void addToElementFluxAndSourceResidual(ElementResidualVector& residual,
                                           const Problem& problem,
                                           const Element& element,
                                           const FVElementGeometry& fvGeometry,
                                           const ElementVariables& elemVars,
                                           const auto& elemBcTypes) const
    {
        if constexpr (!Dumux::Detail::LocalDofs::hasNonCVLocalDofsInterface<FVElementGeometry>())
            return;

        if (nonCVLocalDofs(fvGeometry).empty())
            return;

        // Volume integral: ∫ sigma : ∇φ_i dX  for each non-CV DOF i
        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, element))
        {
            const auto& ipCache = cache(elemVars, qpData.ipData());
            const auto& lame = problem.spatialParams().lameParamsAtPos(qpData.ipData().global());
            const auto sigma = this->stressTensor_(fvGeometry, elemVars, ipCache, lame);

            for (const auto& localDof : nonCVLocalDofs(fvGeometry))
            {
                const auto idx = localDof.index();
                for (int i = 0; i < dim; ++i)
                    for (int j = 0; j < dimWorld; ++j)
                        residual[idx][i] += sigma[i][j] * ipCache.gradN(idx)[j] * qpData.weight();
            }
        }

        // Boundary integral: Neumann traction for non-CV DOFs on Neumann faces
        if (!elemBcTypes.hasNeumann())
            return;

        for (const auto& intersection : intersections(fvGeometry.gridGeometry().gridView(), element))
        {
            if (!intersection.boundary())
                continue;

            const auto isecBcTypes = problem.boundaryTypesAtPos(intersection.geometry().center());
            if (!isecBcTypes.hasNeumann())
                continue;

            for (const auto& qpData : CVFE::quadratureRule(fvGeometry, intersection))
            {
                const auto& ipCache = cache(elemVars, qpData.ipData());
                const auto flux = problem.boundaryFlux(fvGeometry, elemVars, qpData.ipData());
                const auto& shapeValues = ipCache.shapeValues();

                for (const auto& localDof : nonCVLocalDofs(fvGeometry))
                {
                    const auto idx = localDof.index();
                    for (int eqIdx = 0; eqIdx < dim; ++eqIdx)
                        if (isecBcTypes.isNeumann(eqIdx))
                            residual[idx][eqIdx] +=
                                double(shapeValues[idx]) * flux[eqIdx] * qpData.weight();
                }
            }
        }
    }
private:
    template<class ElemVars, class IpCache, class LameParams>
    StressTensor stressTensor_(const FVElementGeometry& fvGeometry,
                               const ElemVars& elemVars,
                               const IpCache& ipCache,
                               const LameParams& lame) const
    {
        StressTensor gradU(0.0);
        for (int dir = 0; dir < dim; ++dir)
            for (const auto& localDof : localDofs(fvGeometry))
                gradU[dir].axpy(elemVars[localDof].displacement(dir),
                                ipCache.gradN(localDof.index()));

        StressTensor epsilon;
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dimWorld; ++j)
                epsilon[i][j] = 0.5*(gradU[i][j] + gradU[j][i]);

        StressTensor sigma(0.0);
        const auto trEps = trace(epsilon);
        for (int i = 0; i < dim; ++i)
        {
            sigma[i][i] = lame.lambda()*trEps;
            for (int j = 0; j < dimWorld; ++j)
                sigma[i][j] += 2.0*lame.mu()*epsilon[i][j];
        }
        return sigma;
    }
};

} // end namespace Dumux

#endif
