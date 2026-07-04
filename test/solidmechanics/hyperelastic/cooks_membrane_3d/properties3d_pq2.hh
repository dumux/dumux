// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
// 3D Cook's membrane — T2 (PQ2 single-field, pure quadratic displacement, no pressure split).
// Analogue of Table 2.5 "T2" in SPP-1748 PDF.
// Uses the new Experimental::GridVariables + HybridCVFEGridVariablesCache so that
// ElementVariables covers all local DOFs (vertices + edge/face midpoints of P2).
#ifndef DUMUX_COOKS_MEMBRANE_3D_PQ2_PROPERTIES_HH
#define DUMUX_COOKS_MEMBRANE_3D_PQ2_PROPERTIES_HH

#include <dune/alugrid/grid.hh>

#include <dumux/discretization/pq2.hh>
#include <dumux/solidmechanics/hyperelastic/model.hh>
#include <dumux/solidmechanics/hyperelastic/volumevariables.hh>

// New variables concept
#include <dumux/discretization/gridvariables.hh>
#include <dumux/discretization/cvfe/hybrid/gridvariablescache.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>

#include "spatialparams.hh"
#include "problem3d.hh"

namespace Dumux {

// Forward declarations needed for the local residual below
#include <dune/common/fmatrix.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

/*!
 * \brief PQ2 hyperelastic local residual — uses ALL local DOFs for F computation.
 * CV part  (vertex DOFs):    fluxIntegral via scvf loop.
 * FE part  (edge/face DOFs): addToElementFluxAndSourceResidual via element QPs.
 */
template<class TypeTag>
class HyperelasticPQ2LocalResidual : public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    static constexpr int dim = GridView::dimension;
    using Tensor = Dune::FieldMatrix<Scalar, dim, dim>;

    using GV = GetPropType<TypeTag, Properties::GridVariables>;
    using GVC = typename GV::GridVariablesCache;
    using ElementVariables = typename GVC::LocalView;

public:
    using ParentType::ParentType;
    using ElementResidualVector = typename ParentType::ElementResidualVector;

    NumEqVector storageIntegral(const FVElementGeometry&, const ElementVariables&, const auto&, bool) const
    { return NumEqVector(0.0); }
    NumEqVector sourceIntegral(const FVElementGeometry&, const ElementVariables&, const auto&) const
    { return NumEqVector(0.0); }

    //! CV flux integral for vertex DOFs: -P·n integrated over scvf QPs.
    NumEqVector fluxIntegral(const FVElementGeometry& fvGeometry,
                             const ElementVariables& elemVars,
                             const SubControlVolumeFace& scvf) const
    {
        const auto& n = scvf.unitOuterNormal();
        const auto& problem = this->asImp().problem();
        NumEqVector f(0.0);
        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scvf))
        {
            const auto& ipCache = cache(elemVars, qpData.ipData());
            Tensor F(0.0);
            for (const auto& localDof : localDofs(fvGeometry))
                for (int d = 0; d < dim; ++d)
                    F[d].axpy(elemVars[localDof].displacement(d), ipCache.gradN(localDof.index()));
            for (int d = 0; d < dim; ++d) F[d][d] += 1.0;
            const auto P = problem.firstPiolaKirchhoffStressTensor(F);
            for (int i = 0; i < dim; ++i)
                for (int j = 0; j < dim; ++j)
                    f[i] -= P[i][j] * n[j] * qpData.weight();
        }
        return f;
    }

    //! FE contribution for the non-CV (edge/face midpoint) DOFs: the standard
    //! Galerkin volume integral int P : grad(phi_i) dX. Signature matches the
    //! current Experimental::LocalResidual base (5 args); boundary tractions for
    //! all dofs are added by the assembler base via the problem's boundaryFlux().
    void addToElementFluxAndSourceResidual(ElementResidualVector& residual,
                                           const Problem& problem,
                                           const Element& element,
                                           const FVElementGeometry& fvGeometry,
                                           const ElementVariables& elemVars) const
    {
        if constexpr (!Dumux::Detail::LocalDofs::hasNonCVLocalDofsInterface<FVElementGeometry>())
            return;
        else
        {
            if (nonCVLocalDofs(fvGeometry).empty())
                return;

            for (const auto& qpData : CVFE::quadratureRule(fvGeometry, element))
            {
                const auto& ipCache = cache(elemVars, qpData.ipData());
                Tensor F(0.0);
                for (const auto& localDof : localDofs(fvGeometry))
                    for (int d = 0; d < dim; ++d)
                        F[d].axpy(elemVars[localDof].displacement(d), ipCache.gradN(localDof.index()));
                for (int d = 0; d < dim; ++d) F[d][d] += 1.0;
                const auto P = problem.firstPiolaKirchhoffStressTensor(F);
                for (const auto& nonCVdof : nonCVLocalDofs(fvGeometry))
                {
                    const auto idx = nonCVdof.index();
                    for (int i = 0; i < dim; ++i)
                        for (int j = 0; j < dim; ++j)
                            residual[idx][i] += P[i][j] * ipCache.gradN(idx)[j] * qpData.weight();
                }
            }
        }
    }
};

//! Wrapper around HyperelasticVolumeVariables adding the new-style update()
//! needed by HybridCVFEGridVariablesCache (covers non-SCV edge/face midpoint DOFs).
template<class Traits>
class HyperelasticVolumeVariablesWithNewUpdate : public HyperelasticVolumeVariables<Traits>
{
public:
    using HyperelasticVolumeVariables<Traits>::HyperelasticVolumeVariables;

    template<class ES, class P, class FVG, class IpD>
    void update(const ES& elemSol, const P& problem, const FVG& fvGeometry, const IpD& ipData)
    {
        const auto idx = ipData.localDofIndex();
        if (idx < fvGeometry.numScv())
        {
            BasicVolumeVariables<Traits>::update(elemSol, problem, fvGeometry.element(), fvGeometry.scv(idx));
        }
        else
        {
            // Non-SCV DOF (P2 edge/face midpoint): proxy remaps index to scv(0)
            const auto& scv0 = fvGeometry.scv(0);
            struct Proxy {
                const ES& s; std::size_t from, to;
                auto operator[](std::size_t i) const { return (i == from) ? s[to] : s[i]; }
            } proxy{elemSol, scv0.localDofIndex(), idx};
            BasicVolumeVariables<Traits>::update(proxy, problem, fvGeometry.element(), scv0);
        }
    }
};
} // end namespace Dumux

namespace Dumux::Properties {

namespace TTag {
struct CooksMembrane3D_PQ2
{
    using InheritsFrom = std::tuple<Hyperelastic, PQ2HybridModel>;
    using Grid = Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming>;
};
} // end namespace TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::CooksMembrane3D_PQ2>
{ using type = CooksMembrane3DProblem<TypeTag>; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::CooksMembrane3D_PQ2>
{ using type = HyperelasticPQ2LocalResidual<TypeTag>; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::CooksMembrane3D_PQ2>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using Traits = HyperelasticVolumeVariablesTraits<PV, MT>;
public:
    using type = HyperelasticVolumeVariablesWithNewUpdate<Traits>;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::CooksMembrane3D_PQ2>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using type = CooksMembraneHyperelasticSpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::CooksMembrane3D_PQ2>
{ static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::CooksMembrane3D_PQ2>
{ static constexpr bool value = true; };

//! PQ2 uses HybridCVFEGridVariablesCache (stores element-level QP data).
template<class TypeTag>
struct GridVariables<TypeTag, TTag::CooksMembrane3D_PQ2>
{
private:
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Variables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using IPDataCache = Dumux::CVFE::LocalBasisInterpolationPointData<GG>;
    using Traits = Dumux::Experimental::CVFE::HybridCVFEDefaultGridVariablesCacheTraits<Problem, Variables, IPDataCache>;
    using GVC = Dumux::Experimental::CVFE::HybridCVFEGridVariablesCache<Traits, enableCache>;
public:
    using type = Dumux::Experimental::GridVariables<GG, GVC>;
};

} // end namespace Dumux::Properties
#endif
