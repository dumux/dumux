// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Thick cylindrical shell under pressure — T2 (PQ2 single-field, pure quadratic
// displacement, no pressure split). Locking-free without a mixed formulation
// because P2 has enough displacement modes to handle near-incompressibility.
// Uses the new Experimental::GridVariables + HybridCVFEGridVariablesCache so the
// ElementVariables cover all local DOFs (vertices + edge/face midpoints of P2).
#ifndef DUMUX_HYPERELASTIC_THICK_CYLINDRICAL_SHELL_PROPERTIES_HH
#define DUMUX_HYPERELASTIC_THICK_CYLINDRICAL_SHELL_PROPERTIES_HH

#include <dune/alugrid/grid.hh>

#include <dumux/common/concepts/ipdata_.hh>
#include <dumux/discretization/pq2.hh>
#include <dumux/solidmechanics/hyperelastic/model.hh>

// New variables concept
#include <dumux/discretization/gridvariables.hh>
#include <dumux/discretization/cvfe/hybrid/gridvariablescache.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux {

#include <dune/common/fmatrix.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

/*!
 * \brief PQ2 hyperelastic momentum local residual using the new BoundaryFace interface.
 *        Assembles −∇_X·P = 0 with P from problem.firstPiolaKirchhoffStressTensor(F).
 * CV part  (vertex DOFs):    fluxIntegral via scvf loop.
 * FE part  (edge/face DOFs): addToElementFluxAndSourceResidual via element QPs.
 * Modeled on Dumux::PoroElasticLargeDefMomentumLocalResidual (without coupling).
 */
template<class TypeTag>
class HyperelasticPQ2LocalResidual : public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    static constexpr int dim = GridView::dimension;
    using Tensor = Dune::FieldMatrix<Scalar, dim, dim>;

    using GV = GetPropType<TypeTag, Properties::GridVariables>;
    using GVC = typename GV::GridVariablesCache;
    using ElementVariables = typename GVC::LocalView;

public:
    using ParentType::ParentType;
    using ElementResidualVector = typename ParentType::ElementResidualVector;

    NumEqVector computeStorage(const auto&, const auto&, const auto&) const
    { return NumEqVector(0.0); }
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

    //! FE contribution for non-CV (edge/face midpoint) DOFs via element QPs.
    //! Only the volume integral ∫ P : ∇φ_i dX; boundary traction is handled
    //! by addBoundaryFluxIntegral via problem.addFEBoundaryFluxIntegral.
    void addToElementFluxAndSourceResidual(ElementResidualVector& residual,
                                           const auto& problem,
                                           const auto& element,
                                           const FVElementGeometry& fvGeometry,
                                           const ElementVariables& elemVars) const
    {
        if constexpr (Dumux::Detail::LocalDofs::hasNonCVLocalDofsInterface<FVElementGeometry>())
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

//! Minimal momentum volume variables: store the displacement at each local DOF.
//! Reading elemSol[ipData.localDofIndex()] works uniformly for vertex (CV) and
//! edge/face-midpoint (non-CV) DOFs of the PQ2 hybrid discretization.
template<class Traits>
class ThickShellMomentumVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;
public:
    using PrimaryVariables = typename Traits::PrimaryVariables;
    using Indices = typename Traits::ModelTraits::Indices;

    template<class ElementSolution, class Problem, class FVElementGeometry, Concept::LocalDofIpData IpData>
    void update(const ElementSolution& elemSol, const Problem&,
                const FVElementGeometry&, const IpData& ipData)
    { priVars_ = elemSol[ipData.localDofIndex()]; }

    PrimaryVariables displacement() const { return priVars_; }
    Scalar displacement(int dirIdx) const { return priVars_[dirIdx]; }
    Scalar priVar(int pvIdx) const { return priVars_[pvIdx]; }
    const PrimaryVariables& priVars() const { return priVars_; }
    Scalar extrusionFactor() const { return 1.0; }

private:
    PrimaryVariables priVars_;
};
} // end namespace Dumux

namespace Dumux::Properties {

namespace TTag {
struct ThickCylindricalShell
{
    using InheritsFrom = std::tuple<Hyperelastic, PQ2HybridModel>;
    using Grid = Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming>;
};
} // end namespace TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::ThickCylindricalShell>
{ using type = ThickCylindricalShellProblem<TypeTag>; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::ThickCylindricalShell>
{ using type = HyperelasticPQ2LocalResidual<TypeTag>; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::ThickCylindricalShell>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using Traits = HyperelasticVolumeVariablesTraits<PV, MT>;
public:
    using type = ThickShellMomentumVariables<Traits>;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::ThickCylindricalShell>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using type = ThickCylindricalShellSpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::ThickCylindricalShell>
{ static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::ThickCylindricalShell>
{ static constexpr bool value = true; };

//! PQ2 uses HybridCVFEGridVariablesCache (stores element-level QP data).
template<class TypeTag>
struct GridVariables<TypeTag, TTag::ThickCylindricalShell>
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
