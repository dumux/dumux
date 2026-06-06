// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoroElasticLargeDef
 * \brief Local residuals for the large-deformation poroelastic mixed u–p_s–p_f model.
 *
 * Three classes are provided:
 *
 * - `PoroElasticLargeDefMomentumLocalResidual`:
 *   Quasi-static momentum balance −∇_X·**P**_eff = **0**.
 *   Identical to `HyperelasticVolIsoMomentumLocalResidual`; the effective stress
 *   tensor (including fluid-pressure Biot coupling) is assembled by the problem.
 *
 * - `PoroElasticLargeDefSolidPressureLocalResidual`:
 *   Algebraic solid-pressure constraint \f$p_s^{\mathrm{eq}}(J) - p_s = 0\f$.
 *   Identical to `HyperelasticVolIsoPressureLocalResidual`.
 *
 * - `PoroElasticLargeDefFluidPressureLocalResidual`:
 *   Transient fluid mass balance (pulled back to reference frame):
 *   \f[ \frac{\partial}{\partial t}\!\left(J + \frac{S_p}{\alpha_B}p_f\right)
 *       - \nabla_X \cdot \!\left(J\,\frac{\kappa}{\mu_f}\,\mathbf{C}^{-1}\nabla_X p_f\right) = 0. \f]
 */
#ifndef DUMUX_POROMECHANICS_POROELASTIC_LARGE_DEFORMATIONS_LOCAL_RESIDUAL_HH
#define DUMUX_POROMECHANICS_POROELASTIC_LARGE_DEFORMATIONS_LOCAL_RESIDUAL_HH

#include <dune/common/fmatrix.hh>

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

namespace Dumux {

/*!
 * \ingroup PoroElasticLargeDef
 * \brief Momentum local residual for the poroelastic large-deformation model.
 *
 * Assembles −∇_X·**P**_eff = **0** where the effective first Piola–Kirchhoff stress
 * **P**_eff is provided by `problem.firstPiolaKirchhoffStressTensor(F, fvGeometry, ip)`.
 * The problem implementation is responsible for incorporating the solid-pressure
 * and fluid-pressure contributions via the coupling manager.
 *
 * This residual is quasi-static (no storage term).
 */
template<class TypeTag>
class PoroElasticLargeDefMomentumLocalResidual
    : public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Tensor = Dune::FieldMatrix<Scalar, GridView::dimension, GridView::dimension>;
    static constexpr int dim = GridView::dimension;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVariables = typename GridVariables::GridVariablesCache::LocalView;

public:
    using ParentType::ParentType;
    using ElementResidualVector = typename ParentType::ElementResidualVector;

    NumEqVector computeStorage(const auto&, const auto&, const auto&) const
    { return NumEqVector(0.0); }

    NumEqVector storageIntegral(const FVElementGeometry&, const ElementVariables&,
                                const auto&, bool) const
    { return NumEqVector(0.0); }

    NumEqVector sourceIntegral(const FVElementGeometry&, const ElementVariables&,
                               const auto&) const
    { return NumEqVector(0.0); }

    //! FE contributions for non-CV DOFs (e.g. PQ1Bubble bubble function)
    void addToElementFluxAndSourceResidual(ElementResidualVector& residual,
                                           const auto& problem,
                                           const auto& element,
                                           const FVElementGeometry& fvGeometry,
                                           const ElementVariables& elemVars,
                                           const auto& elemBcTypes) const
    {
        if constexpr (Dumux::Detail::LocalDofs::hasNonCVLocalDofsInterface<FVElementGeometry>())
        {

        if (nonCVLocalDofs(fvGeometry).empty())
            return;

        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, element))
        {
            const auto& ipCache = cache(elemVars, qpData.ipData());
            const auto& ipGlobal = qpData.ipData().global();

            Tensor F(0.0);
            for (const auto& localDof : localDofs(fvGeometry))
                for (int d = 0; d < dim; ++d)
                    F[d].axpy(elemVars[localDof].displacement(d), ipCache.gradN(localDof.index()));
            for (int d = 0; d < dim; ++d) F[d][d] += 1.0;

            const auto P = problem.firstPiolaKirchhoffStressTensor(F, fvGeometry, ipGlobal);

            for (const auto& nonCVdof : nonCVLocalDofs(fvGeometry))
            {
                const auto idx = nonCVdof.index();
                for (int i = 0; i < dim; ++i)
                    for (int j = 0; j < dim; ++j)
                        residual[idx][i] += P[i][j] * ipCache.gradN(idx)[j] * qpData.weight();
            }
        }

        if (!elemBcTypes.hasFluxBoundary())
            return;

        // Neumann (traction) boundary contribution for the non-CV (FEM) local dofs.
        // Iterate the CVFE boundary faces (works for both PQ1Bubble and PQ2 hybrid),
        // see dumux/freeflow/navierstokes/momentum/cvfe for the same pattern.
        for (const auto& bf : boundaryFaces(fvGeometry))
        {
            const auto& bfBcTypes = elemBcTypes.get(fvGeometry, bf);
            if (!bfBcTypes.hasFluxBoundary()) continue;

            for (const auto& qpData : CVFE::quadratureRule(fvGeometry, bf))
            {
                const auto& ipCache = cache(elemVars, qpData.ipData());
                const auto flux = problem.boundaryFlux(fvGeometry, elemVars, qpData.ipData());
                const auto& shapeValues = ipCache.shapeValues();

                for (const auto& nonCVdof : nonCVLocalDofs(fvGeometry))
                {
                    const auto idx = nonCVdof.index();
                    for (int eqIdx = 0; eqIdx < dim; ++eqIdx)
                        if (bfBcTypes.isFluxBoundary(eqIdx))
                            residual[idx][eqIdx] += double(shapeValues[idx]) * flux[eqIdx] * qpData.weight();
                }
            }
        }
        } // end if constexpr
    }

    //! Interior flux integral: −∫ **P**_eff · **n** dA
    NumEqVector fluxIntegral(const FVElementGeometry& fvGeometry,
                             const ElementVariables& elemVars,
                             const SubControlVolumeFace& scvf) const
    {
        const auto& problem = this->asImp().problem();
        const auto& n = scvf.unitOuterNormal();
        NumEqVector f(0.0);

        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scvf))
        {
            const auto& ipCache = cache(elemVars, qpData.ipData());
            const auto& ipGlobal = qpData.ipData().global();

            Tensor F(0.0);
            for (const auto& localDof : localDofs(fvGeometry))
                for (int d = 0; d < dim; ++d)
                    F[d].axpy(elemVars[localDof].displacement(d), ipCache.gradN(localDof.index()));
            for (int d = 0; d < dim; ++d) F[d][d] += 1.0;

            const auto P = problem.firstPiolaKirchhoffStressTensor(F, fvGeometry, ipGlobal);

            for (int i = 0; i < dim; ++i)
                for (int j = 0; j < dim; ++j)
                    f[Indices::momentum(i)] -= P[i][j] * n[j] * qpData.weight();
        }
        return f;
    }
};

/*!
 * \ingroup PoroElasticLargeDef
 * \brief Solid pressure local residual: algebraic constraint \f$p_s^{\mathrm{eq}}(J) - p_s = 0\f$.
 *
 * The equilibrium solid pressure \f$p_s^{\mathrm{eq}}(J)\f$ is provided by the problem
 * via `problem.volumetricPressure(J)`.  This is identical to
 * `HyperelasticVolIsoPressureLocalResidual` and assembles a purely algebraic
 * (no time derivative) per-SCV constraint.
 */
template<class TypeTag>
class PoroElasticLargeDefSolidPressureLocalResidual
    : public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Tensor = Dune::FieldMatrix<Scalar, GridView::dimension, GridView::dimension>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVariables = typename GridVariables::GridVariablesCache::LocalView;

public:
    using ParentType::ParentType;
    using ElementResidualVector = typename ParentType::ElementResidualVector;

    NumEqVector computeStorage(const auto&, const auto&, const auto&) const
    { return NumEqVector(0.0); }

    NumEqVector storageIntegral(const FVElementGeometry&, const ElementVariables&,
                                const auto&, bool) const
    { return NumEqVector(0.0); }

    NumEqVector fluxIntegral(const FVElementGeometry&, const ElementVariables&,
                             const auto&) const
    { return NumEqVector(0.0); }

    //! Enforces \f$p_s^{\mathrm{eq}}(J) - p_s = 0\f$ averaged over each SCV.
    template<class SubControlVolume_>
    NumEqVector sourceIntegral(const FVElementGeometry& fvGeometry,
                               const ElementVariables& elemVars,
                               const SubControlVolume_& scv) const
    {
        NumEqVector source(0.0);

        const auto& problem = this->asImp().problem();
        const auto& localBasis = fvGeometry.gridGeometry().feCache().get(fvGeometry.element().type()).localBasis();
        // Scratch buffer for the P1 shape values; thread_local so it keeps its
        // capacity across the many calls during numeric-differentiation assembly
        // (evaluateFunction refills it each call) instead of heap-allocating anew.
        // Safe under coloured multithreaded assembly (one element per thread at a time).
        thread_local std::vector< Dune::FieldVector<Scalar, 1> > shapeValues;

        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scv, QuadratureRules::DuneQuadrature<3>{}))
        {
            const Tensor F = problem.couplingManager().deformationGradientAtPoint(
                fvGeometry, qpData.ipData().global());
            const Scalar J = F.determinant();
            const Scalar psEq = problem.volumetricPressure(J);

            localBasis.evaluateFunction(qpData.ipData().local(), shapeValues);
            Scalar ps = 0.0;
            for (const auto& localDof : localDofs(fvGeometry))
                ps += elemVars[localDof].solidPressure() * shapeValues[localDof.index()][0];

            source[0] -= (psEq - ps) * qpData.weight();
        }

        return source;
    }
};

/*!
 * \ingroup PoroElasticLargeDef
 * \brief Fluid pressure local residual for the poroelastic large-deformation model.
 *
 * Assembles the Lagrangian fluid mass balance:
 * \f[ \frac{\partial}{\partial t}\!\left(J + \frac{S_p}{\alpha_B}p_f\right)
 *     - \nabla_X \cdot \!\left(J\,\frac{\kappa}{\mu_f}\,\mathbf{C}^{-1}\nabla_X p_f\right) = 0 \f]
 *
 * **Storage term** (`storageIntegral`):
 * \f$ J + \frac{S_p}{\alpha_B}\,p_f \f$ (for incompressible solid/fluid: just \f$J\f$).
 * \f$J\f$ is obtained from the coupling manager which always uses its current `curSol`.
 * The `isPreviousTimeLevel` flag is intentionally ignored — the multi-stage assembler
 * keeps `curSol` up-to-date for the correct stage via `updateGridVariables`.
 *
 * **Flux term** (`fluxIntegral`):
 * Darcy flux pulled back to the reference frame:
 * \f$ J\,\frac{\kappa}{\mu_f}\,\mathbf{C}^{-1}\nabla_X p_f \cdot \mathbf{n} \f$
 * where \f$\mathbf{C}^{-1} = \mathbf{F}^{-1}\mathbf{F}^{-\top}\f$.
 *
 * Permeability \f$\kappa\f$ and viscosity \f$\mu_f\f$ are provided by
 * `problem.spatialParams()`.
 */
template<class TypeTag>
class PoroElasticLargeDefFluidPressureLocalResidual
    : public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using GlobalPosition = typename GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
    using Tensor = Dune::FieldMatrix<Scalar, GridView::dimension, GridView::dimension>;
    static constexpr int dim = GridView::dimension;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVariables = typename GridVariables::GridVariablesCache::LocalView;

public:
    using ParentType::ParentType;
    using ElementResidualVector = typename ParentType::ElementResidualVector;

    NumEqVector computeStorage(const auto&, const auto&, const auto&) const
    { return NumEqVector(0.0); }

    /*!
     * \brief Storage term: \f$J + S_p/\alpha_B \cdot p_f\f$
     *
     * \f$J\f$ is retrieved from the coupling manager using its current solution.
     *
     * The `isPreviousTimeLevel` flag is intentionally ignored here: the
     * multi-stage assembler ensures that the coupling manager's `curSol` is
     * already set to the correct time level (old or new) for each stage via
     * `MultiStageMultiDomainAssembler::updateGridVariables`. Using the flag
     * would be redundant for multi-stage and incorrect for Jacobian assembly
     * (where the coupling manager must track the perturbed Newton iterate).
     */
    template<class SubControlVolume_>
    NumEqVector storageIntegral(const FVElementGeometry& fvGeometry,
                                const ElementVariables& elemVars,
                                const SubControlVolume_& scv,
                                bool /*isPreviousTimeLevel*/) const
    {
        NumEqVector storage(0.0);
        const auto& problem = this->asImp().problem();

        // Evaluate J at the SCV center (deformation gradient from current displacement).
        // For P1 (Box) elements, the gradient of displacement is constant per element,
        // so any evaluation point in the element gives the same F.
        const Tensor F = problem.couplingManager().deformationGradientAtPoint(
            fvGeometry, scv.center());
        const Scalar J = F.determinant();

        // Fluid pressure at this SCV DOF (direct value, no interpolation needed for P1)
        const Scalar pf = elemVars[scv].fluidPressure();

        // storage = J + Sp/alphaB * pf  (for Sp=0: just J)
        const Scalar alphaB = problem.spatialParams().biotCoefficient(fvGeometry.element(), scv);
        const Scalar Sp = problem.spatialParams().storageCoefficient(fvGeometry.element(), scv);

        // Multiply by SCV volume so the residual is volume-integrated,
        // consistent with fluxIntegral (area-integrated) and sourceIntegral (volume-integrated).
        const Scalar volume = scv.volume();
        storage[0] = (J + (alphaB > 0 ? Sp / alphaB * pf : Scalar(0))) * volume;

        return storage;
    }

    /*!
     * \brief Flux term: Darcy flux in the reference frame.
     *
     * \f[ -\int_\text{SCVF} J\,\frac{\kappa}{\mu_f}\,\mathbf{C}^{-1}\nabla_X p_f \cdot \mathbf{n}\,\mathrm{d}A_0 \f]
     */
    NumEqVector fluxIntegral(const FVElementGeometry& fvGeometry,
                             const ElementVariables& elemVars,
                             const SubControlVolumeFace& scvf) const
    {
        const auto& problem = this->asImp().problem();
        const auto& n = scvf.unitOuterNormal();
        NumEqVector f(0.0);

        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scvf))
        {
            const auto& ipCache = cache(elemVars, qpData.ipData());
            const auto& ipGlobal = qpData.ipData().global();

            // gradient of p_f at quadrature point
            GlobalPosition gradPf(0.0);
            for (const auto& localDof : localDofs(fvGeometry))
                gradPf.axpy(elemVars[localDof].fluidPressure(), ipCache.gradN(localDof.index()));

            // deformation gradient from coupling manager (curSol, managed by assembler)
            const Tensor F = problem.couplingManager().deformationGradientAtPoint(
                fvGeometry, ipGlobal);
            const Scalar J = F.determinant();

            // F^{-1} and C^{-1} = F^{-1} F^{-T}
            Tensor Finv(F);
            Finv.invert();

            // permeability from spatial params (may depend on J, φ)
            const auto kappa = problem.spatialParams().permeabilityAtPoint(
                fvGeometry.element(), fvGeometry, ipGlobal, J);
            const Scalar muInv = 1.0 / problem.spatialParams().fluidViscosity(
                fvGeometry.element(), fvGeometry);

            // pulled-back Darcy flux: q_0 = J * (kappa/muf) * C^{-1} * grad_X p_f
            GlobalPosition FinvT_gradPf(0.0);
            Finv.mtv(gradPf, FinvT_gradPf);         // F^{-T} * gradPf
            GlobalPosition q_ref(0.0);
            Finv.mv(FinvT_gradPf, q_ref);            // F^{-1} * F^{-T} * gradPf = C^{-1} * gradPf
            q_ref *= J * kappa * muInv;

            // outward normal flux contribution
            f[0] -= q_ref * n * qpData.weight();
        }
        return f;
    }

    //! No source term by default; override in a subclass if mass sources are needed.
    template<class SubControlVolume_>
    NumEqVector sourceIntegral(const FVElementGeometry&, const ElementVariables&,
                               const SubControlVolume_&) const
    { return NumEqVector(0.0); }
};

} // end namespace Dumux
#endif
