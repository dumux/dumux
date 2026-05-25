// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup HyperelasticVolIso
 * \brief Local residuals for the volumetric-isochoric split hyperelastic mixed u–p model.
 *
 * Two local residual classes are provided:
 *
 * - `HyperelasticVolIsoMomentumLocalResidual`: Momentum equation
 *   \f[ -\nabla_X \cdot \mathbf{P} = \mathbf{0}, \quad
 *       \mathbf{P} = \mu(\mathbf{F} - \mathbf{F}^{-T}) + J p_s \mathbf{F}^{-T} \f]
 *
 * - `HyperelasticVolIsoPressureLocalResidual`: Algebraic pressure constraint
 *   \f[ p_s^{\mathrm{eq}}(J) - p_s = 0 \f]
 *   where \f$ p_s^{\mathrm{eq}}(J) \f$ is supplied by `problem.volumetricPressure(J)`.
 */
#ifndef DUMUX_SOLIDMECHANICS_HYPERELASTIC_VOLISO_LOCAL_RESIDUAL_HH
#define DUMUX_SOLIDMECHANICS_HYPERELASTIC_VOLISO_LOCAL_RESIDUAL_HH

#include <dune/common/fmatrix.hh>
#include <dune/grid/common/intersectioniterator.hh>

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

namespace Dumux {

/*!
 * \ingroup HyperelasticVolIso
 * \brief Momentum local residual for the mixed u–p hyperelastic model.
 *
 * Uses the Experimental::GridVariables concept so that `ElementVariables`
 * covers all local DOFs (vertex SCVs + bubble / edge-midpoint DOFs).
 * The first Piola–Kirchhoff stress is
 * \f[ \mathbf{P} = \mu(\mathbf{F} - \mathbf{F}^{-T}) + J p_s \mathbf{F}^{-T}, \f]
 * with the pressure \f$ p_s \f$ obtained from the coupling manager.
 *
 * The shear modulus \f$ \mu \f$ is read from `problem.spatialParams().shearModulus()`.
 * Neumann boundary conditions and non-CV (PQ2 edge-midpoint) DOFs are also handled.
 */
template<class TypeTag>
class HyperelasticVolIsoMomentumLocalResidual : public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
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

    NumEqVector computeStorage(const Problem&, const auto&, const auto&) const
    { return NumEqVector(0.0); }

    NumEqVector storageIntegral(const FVElementGeometry&, const ElementVariables&,
                                const auto&, bool) const
    { return NumEqVector(0.0); }

    NumEqVector sourceIntegral(const FVElementGeometry&, const ElementVariables&,
                               const auto&) const
    { return NumEqVector(0.0); }

    /*!
     * \brief FE contributions for non-CV (edge-midpoint) DOFs — only active for hybrid schemes (PQ2).
     *
     * For MINI (PQ1Bubble), `hasNonCVLocalDofsInterface` is false and this is a no-op.
     * For Taylor-Hood (PQ2), assembles:
     * \f[ \int_K \mathbf{P} : \nabla\phi_i \, dX
     *     + \int_{\partial K_N} \mathbf{T} \cdot \phi_i \, dA \f]
     * for each edge-midpoint DOF \f$ i \f$.
     */
    void addToElementFluxAndSourceResidual(ElementResidualVector& residual,
                                           const Problem& problem,
                                           const auto& element,
                                           const FVElementGeometry& fvGeometry,
                                           const ElementVariables& elemVars,
                                           const auto& elemBcTypes) const
    {
        if constexpr (Dumux::Detail::LocalDofs::hasNonCVLocalDofsInterface<FVElementGeometry>())
        {

        if (nonCVLocalDofs(fvGeometry).empty())
            return;

        const Scalar mu = problem.spatialParams().shearModulus();

        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, element))
        {
            const auto& ipCache = cache(elemVars, qpData.ipData());
            const auto& ipGlobal = qpData.ipData().global();

            Tensor F(0.0);
            for (const auto& localDof : localDofs(fvGeometry))
                for (int d = 0; d < dim; ++d)
                    F[d].axpy(elemVars[localDof].displacement(d), ipCache.gradN(localDof.index()));
            for (int d = 0; d < dim; ++d) F[d][d] += 1.0;

            const Scalar J  = F.determinant();
            const Scalar ps = problem.couplingManager().pressureAtPoint(fvGeometry, ipGlobal);

            Tensor Finv = F; Finv.invert();
            const Tensor FinvT = transpose(Finv);

            Tensor P(0.0);
            for (int i = 0; i < dim; ++i)
                for (int j = 0; j < dim; ++j)
                    P[i][j] = mu*(F[i][j] - FinvT[i][j]) + J*ps*FinvT[i][j];

            for (const auto& nonCVdof : nonCVLocalDofs(fvGeometry))
            {
                const auto idx = nonCVdof.index();
                for (int i = 0; i < dim; ++i)
                    for (int j = 0; j < dim; ++j)
                        residual[idx][i] += P[i][j] * ipCache.gradN(idx)[j] * qpData.weight();
            }
        }

        if (!elemBcTypes.hasNeumann())
            return;

        for (const auto& intersection : intersections(fvGeometry.gridGeometry().gridView(), element))
        {
            if (!intersection.boundary())
                continue;

            const auto& isecBcTypes = elemBcTypes.get(fvGeometry, intersection);
            if (!isecBcTypes.hasNeumann())
                continue;

            for (const auto& qpData : CVFE::quadratureRule(fvGeometry, intersection))
            {
                const auto& ipCache = cache(elemVars, qpData.ipData());
                const auto flux = problem.boundaryFlux(fvGeometry, elemVars, qpData.ipData());
                const auto& shapeValues = ipCache.shapeValues();

                for (const auto& nonCVdof : nonCVLocalDofs(fvGeometry))
                {
                    const auto idx = nonCVdof.index();
                    for (int eqIdx = 0; eqIdx < dim; ++eqIdx)
                        if (isecBcTypes.isNeumann(eqIdx))
                            residual[idx][eqIdx] += double(shapeValues[idx]) * flux[eqIdx] * qpData.weight();
                }
            }
        }
        } // end if constexpr
    }

    /*!
     * \brief Interior flux integral: \f$ -\mathbf{P} \cdot \mathbf{n} \f$ integrated over an scvf.
     *
     * Called by Experimental::CVFELocalResidual::evalFlux for interior faces.
     */
    NumEqVector fluxIntegral(const FVElementGeometry& fvGeometry,
                             const ElementVariables& elemVars,
                             const SubControlVolumeFace& scvf) const
    {
        const auto& problem = this->asImp().problem();
        const auto& n = scvf.unitOuterNormal();
        const Scalar mu = problem.spatialParams().shearModulus();
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

            const Scalar J  = F.determinant();
            const Scalar ps = problem.couplingManager().pressureAtPoint(fvGeometry, ipGlobal);

            Tensor Finv = F; Finv.invert();
            const Tensor FinvT = transpose(Finv);

            for (int i = 0; i < dim; ++i)
                for (int j = 0; j < dim; ++j)
                    f[Indices::momentum(i)] -=
                        (mu*(F[i][j] - FinvT[i][j]) + J*ps*FinvT[i][j])
                        * n[j] * qpData.weight();
        }
        return f;
    }
};

/*!
 * \ingroup HyperelasticVolIso
 * \brief Pressure local residual for the mixed u–p hyperelastic model.
 *
 * Enforces the algebraic pressure constraint
 * \f[ p_s^{\mathrm{eq}}(J) - p_s = 0 \f]
 * averaged over each SCV.  The equilibrium pressure \f$ p_s^{\mathrm{eq}}(J) \f$
 * is material-specific and must be supplied by the pressure problem via
 * `problem.volumetricPressure(J)`.  For example:
 * - \f$ W_{\mathrm{vol}} = \tfrac{\lambda}{2}(J-1)^2 \f$:
 *   \f$ p_s^{\mathrm{eq}} = \lambda(J-1) \f$
 * - Ogden \f$ \psi_1 \f$:
 *   \f$ p_s^{\mathrm{eq}} = \lambda(J^2-1)/(2J) \f$
 *
 * The deformation gradient \f$ \mathbf{F} \f$ is obtained from the coupling manager.
 */
template<class TypeTag>
class HyperelasticVolIsoPressureLocalResidual : public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
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

    NumEqVector computeStorage(const Problem&, const auto&, const auto&) const
    { return NumEqVector(0.0); }

    NumEqVector storageIntegral(const FVElementGeometry&, const ElementVariables&,
                                const auto&, bool) const
    { return NumEqVector(0.0); }

    NumEqVector fluxIntegral(const FVElementGeometry&, const ElementVariables&,
                             const auto&) const
    { return NumEqVector(0.0); }

    /*!
     * \brief Pressure constraint averaged over an SCV.
     *
     * Returns \f$ -(p_s^{\mathrm{eq}}(\bar J) - p_s) \f$ (negated because
     * `evalSource` subtracts the source, so the residual contribution is
     * \f$ p_s^{\mathrm{eq}} - p_s = 0 \f$ at equilibrium).
     */
    template<class SubControlVolume_>
    NumEqVector sourceIntegral(const FVElementGeometry& fvGeometry,
                               const ElementVariables& elemVars,
                               const SubControlVolume_& scv) const
    {
        const auto& problem = this->asImp().problem();
        const Scalar ps = elemVars[scv].pressure();

        Scalar integral = 0.0, totalWeight = 0.0;
        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scv))
        {
            const Tensor F = problem.couplingManager().deformationGradientAtPoint(
                fvGeometry, qpData.ipData().global());
            const Scalar J = F.determinant();
            const Scalar psEq = problem.volumetricPressure(J);
            integral    += (psEq - ps) * qpData.weight();
            totalWeight += qpData.weight();
        }

        NumEqVector source(0.0);
        source[0] = (totalWeight > 0.0) ? -integral/totalWeight : 0.0;
        return source;
    }
};

} // end namespace Dumux
#endif
