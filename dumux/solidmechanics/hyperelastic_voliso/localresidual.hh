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
 *   \f[ -\nabla_X \cdot \mathbf{P} = \mathbf{0} \f]
 *   where the stress tensor \f$ \mathbf{P} \f$ is supplied by the problem.
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
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

namespace Dumux {

/*!
 * \ingroup HyperelasticVolIso
 * \brief Momentum local residual for the mixed u–p hyperelastic model.
 *
 * Uses the Experimental::GridVariables concept so that `ElementVariables`
 * covers all local DOFs. The stress tensor is obtained from the problem via
 * `firstPiolaKirchhoffStressTensor(...)`. Neumann boundary conditions and
 * non-CV DOFs are also handled.
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
     * \brief FE contributions for non-CV DOFs.
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
     * \brief Interior flux integral.
     *
     * Called by Experimental::CVFELocalResidual::evalFlux for interior faces.
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
 * \ingroup HyperelasticVolIso
 * \brief Pressure local residual for the mixed u–p hyperelastic model.
 *
 * Enforces the algebraic pressure constraint
 * \f[ p_s^{\mathrm{eq}}(J) - p_s = 0 \f]
 * averaged over each SCV. The equilibrium pressure \f$ p_s^{\mathrm{eq}}(J) \f$
 * is supplied by the pressure problem via `problem.volumetricPressure(J)`.
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
     * \brief Pressure constraint integrated over an SCV.
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
        NumEqVector source(0.0);

        const auto& problem = this->asImp().problem();
        const auto& localBasis = fvGeometry.gridGeometry().feCache().get(fvGeometry.element().type()).localBasis();
        std::vector< Dune::FieldVector<Scalar, 1> > shapeValues;

        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scv, QuadratureRules::DuneQuadrature<3>{}))
        {
            const Tensor F = problem.couplingManager().deformationGradientAtPoint(
                fvGeometry, qpData.ipData().global());
            const Scalar J = F.determinant();
            const Scalar psEq = problem.volumetricPressure(J);

            localBasis.evaluateFunction(qpData.ipData().local(), shapeValues);
            Scalar ps = 0.0;
            for (const auto& localDof : localDofs(fvGeometry))
                ps += elemVars[localDof].pressure() * shapeValues[localDof.index()][0];

            source[0] -= (psEq - ps) * qpData.weight();
        }

        return source;
    }
};

} // end namespace Dumux
#endif
