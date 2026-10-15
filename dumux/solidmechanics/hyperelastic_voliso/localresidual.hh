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
 *   \f[ -\nabla_X \cdot \mathbf{P} = \mathbf{f}, \quad
 *       \mathbf{P} = \mathbf{P}_\text{iso}(\mathbf{F}) + p_s J \mathbf{F}^{-\top}, \f]
 *   where the problem supplies the isochoric stress \f$ \mathbf{P}_\text{iso} \f$.
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
 * Assembles \f$ -\nabla_X \cdot \mathbf{P} = \mathbf{f} \f$ with the first Piola–Kirchhoff stress
 * \f[ \mathbf{P} = \mathbf{P}_\text{iso}(\mathbf{F}) + p\, J\, \mathbf{F}^{-\top}. \f]
 * The problem supplies \f$ \mathbf{P}_\text{iso} \f$ through
 * `isochoricFirstPiolaKirchhoffStressTensor(F, elemDisc, ipGlobal)` and the body force
 * \f$ \mathbf{f} \f$ through `source`. The scalar \f$ p \f$ is `volumetricStressPressure`, here the
 * pressure \f$ p_s \f$ of the pressure subdomain; derived residuals may add further pressures.
 * The residual is quasi-static (no storage term).
 */
template<class TypeTag>
class HyperelasticVolIsoMomentumLocalResidual : public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using GlobalPosition = typename GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
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

    NumEqVector storageIntegral(const ElementDiscretization&, const ElementVariables&,
                                const auto&, bool) const
    { return NumEqVector(0.0); }

    //! FE contributions of stress divergence and body force for non-CV dofs (e.g. PQ2 hybrid)
    void addToElementFluxAndSourceResidual(ElementResidualVector& residual,
                                           const Problem& problem,
                                           const auto& element,
                                           const ElementDiscretization& elemDisc,
                                           const ElementVariables& elemVars) const
    {
        if constexpr (Dumux::Detail::LocalDofs::hasNonCVLocalDofsInterface<ElementDiscretization>())
        {
            if (nonCVLocalDofs(elemDisc).empty())
                return;

            for (const auto& qpData : CVFE::quadratureRule(elemDisc, element))
            {
                const auto& ipCache = cache(elemVars, qpData.ipData());
                const auto& ipGlobal = qpData.ipData().global();

                const auto F = deformationGradient_(elemDisc, elemVars, ipCache);
                const auto P = firstPiolaKirchhoffStressTensor_(problem, F, elemDisc, ipGlobal);
                const auto f = problem.source(elemDisc, elemVars, qpData.ipData());
                const auto& shapeValues = ipCache.shapeValues();

                for (const auto& nonCVdof : nonCVLocalDofs(elemDisc))
                {
                    const auto idx = nonCVdof.index();
                    for (int i = 0; i < dim; ++i)
                    {
                        for (int j = 0; j < dim; ++j)
                            residual[idx][Indices::momentum(i)] += P[i][j] * ipCache.gradN(idx)[j] * qpData.weight();
                        residual[idx][Indices::momentum(i)] -= f[Indices::momentum(i)] * shapeValues[idx][0] * qpData.weight();
                    }
                }
            }
        }
    }

    //! Interior flux integral: \f$ -\int \mathbf{P} \cdot \mathbf{n} \, \mathrm{d}A \f$
    template<class SubControlVolumeFace>
    NumEqVector fluxIntegral(const ElementDiscretization& elemDisc,
                             const ElementVariables& elemVars,
                             const SubControlVolumeFace& scvf) const
    {
        const auto& problem = this->asImp().problem();
        const auto& n = scvf.unitOuterNormal();
        NumEqVector f(0.0);

        for (const auto& qpData : CVFE::quadratureRule(elemDisc, scvf))
        {
            const auto& ipCache = cache(elemVars, qpData.ipData());
            const auto F = deformationGradient_(elemDisc, elemVars, ipCache);
            const auto P = firstPiolaKirchhoffStressTensor_(problem, F, elemDisc, qpData.ipData().global());

            for (int i = 0; i < dim; ++i)
                for (int j = 0; j < dim; ++j)
                    f[Indices::momentum(i)] -= P[i][j] * n[j] * qpData.weight();
        }
        return f;
    }

    //! The scalar \f$ p \f$ of the volumetric stress \f$ p J \mathbf{F}^{-\top} \f$
    Scalar volumetricStressPressure(const Problem& problem,
                                    const ElementDiscretization& elemDisc,
                                    const GlobalPosition& ipGlobal) const
    { return problem.couplingManager().pressureAtPoint(elemDisc, ipGlobal); }

private:
    template<class IpCache>
    Tensor deformationGradient_(const ElementDiscretization& elemDisc,
                                const ElementVariables& elemVars,
                                const IpCache& ipCache) const
    {
        Tensor F(0.0);
        for (const auto& localDof : localDofs(elemDisc))
            for (int d = 0; d < dim; ++d)
                F[d].axpy(elemVars[localDof].displacement(d), ipCache.gradN(localDof.index()));
        for (int d = 0; d < dim; ++d)
            F[d][d] += 1.0;
        return F;
    }

    Tensor firstPiolaKirchhoffStressTensor_(const Problem& problem,
                                            const Tensor& F,
                                            const ElementDiscretization& elemDisc,
                                            const GlobalPosition& ipGlobal) const
    {
        auto P = problem.isochoricFirstPiolaKirchhoffStressTensor(F, elemDisc, ipGlobal);

        Tensor Finv(F);
        Finv.invert();
        const Scalar pJ = this->asImp().volumetricStressPressure(problem, elemDisc, ipGlobal)*F.determinant();
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j)
                P[i][j] += pJ*Finv[j][i];

        return P;
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
    using ElementDiscretization = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Tensor = Dune::FieldMatrix<Scalar, GridView::dimension, GridView::dimension>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVariables = typename GridVariables::GridVariablesCache::LocalView;

public:
    using ParentType::ParentType;
    using ElementResidualVector = typename ParentType::ElementResidualVector;

    NumEqVector computeStorage(const Problem&, const auto&, const auto&) const
    { return NumEqVector(0.0); }

    NumEqVector storageIntegral(const ElementDiscretization&, const ElementVariables&,
                                const auto&, bool) const
    { return NumEqVector(0.0); }

    NumEqVector fluxIntegral(const ElementDiscretization&, const ElementVariables&,
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
    NumEqVector sourceIntegral(const ElementDiscretization& elemDisc,
                               const ElementVariables& elemVars,
                               const SubControlVolume_& scv) const
    {
        NumEqVector source(0.0);

        const auto& problem = this->asImp().problem();
        const auto& localBasis = elemDisc.gridDiscretization().feCache().get(elemDisc.element().type()).localBasis();
        std::vector< Dune::FieldVector<Scalar, 1> > shapeValues;

        for (const auto& qpData : CVFE::quadratureRule(elemDisc, scv, QuadratureRules::DuneQuadrature<3>{}))
        {
            const Tensor F = problem.couplingManager().deformationGradientAtPoint(
                elemDisc, qpData.ipData().global());
            const Scalar J = F.determinant();
            const Scalar psEq = problem.volumetricPressure(J);

            localBasis.evaluateFunction(qpData.ipData().local(), shapeValues);
            Scalar ps = 0.0;
            for (const auto& localDof : localDofs(elemDisc))
                ps += elemVars[localDof].pressure() * shapeValues[localDof.index()][0];

            source[0] -= (psEq - ps) * qpData.weight();
        }

        return source;
    }

    /*!
     * \brief FE contributions for non-CV DOFs.
     */
    void addToElementFluxAndSourceResidual(ElementResidualVector& residual,
                                           const Problem& problem,
                                           const auto& element,
                                           const ElementDiscretization& elemDisc,
                                           const ElementVariables& elemVars) const
    {
        if constexpr (Dumux::Detail::LocalDofs::hasNonCVLocalDofsInterface<ElementDiscretization>())
        {
            if (nonCVLocalDofs(elemDisc).empty())
                return;

            const auto& problem = this->asImp().problem();
            const auto& localBasis = elemDisc.gridDiscretization().feCache().get(elemDisc.element().type()).localBasis();
            std::vector< Dune::FieldVector<Scalar, 1> > shapeValues;

            for (const auto& qpData : CVFE::quadratureRule(elemDisc, element, QuadratureRules::DuneQuadrature<3>{}))
            {
                const Tensor F = problem.couplingManager().deformationGradientAtPoint(
                    elemDisc, qpData.ipData().global());
                const Scalar J = F.determinant();
                const Scalar psEq = problem.volumetricPressure(J);

                localBasis.evaluateFunction(qpData.ipData().local(), shapeValues);
                Scalar ps = 0.0;
                for (const auto& localDof : localDofs(elemDisc))
                    ps += elemVars[localDof].pressure() * shapeValues[localDof.index()][0];

                for (const auto& nonCVdof : nonCVLocalDofs(elemDisc))
                {
                    const auto idx = nonCVdof.index();
                    residual[idx][0] -= (psEq - ps) * shapeValues[idx][0] * qpData.weight();
                }
            }

        }
    }

};

} // end namespace Dumux
#endif
