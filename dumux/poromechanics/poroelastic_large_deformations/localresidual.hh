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
 *   Quasi-static momentum balance \f$ -\nabla_X \cdot \mathbf{P} = \mathbf{f} \f$ with
 *   \f$ \mathbf{P} = \mathbf{P}_\text{iso}(\mathbf{F}) + (p_s - \alpha_B p_f) J \mathbf{F}^{-\top} \f$.
 *
 * - `PoroElasticLargeDefSolidPressureLocalResidual`:
 *   Algebraic solid-pressure constraint \f$p_s^{\mathrm{eq}}(J) - p_s = 0\f$.
 *
 * - `PoroElasticLargeDefFluidPressureLocalResidual`:
 *   Transient fluid mass balance (pulled back to reference frame):
 *   \f[ \frac{\partial}{\partial t}\!\left(\alpha_B J + S_p p_f\right)
 *       - \nabla_X \cdot \!\left(J\,\frac{\kappa}{\mu_f}\,\mathbf{C}^{-1}\nabla_X p_f\right) = q. \f]
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
#include <dumux/solidmechanics/hyperelastic_voliso/localresidual.hh>

namespace Dumux {

/*!
 * \ingroup PoroElasticLargeDef
 * \brief Momentum local residual for the poroelastic large-deformation model.
 *
 * The momentum residual of the mixed u–p hyperelastic model with the first Piola–Kirchhoff stress
 * \f[ \mathbf{P} = \mathbf{P}_\text{iso}(\mathbf{F}) + (p_s - \alpha_B p_f)\, J\, \mathbf{F}^{-\top}. \f]
 * The isochoric part \f$ \mathbf{P}_\text{iso} \f$ is the constitutive law supplied by
 * `problem.isochoricFirstPiolaKirchhoffStressTensor(F, elemDisc, ipGlobal)`.
 * The Biot coefficient is supplied by
 * `problem.spatialParams().biotCoefficient(element, elemDisc, ipGlobal)`.
 * The body force \f$ \mathbf{f} \f$ is the problem's `source`.
 */
template<class TypeTag>
class PoroElasticLargeDefMomentumLocalResidual
    : public HyperelasticVolIsoMomentumLocalResidual<TypeTag>
{
    using ParentType = HyperelasticVolIsoMomentumLocalResidual<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridGeometry::LocalView;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;

public:
    using ParentType::ParentType;

    //! The scalar \f$ p_s - \alpha_B p_f \f$ of the volumetric stress \f$ (p_s - \alpha_B p_f) J \mathbf{F}^{-\top} \f$
    Scalar volumetricStressPressure(const Problem& problem,
                                    const ElementDiscretization& elemDisc,
                                    const GlobalPosition& ipGlobal) const
    {
        const auto& couplingManager = problem.couplingManager();
        const Scalar ps = couplingManager.solidPressureAtPoint(elemDisc, ipGlobal);
        const Scalar pf = couplingManager.fluidPressureAtPoint(elemDisc, ipGlobal);
        const Scalar alphaB = problem.spatialParams().biotCoefficient(elemDisc.element(), elemDisc, ipGlobal);
        return ps - alphaB*pf;
    }
};

/*!
 * \ingroup PoroElasticLargeDef
 * \brief Solid pressure local residual: algebraic constraint \f$p_s^{\mathrm{eq}}(J) - p_s = 0\f$.
 *
 * The equilibrium solid pressure \f$p_s^{\mathrm{eq}}(J)\f$ is provided by the problem
 * via `problem.volumetricPressure(J)`. The constraint is integrated over each SCV, with
 * \f$p_s\f$ interpolated by the P1 shape functions at the quadrature points.
 */
template<class TypeTag>
class PoroElasticLargeDefSolidPressureLocalResidual
    : public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
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

    NumEqVector computeStorage(const auto&, const auto&, const auto&) const
    { return NumEqVector(0.0); }

    NumEqVector storageIntegral(const ElementDiscretization&, const ElementVariables&,
                                const auto&, bool) const
    { return NumEqVector(0.0); }

    NumEqVector fluxIntegral(const ElementDiscretization&, const ElementVariables&,
                             const auto&) const
    { return NumEqVector(0.0); }

    //! Enforces \f$p_s^{\mathrm{eq}}(J) - p_s = 0\f$ integrated over each SCV.
    template<class SubControlVolume_>
    NumEqVector sourceIntegral(const ElementDiscretization& elemDisc,
                               const ElementVariables& elemVars,
                               const SubControlVolume_& scv) const
    {
        NumEqVector source(0.0);

        const auto& problem = this->asImp().problem();
        const auto& localBasis = elemDisc.gridDiscretization().feCache().get(elemDisc.element().type()).localBasis();
        // thread_local keeps the capacity across the many calls of numeric differentiation
        thread_local std::vector< Dune::FieldVector<Scalar, 1> > shapeValues;

        for (const auto& qpData : CVFE::quadratureRule(elemDisc, scv, QuadratureRules::DuneQuadrature<3>{}))
        {
            const Tensor F = problem.couplingManager().deformationGradientAtPoint(
                elemDisc, qpData.ipData().global());
            const Scalar J = F.determinant();
            const Scalar psEq = problem.volumetricPressure(J);

            localBasis.evaluateFunction(qpData.ipData().local(), shapeValues);
            Scalar ps = 0.0;
            for (const auto& localDof : localDofs(elemDisc))
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
 * \f[ \frac{\partial}{\partial t}\!\left(\alpha_B J + S_p p_f\right)
 *     - \nabla_X \cdot \!\left(J\,\frac{\kappa}{\mu_f}\,\mathbf{C}^{-1}\nabla_X p_f\right) = q \f]
 * with the Biot coefficient \f$\alpha_B\f$, the storage coefficient \f$S_p = 1/M\f$
 * (inverse Biot modulus, zero for incompressible constituents) and the fluid source \f$q\f$
 * given by the problem's `source`.
 *
 * The spatial parameters provide
 * `biotCoefficient(element, elemDisc, ipGlobal)`,
 * `storageCoefficient(element, elemDisc, ipGlobal)`,
 * `permeability(element, elemDisc, ipGlobal, J)` and
 * `fluidViscosity(element, elemDisc)`.
 */
template<class TypeTag>
class PoroElasticLargeDefFluidPressureLocalResidual
    : public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridGeometry::LocalView;
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
     * \brief Storage term \f$ \int_\text{SCV} (\alpha_B J + S_p p_f) \, \mathrm{d}V_0 \f$
     *
     * Midpoint rule: \f$J\f$ and the coefficients are evaluated at the SCV center and
     * \f$p_f\f$ is lumped to the SCV dof. This is exact for piecewise linear displacements
     * and second-order accurate for higher-order displacements.
     * \f$J\f$ is evaluated from the coupling manager's current solution, which the
     * multi-stage assembler sets to the solution of the stage being evaluated, so the
     * `isPreviousTimeLevel` flag is not needed here.
     */
    template<class SubControlVolume_>
    NumEqVector storageIntegral(const ElementDiscretization& elemDisc,
                                const ElementVariables& elemVars,
                                const SubControlVolume_& scv,
                                bool /*isPreviousTimeLevel*/) const
    {
        const auto& problem = this->asImp().problem();
        const auto& spatialParams = problem.spatialParams();
        const auto& element = elemDisc.element();
        const auto& center = scv.center();

        const Scalar J = problem.couplingManager().deformationGradientAtPoint(elemDisc, center).determinant();
        const Scalar alphaB = spatialParams.biotCoefficient(element, elemDisc, center);
        const Scalar Sp = spatialParams.storageCoefficient(element, elemDisc, center);
        return NumEqVector((alphaB*J + Sp*elemVars[scv].fluidPressure())*scv.volume());
    }

    /*!
     * \brief Flux term: Darcy flux in the reference frame.
     *
     * \f[ -\int_\text{SCVF} J\,\frac{\kappa}{\mu_f}\,\mathbf{C}^{-1}\nabla_X p_f \cdot \mathbf{n}\,\mathrm{d}A_0 \f]
     */
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
            const auto& ipGlobal = qpData.ipData().global();

            GlobalPosition gradPf(0.0);
            for (const auto& localDof : localDofs(elemDisc))
                gradPf.axpy(elemVars[localDof].fluidPressure(), ipCache.gradN(localDof.index()));

            const Tensor F = problem.couplingManager().deformationGradientAtPoint(
                elemDisc, ipGlobal);
            const Scalar J = F.determinant();

            Tensor Finv(F);
            Finv.invert();

            const auto kappa = problem.spatialParams().permeability(
                elemDisc.element(), elemDisc, ipGlobal, J);
            const Scalar muInv = 1.0 / problem.spatialParams().fluidViscosity(
                elemDisc.element(), elemDisc);

            // q_0 = J (kappa/mu_f) C^{-1} grad_X p_f with C^{-1} = F^{-1} F^{-T}
            GlobalPosition FinvT_gradPf(0.0);
            Finv.mtv(gradPf, FinvT_gradPf);
            GlobalPosition q_ref(0.0);
            Finv.mv(FinvT_gradPf, q_ref);
            q_ref *= J * kappa * muInv;

            f[0] -= q_ref * n * qpData.weight();
        }
        return f;
    }
};

} // end namespace Dumux
#endif
