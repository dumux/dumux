// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoroElasticLargeDef
 * \brief Local residuals for the large-deformation poroelastic mixed
 *        \f$ \mathbf{d} \f$–\f$ p_s \f$–\f$ p_f \f$ model.
 */
#ifndef DUMUX_POROMECHANICS_POROELASTIC_LARGE_DEFORMATIONS_LOCAL_RESIDUAL_HH
#define DUMUX_POROMECHANICS_POROELASTIC_LARGE_DEFORMATIONS_LOCAL_RESIDUAL_HH

#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/solidmechanics/hyperelastic_voliso/localresidual.hh>

namespace Dumux {

/*!
 * \ingroup PoroElasticLargeDef
 * \brief Momentum local residual for the poroelastic large-deformation model.
 *
 * The momentum residual of the mixed \f$ \mathbf{d} \f$–\f$ p_s \f$ hyperelastic model with the first
 * Piola–Kirchhoff stress
 * \f[ \mathbf{P} = \mathbf{P}_\text{iso}(\mathbf{F}) + (p_s - \alpha_B p_f)\, J\, \mathbf{F}^{-\top}, \f]
 * where the problem supplies \f$ \mathbf{P}_\text{iso} \f$ via
 * `isochoricFirstPiolaKirchhoffStressTensor(F, elemDisc, ipGlobal)` and the spatial parameters
 * supply \f$ \alpha_B \f$ via `biotCoefficient(element, elemDisc, ipGlobal)`.
 */
template<class TypeTag>
class PoroElasticLargeDefMomentumLocalResidual
    : public HyperelasticVolIsoMomentumLocalResidual<TypeTag>
{
    using ParentType = HyperelasticVolIsoMomentumLocalResidual<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridDiscretization::LocalView;

public:
    using ParentType::ParentType;

    //! The scalar \f$ p_s - \alpha_B p_f \f$ of the volumetric stress \f$ (p_s - \alpha_B p_f) J \mathbf{F}^{-\top} \f$
    template<class IpData>
    Scalar volumetricStressPressure(const Problem& problem,
                                    const ElementDiscretization& elemDisc,
                                    const IpData& ipData) const
    {
        const auto& couplingManager = problem.couplingManager();
        const Scalar ps = couplingManager.solidPressure(elemDisc, ipData);
        const Scalar pf = couplingManager.fluidPressure(elemDisc, ipData);
        const Scalar alphaB = problem.spatialParams().biotCoefficient(elemDisc.element(), elemDisc, ipData.global());
        return ps - alphaB*pf;
    }
};

/*!
 * \ingroup PoroElasticLargeDef
 * \brief Solid pressure local residual: algebraic constraint \f$p_s^{\mathrm{eq}}(J) - p_s = 0\f$.
 *
 * The equilibrium solid pressure \f$p_s^{\mathrm{eq}}(J)\f$ is provided by the problem
 * via `problem.volumetricPressure(J)`. The constraint is integrated over each SCV with the
 * SCV quadrature rule of the grid geometry, with \f$p_s\f$ interpolated by the shape
 * functions at the quadrature points.
 */
template<class TypeTag>
class PoroElasticLargeDefSolidPressureLocalResidual
    : public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridDiscretization::LocalView;
    using GridView = typename GridDiscretization::GridView;
    using Tensor = Dune::FieldMatrix<Scalar, GridView::dimension, GridView::dimension>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVariables = typename GridVariables::GridVariablesCache::LocalView;

public:
    using ParentType::ParentType;
    using ElementResidualVector = typename ParentType::ElementResidualVector;

    NumEqVector storageIntegral(const ElementDiscretization&, const ElementVariables&,
                                const auto&, bool) const
    { return NumEqVector(0.0); }

    NumEqVector fluxIntegral(const ElementDiscretization&, const ElementVariables&,
                             const auto&) const
    { return NumEqVector(0.0); }

    NumEqVector sourceIntegral(const ElementDiscretization& elemDisc,
                               const ElementVariables& elemVars,
                               const auto& scv) const
    {
        NumEqVector source(0.0);

        const auto& problem = this->asImp().problem();
        const auto& localBasis = elemDisc.gridDiscretization().feCache().get(elemDisc.element().type()).localBasis();
        thread_local std::vector< Dune::FieldVector<Scalar, 1> > shapeValues;

        for (const auto& qpData : CVFE::quadratureRule(elemDisc, scv))
        {
            const Tensor F = problem.couplingManager().deformationGradient(elemDisc, qpData.ipData());
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
    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridDiscretization::LocalView;
    using SubControlVolumeFace = typename GridDiscretization::SubControlVolumeFace;
    using GridView = typename GridDiscretization::GridView;
    using GlobalPosition = typename GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
    using Tensor = Dune::FieldMatrix<Scalar, GridView::dimension, GridView::dimension>;
    static constexpr int dim = GridView::dimension;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVariables = typename GridVariables::GridVariablesCache::LocalView;

public:
    using ParentType::ParentType;
    using ElementResidualVector = typename ParentType::ElementResidualVector;

    /*!
     * \brief Storage term \f$ \int_\text{SCV} (\alpha_B J + S_p p_f) \, \mathrm{d}V_0 \f$
     *
     * \f$J\f$ and the coefficients are evaluated at the quadrature points of the SCV rule of the
     * grid geometry and \f$p_f\f$ is lumped to the SCV dof. \f$J\f$ comes from the coupling manager's
     * solution, which is the one of the evaluated stage, so `isPreviousTimeLevel` is not needed.
     */
    NumEqVector storageIntegral(const ElementDiscretization& elemDisc,
                                const ElementVariables& elemVars,
                                const auto& scv,
                                bool /*isPreviousTimeLevel*/) const
    {
        const auto& problem = this->asImp().problem();
        const auto& spatialParams = problem.spatialParams();
        const auto& element = elemDisc.element();
        const Scalar pf = elemVars[scv].fluidPressure();

        NumEqVector storage(0.0);
        for (const auto& qpData : CVFE::quadratureRule(elemDisc, scv))
        {
            const auto& ipGlobal = qpData.ipData().global();
            const Scalar J = problem.couplingManager().deformationGradient(elemDisc, qpData.ipData()).determinant();
            const Scalar alphaB = spatialParams.biotCoefficient(element, elemDisc, ipGlobal);
            const Scalar Sp = spatialParams.storageCoefficient(element, elemDisc, ipGlobal);
            storage[0] += (alphaB*J + Sp*pf)*qpData.weight();
        }
        return storage;
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

            const Tensor F = problem.couplingManager().deformationGradient(elemDisc, qpData.ipData());
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
