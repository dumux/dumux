// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Hyperelastic
 * \brief Local residual of the hyperelastic model
 */
#ifndef DUMUX_SOLIDMECHANICS_HYPERELASTIC_LOCAL_RESIDUAL_HH
#define DUMUX_SOLIDMECHANICS_HYPERELASTIC_LOCAL_RESIDUAL_HH

#include <ranges>
#include <utility>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

namespace Dumux {

/*!
 * \ingroup Hyperelastic
 * \brief Local residual of the hyperelastic model
 *
 * Balance of momentum in the reference configuration,
 * \f[ \varrho_0 \ddot{\mathbf{u}} - \nabla \cdot \mathbf{P} = \mathbf{f}, \f]
 * with the first Piola-Kirchhoff stress \f$ \mathbf{P} \f$ and the body force \f$ \mathbf{f} \f$.
 * Degrees of freedom with a control volume are tested with its indicator function, all others
 * with their basis function.
 *
 * The problem provides
 * - the stress, as `firstPiolaKirchhoffStressTensor(F, elemDisc, ipData)` or
 *   `firstPiolaKirchhoffStressTensor(F)`,
 * - the body force as `source(elemDisc, elemVars, ipData)`,
 * - and the boundary terms as `boundaryFlux(elemDisc, elemVars, ipData)`, which enters the residual
 *   as \f$ \int \mathbf{g} \cdot \mathbf{w} \f$, i.e. \f$ \mathbf{g} = -\mathbf{P}\mathbf{N} \f$
 *   for a prescribed traction \f$ \mathbf{P}\mathbf{N} \f$.
 *
 * For dynamics, the problem additionally provides the acceleration and velocity of a local dof,
 * `acceleration(elemDisc, elemVars, localDof)` and `velocity(elemDisc, elemVars, localDof)`, as
 * functions of its current displacement (for instance from `Experimental::NewmarkBeta`), and the
 * density `solidDensity(elemDisc, ipData)`. The residual is then assembled as a stationary one.
 * With a velocity, the stress is requested as
 * `firstPiolaKirchhoffStressTensor(F, dF/dt, elemDisc, ipData)`, which allows rate-dependent laws.
 * Degrees of freedom with a control volume use the acceleration of the degree of freedom (lumped
 * mass), all others the interpolated acceleration at the element quadrature points (consistent mass).
 *
 * Alternatively, a problem enables the inertia operator \f$ \varrho_0 \mathbf{u} \f$ with
 * `static constexpr bool enableInertia()` for the time integration methods with second time
 * derivatives (see Experimental::ImplicitSolveParams), with the same mass lumping and the density
 * `solidDensity(elemDisc, ipData)`. Such problems do not provide the acceleration.
 */
template<class TypeTag>
class HyperelasticLocalResidual
: public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridDiscretization::LocalView;
    using GridView = typename GridDiscretization::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVariables = typename GridVariables::GridVariablesCache::LocalView;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using LocalDof = std::ranges::range_value_t<decltype(localDofs(std::declval<const ElementDiscretization&>()))>;

    static constexpr int dim = GridView::dimension;
    using Tensor = Dune::FieldMatrix<Scalar, dim, dim>;
    using Vector = Dune::FieldVector<Scalar, dim>;

    static constexpr bool hasInertia = requires(const Problem& p, const ElementDiscretization& ed,
                                                const ElementVariables& ev, const LocalDof& localDof)
    { p.acceleration(ed, ev, localDof); };

    static constexpr bool hasVelocity = requires(const Problem& p, const ElementDiscretization& ed,
                                                 const ElementVariables& ev, const LocalDof& localDof)
    { p.velocity(ed, ev, localDof); };

public:
    using ParentType::ParentType;
    using ElementResidualVector = typename ParentType::ElementResidualVector;

    template<class SubControlVolume>
    NumEqVector storageIntegral(const ElementDiscretization&, const ElementVariables&,
                                const SubControlVolume&, bool) const
    { return NumEqVector(0.0); }

    /*!
     * \brief The inertia \f$ \int \varrho_0 \mathbf{u} \f$ over a control volume, mass lumped
     * \note Used if the problem enables inertia with `static constexpr bool enableInertia()`.
     */
    template<class SubControlVolume>
    NumEqVector inertiaIntegral(const ElementDiscretization& elemDisc,
                                const ElementVariables& elemVars,
                                const SubControlVolume& scv) const
    {
        static_assert(!hasInertia, "Problems that enable the inertia operator do not provide the acceleration");
        const auto& problem = this->asImp().problem();
        Scalar mass = 0.0;
        for (const auto& qpData : CVFE::quadratureRule(elemDisc, scv))
            mass += qpData.weight()*problem.solidDensity(elemDisc, qpData.ipData());

        NumEqVector inertia(elemVars[scv].displacement());
        inertia *= mass*elemVars[scv].extrusionFactor();
        return inertia;
    }

    /*!
     * \brief The inertia \f$ \int \varrho_0 \mathbf{u} \cdot \mathbf{w} \f$ of the degrees of
     *        freedom without control volume, with the consistent mass
     */
    void addToElementInertia(ElementResidualVector& inertia,
                             const Problem& problem,
                             const Element& element,
                             const ElementDiscretization& elemDisc,
                             const ElementVariables& elemVars) const
    {
        if constexpr (Detail::LocalDofs::hasNonCVLocalDofsInterface<ElementDiscretization>())
        {
            if (std::ranges::empty(nonCVLocalDofs(elemDisc)))
                return;

            for (const auto& qpData : CVFE::quadratureRule(elemDisc, element))
            {
                const auto& ipData = qpData.ipData();
                const auto& shapeValues = cache(elemVars, ipData).shapeValues();
                Vector u(0.0);
                for (const auto& localDof : localDofs(elemDisc))
                    u.axpy(shapeValues[localDof.index()][0], elemVars[localDof].displacement());

                const auto rho = problem.solidDensity(elemDisc, ipData);
                for (const auto& localDof : nonCVLocalDofs(elemDisc))
                    inertia[localDof.index()].axpy(qpData.weight()*rho*shapeValues[localDof.index()][0], u);
            }
        }
    }

    //! Body force minus inertia (lumped) over a control volume
    template<class SubControlVolume>
    NumEqVector sourceIntegral(const ElementDiscretization& elemDisc,
                               const ElementVariables& elemVars,
                               const SubControlVolume& scv) const
    {
        const auto& problem = this->asImp().problem();
        NumEqVector source(0.0);
        for (const auto& qpData : CVFE::quadratureRule(elemDisc, scv))
        {
            const auto& ipData = qpData.ipData();
            NumEqVector value = problem.source(elemDisc, elemVars, ipData);
            if constexpr (hasInertia)
            {
                const LocalDof localDof(scv.localDofIndex(), scv.dofIndex(), scv.elementIndex());
                value.axpy(-problem.solidDensity(elemDisc, ipData), problem.acceleration(elemDisc, elemVars, localDof));
            }
            source.axpy(qpData.weight(), value);
        }
        source *= elemVars[scv].extrusionFactor();
        return source;
    }

    //! Stress flux \f$ -\int \mathbf{P}\mathbf{N} \f$ over a control volume face
    template<class SubControlVolumeFace>
    NumEqVector fluxIntegral(const ElementDiscretization& elemDisc,
                             const ElementVariables& elemVars,
                             const SubControlVolumeFace& scvf) const
    {
        const auto& problem = this->asImp().problem();
        NumEqVector flux(0.0);
        for (const auto& qpData : CVFE::quadratureRule(elemDisc, scvf))
        {
            const auto& ipData = qpData.ipData();
            const auto P = stress_(problem, elemDisc, elemVars, cache(elemVars, ipData), ipData);
            Vector PN;
            P.mv(scvf.unitOuterNormal(), PN);
            flux.axpy(-qpData.weight(), PN);
        }
        flux *= elemVars[elemDisc.scv(scvf.insideScvIdx())].extrusionFactor();
        return flux;
    }

    //! Galerkin terms \f$ \int \mathbf{P} : \nabla \mathbf{w} + (\varrho_0 \ddot{\mathbf{u}} - \mathbf{f}) \cdot \mathbf{w} \f$ of the dofs without control volume
    void addToElementFluxAndSourceResidual(ElementResidualVector& residual,
                                           const Problem& problem,
                                           const Element& element,
                                           const ElementDiscretization& elemDisc,
                                           const ElementVariables& elemVars) const
    {
        if constexpr (Detail::LocalDofs::hasNonCVLocalDofsInterface<ElementDiscretization>())
        {
            if (std::ranges::empty(nonCVLocalDofs(elemDisc)))
                return;

            for (const auto& qpData : CVFE::quadratureRule(elemDisc, element))
            {
                const auto& ipData = qpData.ipData();
                const auto& ipCache = cache(elemVars, ipData);
                const auto P = stress_(problem, elemDisc, elemVars, ipCache, ipData);
                const auto& shapeValues = ipCache.shapeValues();

                NumEqVector bodyForce = problem.source(elemDisc, elemVars, ipData);
                if constexpr (hasInertia)
                {
                    const auto rho = problem.solidDensity(elemDisc, ipData);
                    for (const auto& localDof : localDofs(elemDisc))
                        bodyForce.axpy(-rho*shapeValues[localDof.index()][0], problem.acceleration(elemDisc, elemVars, localDof));
                }

                for (const auto& localDof : nonCVLocalDofs(elemDisc))
                {
                    const auto idx = localDof.index();
                    Vector value;
                    P.mv(ipCache.gradN(idx), value);
                    value.axpy(-shapeValues[idx][0], bodyForce);
                    residual[idx].axpy(qpData.weight(), value);
                }
            }
        }
    }

private:
    template<class ShapeData, class IpData>
    Tensor stress_(const Problem& problem,
                   const ElementDiscretization& elemDisc,
                   const ElementVariables& elemVars,
                   const ShapeData& shapeData,
                   const IpData& ipData) const
    {
        Tensor F(0.0);
        for (const auto& localDof : localDofs(elemDisc))
            for (int i = 0; i < dim; ++i)
                F[i].axpy(elemVars[localDof].displacement(i), shapeData.gradN(localDof.index()));
        for (int i = 0; i < dim; ++i)
            F[i][i] += 1.0;

        if constexpr (hasVelocity)
        {
            Tensor Fdot(0.0);
            for (const auto& localDof : localDofs(elemDisc))
            {
                const auto v = problem.velocity(elemDisc, elemVars, localDof);
                for (int i = 0; i < dim; ++i)
                    Fdot[i].axpy(v[i], shapeData.gradN(localDof.index()));
            }
            return problem.firstPiolaKirchhoffStressTensor(F, Fdot, elemDisc, ipData);
        }
        else if constexpr (requires { problem.firstPiolaKirchhoffStressTensor(F, elemDisc, ipData); })
            return problem.firstPiolaKirchhoffStressTensor(F, elemDisc, ipData);
        else
            return problem.firstPiolaKirchhoffStressTensor(F);
    }
};

} // end namespace Dumux

#endif
