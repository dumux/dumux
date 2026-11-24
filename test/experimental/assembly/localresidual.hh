// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief A nonlinear diffusion residual with a hand-written element Jacobian
 */
#ifndef DUMUX_TEST_EXPERIMENTAL_ASSEMBLY_LOCAL_RESIDUAL_HH
#define DUMUX_TEST_EXPERIMENTAL_ASSEMBLY_LOCAL_RESIDUAL_HH

#include <dumux/common/numeqvector.hh>
#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

namespace Dumux {

/*!
 * \brief Nonlinear diffusion, \f$ -\nabla\cdot(D(u)\nabla u) = q \f$ with \f$ D(u) = 1 + u^2 \f$
 *
 * The point of the nonlinearity is that the Jacobian is not the residual's own matrix, so an
 * error in `addElementJacobian` cannot hide.
 */
template<class TypeTag>
class NonlinearDiffusionLocalResidual : public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridVariablesCache = typename GridVariables::GridVariablesCache;
    using ElementVariables = typename GridVariablesCache::LocalView;

public:
    using ElementResidualVector = typename ParentType::ElementResidualVector;
    using ParentType::ParentType;

    static Scalar diffusionCoefficient(Scalar u)
    { return 1.0 + u*u; }

    static Scalar diffusionCoefficientDerivative(Scalar u)
    { return 2.0*u; }

    //! \f$ \int_{\sigma} -D(u)\nabla u \cdot \mathbf{n} \f$
    NumEqVector fluxIntegral(const FVElementGeometry& fvGeometry,
                             const ElementVariables& elemVars,
                             const SubControlVolumeFace& scvf) const
    {
        NumEqVector flux(0.0);
        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scvf))
        {
            const auto& ipData = qpData.ipData();
            const auto& ipCache = cache(elemVars, ipData);

            Scalar u = 0.0;
            GlobalPosition gradU(0.0);
            interpolate_(fvGeometry, elemVars, ipCache, u, gradU);

            flux[0] -= qpData.weight()*diffusionCoefficient(u)*(gradU*ipData.unitOuterNormal());
        }
        return flux;
    }

    /*!
     * \brief The Galerkin-tested dofs of a hybrid scheme: \f$ \int_E D(u)\nabla u\cdot\nabla N_i
     *        - \int_E q N_i \f$
     */
    void addToElementFluxAndSourceResidual(ElementResidualVector& residual,
                                           const Problem& problem,
                                           const Element& element,
                                           const FVElementGeometry& fvGeometry,
                                           const ElementVariables& elemVars) const
    {
        if constexpr (Dumux::Detail::LocalDofs::hasNonCVLocalDofsInterface<FVElementGeometry>())
        {
            if (nonCVLocalDofs(fvGeometry).empty())
                return;

            for (const auto& qpData : CVFE::quadratureRule(fvGeometry, element))
            {
                const auto& ipData = qpData.ipData();
                const auto& ipCache = cache(elemVars, ipData);

                Scalar u = 0.0;
                GlobalPosition gradU(0.0);
                interpolate_(fvGeometry, elemVars, ipCache, u, gradU);
                const auto source = problem.source(fvGeometry, elemVars, ipData);

                for (const auto& localDof : nonCVLocalDofs(fvGeometry))
                {
                    const auto i = localDof.index();
                    residual[i][0] += qpData.weight()*(
                        diffusionCoefficient(u)*(gradU*ipCache.gradN(i))
                        - ipCache.shapeValues()[i][0]*source[0]
                    );
                }
            }
        }
    }

    /*!
     * \brief \f$ \partial R_i / \partial u_j \f$ for every pair of local dofs of this element
     *
     * The block is the assembler's — zeroed before the call, scattered into the global matrix
     * after it — and is indexed by *local* dof index.
     *
     * Differentiating the flux integral above,
     * \f[
     *   \frac{\partial F}{\partial u_j} = \int_\sigma
     *      -D'(u) N_j (\nabla u \cdot \mathbf{n}) - D(u) (\nabla N_j \cdot \mathbf{n}),
     * \f]
     * and the two terms are scattered exactly as `CVFELocalResidual::evalFlux` scatters the flux
     * itself: onto the inside dof, and negated onto the outside dof unless the face lies inside
     * an overlapping control volume.
     *
     * The source here depends on position only, so it contributes nothing; a solution-dependent
     * source would be differentiated in the problem's `addSourceDerivatives`.
     */
    template<class ElementJacobianMatrix>
    void addFluxAndSourceJacobian(ElementJacobianMatrix& elemJac,
                                  const Problem& problem,
                                  const Element& element,
                                  const FVElementGeometry& fvGeometry,
                                  const ElementVariables& elemVars) const
    {
        ElementResidualVector dFlux(Dumux::Detail::LocalDofs::numLocalDofs(fvGeometry));

        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (scvf.boundary())
                continue;

            dFlux = 0.0;
            for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scvf))
            {
                const auto& ipData = qpData.ipData();
                const auto& ipCache = cache(elemVars, ipData);

                Scalar u = 0.0;
                GlobalPosition gradU(0.0);
                interpolate_(fvGeometry, elemVars, ipCache, u, gradU);

                const auto& normal = ipData.unitOuterNormal();
                const auto gradUn = gradU*normal;
                const auto d = diffusionCoefficient(u);
                const auto dd = diffusionCoefficientDerivative(u);

                for (const auto& localDof : localDofs(fvGeometry))
                {
                    const auto shape = ipCache.shapeValues()[localDof.index()][0];
                    const auto gradNn = ipCache.gradN(localDof.index())*normal;
                    dFlux[localDof.index()][0] -= qpData.weight()*(dd*shape*gradUn + d*gradNn);
                }
            }

            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
            const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());

            bool twoSided = true;
            if constexpr (requires { scvf.isOverlapping(); })
                twoSided = !scvf.isOverlapping();

            for (const auto& localDof : localDofs(fvGeometry))
            {
                const auto derivative = dFlux[localDof.index()][0];
                elemJac[insideScv.localDofIndex()][localDof.index()][0][0] += derivative;
                if (twoSided)
                    elemJac[outsideScv.localDofIndex()][localDof.index()][0][0] -= derivative;
            }
        }

        addGalerkinJacobian_(elemJac, fvGeometry, elemVars);
    }

private:
    /*!
     * \brief The rows a per-sub-control-volume hook could never reach
     *
     * \f$ \partial/\partial u_j \int_E D(u)\nabla u\cdot\nabla N_i
     *   = \int_E D'(u) N_j (\nabla u\cdot\nabla N_i) + D(u)(\nabla N_j\cdot\nabla N_i) \f$
     */
    template<class ElementJacobianMatrix>
    void addGalerkinJacobian_(ElementJacobianMatrix& elemJac,
                              const FVElementGeometry& fvGeometry,
                              const ElementVariables& elemVars) const
    {
        if constexpr (Dumux::Detail::LocalDofs::hasNonCVLocalDofsInterface<FVElementGeometry>())
        {
            if (nonCVLocalDofs(fvGeometry).empty())
                return;

            for (const auto& qpData : CVFE::quadratureRule(fvGeometry, fvGeometry.element()))
            {
                const auto& ipCache = cache(elemVars, qpData.ipData());

                Scalar u = 0.0;
                GlobalPosition gradU(0.0);
                interpolate_(fvGeometry, elemVars, ipCache, u, gradU);
                const auto d = diffusionCoefficient(u);
                const auto dd = diffusionCoefficientDerivative(u);

                for (const auto& testDof : nonCVLocalDofs(fvGeometry))
                {
                    const auto& gradNi = ipCache.gradN(testDof.index());
                    for (const auto& localDof : localDofs(fvGeometry))
                    {
                        const auto j = localDof.index();
                        elemJac[testDof.index()][j][0][0] += qpData.weight()*(
                            dd*ipCache.shapeValues()[j][0]*(gradU*gradNi)
                            + d*(ipCache.gradN(j)*gradNi)
                        );
                    }
                }
            }
        }
    }

    template<class IpCache>
    void interpolate_(const FVElementGeometry& fvGeometry,
                      const ElementVariables& elemVars,
                      const IpCache& ipCache,
                      Scalar& u, GlobalPosition& gradU) const
    {
        for (const auto& localDof : localDofs(fvGeometry))
        {
            const auto value = elemVars[localDof].priVar(0);
            u += ipCache.shapeValues()[localDof.index()][0]*value;
            gradU.axpy(value, ipCache.gradN(localDof.index()));
        }
    }
};

} // end namespace Dumux

#endif
