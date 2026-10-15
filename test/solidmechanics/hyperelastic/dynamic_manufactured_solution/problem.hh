// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Plane-strain manufactured solution for the dynamic hyperelastic model with a viscous
 *        stress, spring-dashpot boundary conditions and a follower pressure load.
 *
 * On the unit square, the displacement
 * \f[ \mathbf{u} = A t^2 \mathbf{g}, \quad
 *     \mathbf{g} = \begin{pmatrix} \sin x \cos y \\ \cos x \cos y \end{pmatrix}, \f]
 * with \f$ A = 0.3 \f$ has a constant acceleration, which the average-acceleration Newmark
 * method integrates exactly. The stress is compressible neo-Hookean with a viscous part,
 * \f[ \mathbf{P} = \mu \mathbf{F} + (\lambda \ln J - \mu) \mathbf{F}^{-\top} + \eta \mathbf{F}\mathbf{D},
 *     \quad \mathbf{D} = \operatorname{sym}(\mathbf{F}^\top \dot{\mathbf{F}}), \f]
 * and the body force is \f$ \varrho \ddot{\mathbf{u}} - \nabla \cdot \mathbf{P} \f$.
 *
 * On the whole boundary, the residual contains
 * \f[ \int \left(\alpha \mathbf{u} + \beta \dot{\mathbf{u}} + p J \mathbf{F}^{-\top}\mathbf{N}
 *     + \mathbf{h}\right) \cdot \mathbf{w}, \f]
 * with \f$ \mathbf{h} \f$ chosen such that the exact solution satisfies the boundary condition.
 */
#ifndef DUMUX_TEST_HYPERELASTIC_DYNAMIC_MMS_PROBLEM_HH
#define DUMUX_TEST_HYPERELASTIC_DYNAMIC_MMS_PROBLEM_HH

#include <array>
#include <cmath>
#include <memory>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/problemwithspatialparams.hh>
#include <dumux/common/properties.hh>
#include <dumux/experimental/timestepping/newmarkbeta.hh>

namespace Dumux {

template<class TypeTag>
class HyperelasticDynamicMMSProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridGeometry::LocalView;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = typename ParentType::Traits::NumEqVector;
    using BoundaryTypes = Dumux::Experimental::BoundaryTypes<PrimaryVariables::size()>;
    using Newmark = Experimental::NewmarkBeta<Scalar, SolutionVector>;

    static constexpr int dim = 2;
    using Tensor = Dune::FieldMatrix<Scalar, dim, dim>;
    using Vector = Dune::FieldVector<Scalar, dim>;

public:
    HyperelasticDynamicMMSProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        mu_ = getParam<Scalar>("Problem.ShearModulus");
        lambda_ = getParam<Scalar>("Problem.LameLambda");
        eta_ = getParam<Scalar>("Problem.Viscosity");
        rho_ = getParam<Scalar>("Problem.Density");
        springStiffness_ = getParam<Scalar>("Problem.SpringStiffness");
        dashpotCoefficient_ = getParam<Scalar>("Problem.DashpotCoefficient");
        pressure_ = getParam<Scalar>("Problem.Pressure");
    }

    void setNewmarkScheme(std::shared_ptr<const Newmark> newmark)
    { newmark_ = std::move(newmark); }

    void setTime(const Scalar time)
    { time_ = time; }

    void setTimeStepSize(const Scalar dt)
    { dt_ = dt; }

    Scalar time() const
    { return time_; }

    Vector exactDisplacement(const GlobalPosition& pos, const Scalar t) const
    { return amplitude_*t*t*g_(pos); }

    Tensor exactDisplacementGradient(const GlobalPosition& pos, const Scalar t) const
    { return amplitude_*t*t*gradG_(pos); }

    Vector exactAcceleration(const GlobalPosition& pos) const
    { return 2.0*amplitude_*g_(pos); }

    template<class ElementVariables, class IpData>
    NumEqVector source(const ElementDiscretization&, const ElementVariables&, const IpData& ipData) const
    {
        const auto& pos = ipData.global();
        const auto f = rho_*exactAcceleration(pos) - exactStressDivergence_(pos);
        return {f[0], f[1]};
    }

    template<class IpData>
    Tensor firstPiolaKirchhoffStressTensor(const Tensor& F, const Tensor& Fdot,
                                           const ElementDiscretization&, const IpData&) const
    { return stress_(F, Fdot); }

    template<class IpData>
    Scalar solidDensity(const ElementDiscretization&, const IpData&) const
    { return rho_; }

    template<class ElementVariables, class LocalDof>
    Vector acceleration(const ElementDiscretization&, const ElementVariables& elemVars, const LocalDof& localDof) const
    { return newmark_->acceleration(localDof.dofIndex(), dt_, elemVars[localDof].displacement()); }

    template<class ElementVariables, class LocalDof>
    Vector velocity(const ElementDiscretization&, const ElementVariables& elemVars, const LocalDof& localDof) const
    { return newmark_->velocity(localDof.dofIndex(), dt_, elemVars[localDof].displacement()); }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition&) const
    {
        BoundaryTypes values;
        values.setAllFluxBoundary();
        return values;
    }

    template<class ElementVariables, class FaceIpData>
    NumEqVector boundaryFlux(const ElementDiscretization& elemDisc,
                             const ElementVariables& elemVars,
                             const FaceIpData& ipData) const
    {
        const auto& ipCache = cache(elemVars, ipData);
        Vector u(0.0), v(0.0);
        Tensor F(0.0);
        for (const auto& localDof : localDofs(elemDisc))
        {
            const auto& displacement = elemVars[localDof].displacement();
            const auto shapeValue = ipCache.shapeValues()[localDof.index()][0];
            u.axpy(shapeValue, displacement);
            v.axpy(shapeValue, velocity(elemDisc, elemVars, localDof));
            for (int i = 0; i < dim; ++i)
                F[i].axpy(displacement[i], ipCache.gradN(localDof.index()));
        }
        for (int i = 0; i < dim; ++i)
            F[i][i] += 1.0;
        const auto& N = ipData.unitOuterNormal();

        const auto& pos = ipData.global();
        const auto FExact = exactDeformationGradient_(pos);
        auto flux = springDashpotPressure_(u, v, F, N);
        flux -= springDashpotPressure_(exactDisplacement(pos, time_), exactVelocity_(pos), FExact, N);
        Vector PN;
        stress_(FExact, exactVelocityGradient_(pos)).mv(N, PN);
        flux -= PN;
        return {flux[0], flux[1]};
    }

private:
    Tensor stress_(const Tensor& F, const Tensor& Fdot) const
    {
        using std::log;
        auto FinvT = F;
        FinvT.invert();
        FinvT = getTransposed(FinvT);
        auto D = getTransposed(F)*Fdot;
        D = 0.5*(D + getTransposed(D));
        return mu_*F + (lambda_*log(F.determinant()) - mu_)*FinvT + eta_*F*D;
    }

    /*!
     * \brief The divergence \f$ \sum_k \partial_k \mathbf{P} \mathbf{e}_k \f$ of the stress of the
     *        exact solution, with
     * \f[ \partial_k \mathbf{P} = \mu \partial_k \mathbf{F}
     *     - (\lambda \ln J - \mu) \mathbf{F}^{-\top} \partial_k \mathbf{F}^\top \mathbf{F}^{-\top}
     *     + \lambda \operatorname{tr}(\mathbf{F}^{-1} \partial_k \mathbf{F}) \mathbf{F}^{-\top}
     *     + \eta (\partial_k \mathbf{F} \mathbf{D} + \mathbf{F} \partial_k \mathbf{D}) \f]
     */
    Vector exactStressDivergence_(const GlobalPosition& pos) const
    {
        using std::log;
        const auto F = exactDeformationGradient_(pos);
        const auto Fdot = exactVelocityGradient_(pos);
        auto Finv = F;
        Finv.invert();
        const auto FinvT = getTransposed(Finv);
        auto D = getTransposed(F)*Fdot;
        D = 0.5*(D + getTransposed(D));
        const auto lnJ = log(F.determinant());

        const auto gradGradG = gradGradG_(pos);
        Vector divP(0.0);
        for (int k = 0; k < dim; ++k)
        {
            const auto dF = amplitude_*time_*time_*gradGradG[k];
            const auto dFdot = 2.0*amplitude_*time_*gradGradG[k];
            auto dD = getTransposed(dF)*Fdot + getTransposed(F)*dFdot;
            dD = 0.5*(dD + getTransposed(dD));
            const auto dP = mu_*dF - (lambda_*lnJ - mu_)*FinvT*getTransposed(dF)*FinvT
                            + lambda_*trace(Finv*dF)*FinvT + eta_*(dF*D + F*dD);
            for (int i = 0; i < dim; ++i)
                divP[i] += dP[i][k];
        }
        return divP;
    }

    Vector exactVelocity_(const GlobalPosition& pos) const
    { return 2.0*amplitude_*time_*g_(pos); }

    Tensor exactDeformationGradient_(const GlobalPosition& pos) const
    {
        auto F = exactDisplacementGradient(pos, time_);
        for (int i = 0; i < dim; ++i)
            F[i][i] += 1.0;
        return F;
    }

    Tensor exactVelocityGradient_(const GlobalPosition& pos) const
    { return 2.0*amplitude_*time_*gradG_(pos); }

    static Vector g_(const GlobalPosition& pos)
    {
        using std::sin; using std::cos;
        const auto x = pos[0], y = pos[1];
        return {sin(x)*cos(y), cos(x)*cos(y)};
    }

    static Tensor gradG_(const GlobalPosition& pos)
    {
        using std::sin; using std::cos;
        const auto x = pos[0], y = pos[1];
        return {{cos(x)*cos(y), -sin(x)*sin(y)},
                {-sin(x)*cos(y), -cos(x)*sin(y)}};
    }

    //! The derivatives \f$ \partial_k \nabla \mathbf{g} \f$ in the directions \f$ k = x, y \f$
    static std::array<Tensor, dim> gradGradG_(const GlobalPosition& pos)
    {
        using std::sin; using std::cos;
        const auto x = pos[0], y = pos[1];
        return {{{{-sin(x)*cos(y), -cos(x)*sin(y)},
                  {-cos(x)*cos(y), sin(x)*sin(y)}},
                 {{-cos(x)*sin(y), -sin(x)*cos(y)},
                  {sin(x)*sin(y), -cos(x)*cos(y)}}}};
    }

    Vector springDashpotPressure_(const Vector& u, const Vector& v, const Tensor& F, const Vector& N) const
    {
        auto FinvT = F;
        FinvT.invert();
        FinvT = getTransposed(FinvT);
        Vector value;
        FinvT.mv(N, value);
        value *= pressure_*F.determinant();
        value.axpy(springStiffness_, u);
        value.axpy(dashpotCoefficient_, v);
        return value;
    }

    static constexpr Scalar amplitude_ = 0.3;
    Scalar mu_, lambda_, eta_, rho_;
    Scalar springStiffness_, dashpotCoefficient_, pressure_;
    Scalar time_ = 0.0;
    Scalar dt_ = 1.0;
    std::shared_ptr<const Newmark> newmark_;
};

} // end namespace Dumux

#endif
