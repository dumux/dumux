// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Taylor-Green vortex test for the (hybrid) CVFE Navier-Stokes models
 *        (Taylor & Green 1937 \cite Taylor1937 in 2D, Antuono 2020 \cite Antuono2020 in 3D).
 */
#ifndef DUMUX_TAYLOR_GREEN_VORTEX_TEST_PROBLEM_HH
#define DUMUX_TAYLOR_GREEN_VORTEX_TEST_PROBLEM_HH

#include <cmath>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/math.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/discretization/cvfe/interpolate.hh>
#include <dumux/discretization/dirichletconstraints.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief Taylor-Green vortex test problem for the (hybrid) CVFE schemes.
 *
 * In 2D, the classic Taylor-Green vortex \cite Taylor1937 is considered, in 3D the tri-periodic
 * Beltrami flow of Antuono \cite Antuono2020. Both are exact solutions of the incompressible
 * Navier-Stokes equations which decay in time with the factor \f$ F(t) = \exp(-d \nu k^2 t) \f$,
 * where \f$ d \f$ is the dimension. For the stationary variant, \f$ F \equiv 1 \f$ and the
 * vortex is sustained by a manufactured source term. The analytical velocity is prescribed
 * as Dirichlet boundary condition. See README.md for details.
 */
template <class TypeTag, class BaseProblem>
class TaylorGreenTestProblem : public BaseProblem
{
    using ParentType = BaseProblem;

    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridDiscretization::LocalView;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Sources = typename ParentType::Sources;
    using DirichletValues = typename ParentType::DirichletValues;
    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ConstraintInfo = Dumux::DirichletConstraintInfo<ModelTraits::numEq()>;
    using ConstraintValues = Dune::FieldVector<Scalar, ModelTraits::numEq()>;
    using GridIndexType = typename IndexTraits<typename GridDiscretization::GridView>::GridIndex;
    using DirichletConstraintData = Dumux::DirichletConstraintData<ConstraintInfo, ConstraintValues, GridIndexType>;

    static constexpr int dimWorld = GridDiscretization::GridView::dimensionworld;
    static_assert(dimWorld == 2 || dimWorld == 3, "The Taylor-Green vortex test is only implemented in 2D and 3D");

    using Element = typename ElementDiscretization::Element;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Velocity = Dune::FieldVector<Scalar, dimWorld>;
    using VelocityGradient = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    using Indices = typename ModelTraits::Indices;

    TaylorGreenTestProblem(std::shared_ptr<const GridDiscretization> gridDiscretization,
                           std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridDiscretization, couplingManager)
    , time_(0.0)
    {
        isStationary_ = getParam<bool>("Problem.IsStationary");
        rho_ = getParam<Scalar>("Component.LiquidDensity");
        nu_ = getParam<Scalar>("Component.LiquidKinematicViscosity");
        mu_ = rho_*nu_;
        k_ = getParam<Scalar>("Problem.WaveNumber", 2.0*M_PI);
        u0_ = getParam<Scalar>("Problem.ReferenceVelocity", 1.0);
        useNeumann_ = getParam<bool>("Problem.UseNeumann", false);
        useUnsymmetrizedVelocityGradient_ = getParam<bool>("FreeFlow.EnableUnsymmetrizedVelocityGradient", false);

        updateConstraints_();
    }

    /*!
     * \brief Return the sources within the domain.
     *
     * The source is computed from the analytical solution, using that it is divergence free,
     * \f$ \mathbf{f} = \rho \partial_t \mathbf{u} + \rho (\mathbf{u} \cdot \nabla) \mathbf{u} - \mu \Delta \mathbf{u} + \nabla p \f$.
     * It vanishes for the instationary Navier-Stokes problem.
     *
     * \param globalPos The global position
     */
    Sources sourceAtPos(const GlobalPosition& globalPos) const
    {
        Sources source(0.0);
        if constexpr (ParentType::isMomentumProblem())
        {
            const auto u = velocity_(globalPos, time_);

            // time derivative and viscous term (the velocity is an eigenfunction of the Laplacian)
            const Scalar laplacianFactor = -dimWorld*k_*k_;
            const Scalar timeDerivativeFactor = isStationary_ ? 0.0 : laplacianFactor*nu_;
            for (int i = 0; i < dimWorld; ++i)
                source[i] = (rho_*timeDerivativeFactor - mu_*laplacianFactor)*u[i];

            source += pressureGradient_(globalPos, time_);

            // convective term (u.grad)u
            if (this->enableInertiaTerms())
            {
                Velocity convection(0.0);
                velocityGradient_(globalPos, time_).mv(u, convection);
                source.axpy(rho_, convection);
            }
        }

        return source;
    }

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on the boundary.
     *
     * \param globalPos The position at the boundary
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;

        if constexpr (ParentType::isMomentumProblem())
        {
            if (isMomentumFluxBoundary_(globalPos))
                values.setAllFluxBoundary();
        }
        else
            values.setFluxBoundary(Indices::conti0EqIdx);

        return values;
    }

    /*!
     * \brief Return Dirichlet boundary values at a given position
     *
     * \param globalPos The global position
     */
    DirichletValues dirichletAtPos(const GlobalPosition& globalPos) const
    { return analyticalSolution(globalPos, time_); }

    /*!
     * \brief Return Dirichlet boundary constraints and internal constraints.
     */
    const auto& constraints() const
    { return constraints_; }

    /*!
     * \brief Evaluates the boundary flux related to a localDof at a given interpolation point.
     *
     * \param elemDisc The element discretization
     * \param elemVars All variables related to the element
     * \param faceIpData Interpolation point data
     */
    template<class ElementVariables, class FaceIpData>
    BoundaryFluxes boundaryFlux(const ElementDiscretization& elemDisc,
                                const ElementVariables& elemVars,
                                const FaceIpData& faceIpData) const
    {
        BoundaryFluxes values(0.0);

        if constexpr (ParentType::isMomentumProblem())
        {
            const auto& globalPos = faceIpData.global();
            const auto u = velocity_(globalPos, time_);
            const auto gradU = velocityGradient_(globalPos, time_);

            VelocityGradient momentumFlux(0.0);
            for (int i = 0; i < dimWorld; ++i)
            {
                for (int j = 0; j < dimWorld; ++j)
                {
                    momentumFlux[i][j] = useUnsymmetrizedVelocityGradient_ ? -mu_*gradU[i][j] : -mu_*(gradU[i][j] + gradU[j][i]);
                    if (this->enableInertiaTerms())
                        momentumFlux[i][j] += rho_*u[i]*u[j];
                }
                momentumFlux[i][i] += pressure_(globalPos, time_);
            }

            momentumFlux.mv(faceIpData.unitOuterNormal(), values);
        }
        else
        {
            const auto& scvf = elemDisc.scvf(faceIpData.scvfIndex());
            const auto insideDensity = elemVars[scvf.insideScvIdx()].density();
            values[Indices::conti0EqIdx] = this->velocity(elemDisc, faceIpData) * insideDensity * scvf.unitOuterNormal();
        }

        return values;
    }

    /*!
     * \brief Evaluates the boundary flux related to a localDof at a given interpolation point.
     *
     * \param elemDisc The element discretization
     * \param elemVars All variables related to the element
     * \param elemFluxVarsCache The element flux variables cache
     * \param faceIpData Interpolation point data
     */
    template<class ElementVariables, class ElementFluxVariablesCache, class FaceIpData>
    BoundaryFluxes boundaryFlux(const ElementDiscretization& elemDisc,
                                const ElementVariables& elemVars,
                                const ElementFluxVariablesCache& elemFluxVarsCache,
                                const FaceIpData& faceIpData) const
    { return boundaryFlux(elemDisc, elemVars, faceIpData); }

    // \}

    /*!
     * \brief Return the analytical solution of the problem at a given time and position
     *
     * \param globalPos The global position
     * \param time The time
     */
    DirichletValues analyticalSolution(const GlobalPosition& globalPos, const Scalar time) const
    {
        DirichletValues values(0.0);

        if constexpr (ParentType::isMomentumProblem())
        {
            const auto u = velocity_(globalPos, time);
            for (int i = 0; i < dimWorld; ++i)
                values[Indices::velocity(i)] = u[i];
        }
        else
            values[Indices::pressureIdx] = pressure_(globalPos, time);

        return values;
    }

    /*!
     * \brief Return the analytical solution of the problem at a given position and the current time
     *
     * \param globalPos The global position
     */
    DirichletValues analyticalSolution(const GlobalPosition& globalPos) const
    { return analyticalSolution(globalPos, time_); }

    /*!
     * \brief Return the gradient of the analytical solution at a given position and the current time
     *
     * \param globalPos The global position
     */
    Dune::FieldVector<GlobalPosition, DirichletValues::dimension> gradAnalyticalSolution(const GlobalPosition& globalPos) const
    {
        Dune::FieldVector<GlobalPosition, DirichletValues::dimension> values;

        if constexpr (ParentType::isMomentumProblem())
        {
            const auto gradU = velocityGradient_(globalPos, time_);
            for (int i = 0; i < dimWorld; ++i)
                values[Indices::velocity(i)] = gradU[i];
        }
        else
            values[Indices::pressureIdx] = pressureGradient_(globalPos, time_);

        return values;
    }

    /*!
     * \brief Applies the initial solution for all degrees of freedom
     * \note The stationary problem starts from zero, the instationary one from
     *       the (projected) analytical solution at the current time.
     */
    template<class SolutionVector>
    void applyInitialSolution(SolutionVector& sol) const
    {
        sol.resize(this->gridDiscretization().numDofs());
        sol = 0.0;

        if (isStationary_)
            return;

        const auto initialSolution = [&](const GlobalPosition& globalPos)
        { return typename SolutionVector::value_type(analyticalSolution(globalPos, time_)); };

        // the L2 projection also yields the correct coefficients for enriched (e.g. bubble) basis functions
        if constexpr (ParentType::isMomentumProblem())
            CVFE::interpolate(this->gridDiscretization(), sol, initialSolution, CVFE::InterpolationPolicy::L2Projection{});
        else
            CVFE::interpolate(this->gridDiscretization(), sol, initialSolution);
    }

    /*!
     * \brief Updates the time and the time-dependent Dirichlet constraints
     */
    void updateTime(const Scalar time)
    { setTime(time); }

    /*!
     * \brief Updates the time and the time-dependent Dirichlet constraints
     * \note This is called by the multi-stage assembler (on a const problem) for each stage.
     */
    void setTime(const Scalar time) const
    {
        time_ = time;
        updateConstraints_();
    }

    //! The current time
    Scalar time() const
    { return time_; }

private:
    //! The momentum flux is prescribed at the upper boundaries (x_i = max) if Neumann conditions are enabled
    bool isMomentumFluxBoundary_(const GlobalPosition& globalPos) const
    {
        if (!useNeumann_)
            return false;

        static constexpr Scalar eps = 1e-6;
        const auto& bBoxMax = this->gridDiscretization().bBoxMax();
        for (int i = 0; i < dimWorld; ++i)
            if (globalPos[i] > bBoxMax[i] - eps)
                return true;

        return false;
    }

    void updateConstraints_() const
    {
        constraints_.clear();
        if constexpr (ParentType::isMomentumProblem())
            appendDirichletConstraints_();
        else if (!useNeumann_)
            appendPressureConstraint_();
    }

    void appendDirichletConstraints_() const
    {
        auto elemDisc = localView(this->gridDiscretization());
        for (const auto& element : elements(this->gridDiscretization().gridView()))
        {
            elemDisc.bind(element);
            for (const auto& boundaryFace : boundaryFaces(elemDisc))
            {
                if (isMomentumFluxBoundary_(boundaryFace.center()))
                    continue;

                for (const auto& localDof : localDofs(elemDisc, boundaryFace))
                {
                    ConstraintInfo info;
                    info.setAll();
                    ConstraintValues values(dirichletAtPos(ipData(elemDisc, localDof).global()));
                    constraints_.push_back(DirichletConstraintData{std::move(info), std::move(values), localDof.dofIndex()});
                }
            }
        }
    }

    //! Without flux boundaries, the pressure is only determined up to a constant, so fix it at one dof
    void appendPressureConstraint_() const
    {
        static_assert(GridDiscretization::discMethod == DiscretizationMethods::box,
                      "The pressure constraint is only implemented for the Box mass discretization.");

        static constexpr Scalar eps = 1e-8;
        auto elemDisc = localView(this->gridDiscretization());
        for (const auto& element : elements(this->gridDiscretization().gridView()))
        {
            elemDisc.bind(element);
            for (const auto& scv : scvs(elemDisc))
            {
                if ((scv.dofPosition() - this->gridDiscretization().bBoxMin()).two_norm() < eps)
                {
                    ConstraintInfo info;
                    info.setAll();
                    ConstraintValues values(pressure_(scv.dofPosition(), time_));
                    constraints_.push_back(DirichletConstraintData{std::move(info), std::move(values), scv.dofIndex()});
                    return;
                }
            }
        }
    }

    //! The temporal decay factor of the velocity
    Scalar decayFactor_(const Scalar t) const
    {
        using std::exp;
        return isStationary_ ? 1.0 : exp(-dimWorld*nu_*k_*k_*t);
    }

    //! The velocity (without decay factor)
    Velocity velocityShape_(const GlobalPosition& globalPos) const
    {
        using std::sin; using std::cos;
        Velocity u(0.0);

        if constexpr (dimWorld == 2)
        {
            const Scalar kx = k_*globalPos[0];
            const Scalar ky = k_*globalPos[1];
            u[0] = u0_*sin(kx)*cos(ky);
            u[1] = -u0_*cos(kx)*sin(ky);
        }
        else
        {
            // cyclic permutations (x, y, z) -> (y, z, x) -> (z, x, y)
            for (int i = 0; i < 3; ++i)
            {
                const Scalar ka = k_*globalPos[i];
                const Scalar kb = k_*globalPos[(i+1)%3];
                const Scalar kc = k_*globalPos[(i+2)%3];
                u[i] = antuonoScale_()*(sin(ka - alpha_)*cos(kb - beta_)*sin(kc)
                                        - cos(kc - alpha_)*sin(ka - beta_)*sin(kb));
            }
        }

        return u;
    }

    //! The velocity gradient (without decay factor), gradU[i][j] = du_i/dx_j
    VelocityGradient velocityGradientShape_(const GlobalPosition& globalPos) const
    {
        using std::sin; using std::cos;
        VelocityGradient gradU(0.0);

        if constexpr (dimWorld == 2)
        {
            const Scalar kx = k_*globalPos[0];
            const Scalar ky = k_*globalPos[1];
            gradU[0][0] = u0_*k_*cos(kx)*cos(ky);
            gradU[0][1] = -u0_*k_*sin(kx)*sin(ky);
            gradU[1][0] = u0_*k_*sin(kx)*sin(ky);
            gradU[1][1] = -u0_*k_*cos(kx)*cos(ky);
        }
        else
        {
            for (int i = 0; i < 3; ++i)
            {
                const int a = i, b = (i+1)%3, c = (i+2)%3;
                const Scalar ka = k_*globalPos[a];
                const Scalar kb = k_*globalPos[b];
                const Scalar kc = k_*globalPos[c];
                const Scalar scale = antuonoScale_()*k_;
                gradU[i][a] = scale*(cos(ka - alpha_)*cos(kb - beta_)*sin(kc)
                                     - cos(kc - alpha_)*cos(ka - beta_)*sin(kb));
                gradU[i][b] = scale*(-sin(ka - alpha_)*sin(kb - beta_)*sin(kc)
                                     - cos(kc - alpha_)*sin(ka - beta_)*cos(kb));
                gradU[i][c] = scale*(sin(ka - alpha_)*cos(kb - beta_)*cos(kc)
                                     + sin(kc - alpha_)*sin(ka - beta_)*sin(kb));
            }
        }

        return gradU;
    }

    Velocity velocity_(const GlobalPosition& globalPos, const Scalar t) const
    { return velocityShape_(globalPos) *= decayFactor_(t); }

    VelocityGradient velocityGradient_(const GlobalPosition& globalPos, const Scalar t) const
    { return velocityGradientShape_(globalPos) *= decayFactor_(t); }

    Scalar pressure_(const GlobalPosition& globalPos, const Scalar t) const
    {
        using std::cos;
        const Scalar f = decayFactor_(t);

        if constexpr (dimWorld == 2)
            return 0.25*rho_*u0_*u0_*(cos(2.0*k_*globalPos[0]) + cos(2.0*k_*globalPos[1]))*f*f;
        else
            // Bernoulli pressure p = -rho/2 |u|^2 of the Beltrami flow
            return -0.5*rho_*velocity_(globalPos, t).two_norm2();
    }

    GlobalPosition pressureGradient_(const GlobalPosition& globalPos, const Scalar t) const
    {
        using std::sin;
        GlobalPosition gradP(0.0);

        if constexpr (dimWorld == 2)
        {
            const Scalar f = decayFactor_(t);
            for (int i = 0; i < 2; ++i)
                gradP[i] = -0.5*rho_*u0_*u0_*k_*sin(2.0*k_*globalPos[i])*f*f;
        }
        else
        {
            // grad p = -rho (grad u)^T u
            velocityGradient_(globalPos, t).mtv(velocity_(globalPos, t), gradP);
            gradP *= -rho_;
        }

        return gradP;
    }

    //! Normalization of the 3D solution (such that the mean kinetic energy is rho U_0^2/2)
    static Scalar antuonoScale_()
    { using std::sqrt; return 4.0*sqrt(2.0)/(3.0*sqrt(3.0)); }

    static constexpr Scalar alpha_ = 5.0*M_PI/6.0;
    static constexpr Scalar beta_ = M_PI/6.0;

    Scalar rho_;
    Scalar nu_;
    Scalar mu_;
    Scalar k_;
    Scalar u0_;
    mutable Scalar time_;
    bool isStationary_;
    bool useNeumann_;
    bool useUnsymmetrizedVelocityGradient_;
    mutable std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux

#endif
