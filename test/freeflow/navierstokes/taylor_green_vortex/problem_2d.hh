// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief The two-dimensional Taylor-Green vortex (Taylor 1923 \cite Taylor1923Decay), see README.md.
 */
#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_TAYLOR_GREEN_VORTEX_PROBLEM_2D_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_TAYLOR_GREEN_VORTEX_PROBLEM_2D_HH

#include <cmath>
#include <numbers>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

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
 * \brief The two-dimensional Taylor-Green vortex for the (hybrid) CVFE schemes
 *
 * The bi-periodic solution of Taylor \cite Taylor1923Decay (in the form of Antuono \cite Antuono2020, Sec. 2.1)
 * \f[
 * u = -U_0 \cos(kx) \sin(ky) F(t), \quad
 * v = U_0 \sin(kx) \cos(ky) F(t), \quad
 * p = -\frac{\rho U_0^2}{4} \left( \cos(2kx) + \cos(2ky) \right) F^2(t),
 * \f]
 * with the wave number \f$ k = 2\pi/L \f$, the wavelength \f$ L \f$ (Problem.WaveLength),
 * the amplitude \f$ U_0 \f$ (Problem.ReferenceVelocity) and the decay factor \f$ F(t) = e^{-2 \nu k^2 t} \f$,
 * is an exact solution of the instationary Navier-Stokes equations without source term.
 * The analytical solution is used for the initial condition, the Dirichlet velocity on the entire boundary,
 * the pressure at one degree of freedom and the source terms of the stationary and Stokes variants.
 */
template <class TypeTag, class BaseProblem>
class TaylorGreenVortexProblem2D : public BaseProblem
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

    static constexpr int dim = 2;
    static_assert(GridDiscretization::GridView::dimensionworld == dim, "The 2D Taylor-Green vortex requires a 2D grid");

    using Element = typename ElementDiscretization::Element;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Velocity = Dune::FieldVector<Scalar, dim>;
    using VelocityGradient = Dune::FieldMatrix<Scalar, dim, dim>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    using Indices = typename ModelTraits::Indices;
    //! The vorticity of a planar flow is a scalar
    using Vorticity = Scalar;

    TaylorGreenVortexProblem2D(std::shared_ptr<const GridDiscretization> gridDiscretization,
                               std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridDiscretization, couplingManager)
    , time_(0.0)
    {
        isStationary_ = getParam<bool>("Problem.IsStationary");
        rho_ = getParam<Scalar>("Component.LiquidDensity");
        nu_ = getParam<Scalar>("Component.LiquidKinematicViscosity");
        mu_ = rho_*nu_;
        // by default, the domain contains exactly one wavelength
        const Scalar domainLength = gridDiscretization->bBoxMax()[0] - gridDiscretization->bBoxMin()[0];
        k_ = 2.0*std::numbers::pi_v<Scalar>/getParam<Scalar>("Problem.WaveLength", domainLength);
        u0_ = getParam<Scalar>("Problem.ReferenceVelocity", 1.0);

        updateConstraints_();
    }

    /*!
     * \brief Return the sources within the domain.
     *
     * The instationary Navier-Stokes problem has no source term. The stationary variant (the vortex at t = 0) adds
     * the storage term \f$ -\rho \partial_t \mathbf{u} = 2 \mu k^2 \mathbf{u} \f$, the Stokes variant the
     * convective term \f$ -\nabla \cdot (\rho \mathbf{u} \otimes \mathbf{u}) = \nabla p \f$ (both identities hold
     * for the analytical solution, see README.md).
     *
     * \param globalPos The global position
     */
    Sources sourceAtPos(const GlobalPosition& globalPos) const
    {
        Sources source(0.0);
        if constexpr (ParentType::isMomentumProblem())
        {
            const Scalar t = solutionTime_(time_);
            if (isStationary_)
            {
                const auto u = velocity_(globalPos, t);
                for (int i = 0; i < dim; ++i)
                    source[Indices::momentumBalanceIdx(i)] += 2.0*mu_*k_*k_*u[i];
            }

            if (!this->enableInertiaTerms())
            {
                const auto gradP = pressureGradient_(globalPos, t);
                for (int i = 0; i < dim; ++i)
                    source[Indices::momentumBalanceIdx(i)] += gradP[i];
            }
        }

        return source;
    }

    /*!
     * \brief The velocity is prescribed on the entire boundary (by constraints),
     *        the mass balance has a flux boundary
     *
     * \param globalPos The position at the boundary
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        if constexpr (!ParentType::isMomentumProblem())
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
     * \brief The Dirichlet constraints of the velocity on the boundary and of the pressure at one dof
     */
    const auto& constraints() const
    { return constraints_; }

    /*!
     * \brief The boundary flux: the mass flux of the discrete velocity
     *        (the momentum balance has no flux boundaries)
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
        if constexpr (!ParentType::isMomentumProblem())
        {
            const auto& scvf = elemDisc.scvf(faceIpData.scvfIndex());
            const auto insideDensity = elemVars[scvf.insideScvIdx()].density();
            values[Indices::conti0EqIdx] = this->velocity(elemDisc, faceIpData) * insideDensity * scvf.unitOuterNormal();
        }
        return values;
    }

    /*!
     * \brief The boundary flux (overload with flux variables cache)
     */
    template<class ElementVariables, class ElementFluxVariablesCache, class FaceIpData>
    BoundaryFluxes boundaryFlux(const ElementDiscretization& elemDisc,
                                const ElementVariables& elemVars,
                                const ElementFluxVariablesCache& elemFluxVarsCache,
                                const FaceIpData& faceIpData) const
    { return boundaryFlux(elemDisc, elemVars, faceIpData); }

    /*!
     * \brief Return the analytical solution of the problem at a given time and position
     *
     * \param globalPos The global position
     * \param time The time
     */
    DirichletValues analyticalSolution(const GlobalPosition& globalPos, const Scalar time) const
    {
        DirichletValues values(0.0);
        const Scalar t = solutionTime_(time);

        if constexpr (ParentType::isMomentumProblem())
        {
            const auto u = velocity_(globalPos, t);
            for (int i = 0; i < dim; ++i)
                values[Indices::velocity(i)] = u[i];
        }
        else
            values[Indices::pressureIdx] = pressure_(globalPos, t);

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
        const Scalar t = solutionTime_(time_);

        if constexpr (ParentType::isMomentumProblem())
        {
            const auto gradU = velocityGradient_(globalPos, t);
            for (int i = 0; i < dim; ++i)
                values[Indices::velocity(i)] = gradU[i];
        }
        else
            values[Indices::pressureIdx] = pressureGradient_(globalPos, t);

        return values;
    }

    /*!
     * \brief Return the analytical vorticity \f$ \omega = \partial_x v - \partial_y u = 2 k U_0 \cos(kx) \cos(ky) F(t) \f$
     *
     * \param globalPos The global position
     * \param time The time
     */
    Vorticity analyticalVorticity(const GlobalPosition& globalPos, const Scalar time) const
    {
        using std::cos;
        return 2.0*k_*u0_*cos(k_*globalPos[0])*cos(k_*globalPos[1])*decayFactor_(solutionTime_(time));
    }

    /*!
     * \brief Return the analytical vorticity at a given position and the current time
     *
     * \param globalPos The global position
     */
    Vorticity analyticalVorticity(const GlobalPosition& globalPos) const
    { return analyticalVorticity(globalPos, time_); }

    /*!
     * \brief Applies the (L2-projected) analytical solution at the current time, which is the initial condition
     *        of the instationary problem and the initial guess of the Newton solver for the stationary problem
     */
    template<class SolutionVector>
    void applyInitialSolution(SolutionVector& sol) const
    {
        sol.resize(this->gridDiscretization().numDofs());
        sol = 0.0;

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
    //! The stationary problem is the vortex at t = 0
    Scalar solutionTime_(const Scalar time) const
    { return isStationary_ ? 0.0 : time; }

    void updateConstraints_() const
    {
        constraints_.clear();
        if constexpr (ParentType::isMomentumProblem())
            appendVelocityConstraints_();
        else
            appendPressureConstraint_();
    }

    //! The analytical velocity on all boundary dofs
    void appendVelocityConstraints_() const
    {
        appendBoundaryDirichletConstraints(
            constraints_, this->gridDiscretization(),
            [] (const auto& elemDisc, const auto& boundaryFace, const auto& localDof)
            {
                ConstraintInfo info;
                info.setAll();
                return info;
            },
            [&] (const auto& elemDisc, const auto& localDof)
            { return dirichletAtPos(ipData(elemDisc, localDof).global()); }
        );
    }

    //! With the velocity prescribed on the entire boundary, the pressure is only determined up to a constant,
    //! so it is fixed to the analytical value at the dof of the lower left corner
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
                    ConstraintValues values(pressure_(scv.dofPosition(), solutionTime_(time_)));
                    constraints_.push_back(DirichletConstraintData{std::move(info), std::move(values), scv.dofIndex()});
                    return;
                }
            }
        }
    }

    //! The temporal decay factor F(t) = exp(-2 nu k^2 t)
    Scalar decayFactor_(Scalar t) const
    { using std::exp; return exp(-2.0*nu_*k_*k_*t); }

    Velocity velocity_(const GlobalPosition& x, Scalar t) const
    {
        using std::sin; using std::cos;
        const Scalar kx = k_*x[0], ky = k_*x[1];
        const Scalar a = u0_*decayFactor_(t);
        return { -a*cos(kx)*sin(ky), a*sin(kx)*cos(ky) };
    }

    //! The velocity gradient, gradU[i][j] = du_i/dx_j
    VelocityGradient velocityGradient_(const GlobalPosition& x, Scalar t) const
    {
        using std::sin; using std::cos;
        const Scalar kx = k_*x[0], ky = k_*x[1];
        const Scalar a = u0_*k_*decayFactor_(t);
        VelocityGradient gradU;
        gradU[0][0] = a*sin(kx)*sin(ky);
        gradU[0][1] = -a*cos(kx)*cos(ky);
        gradU[1][0] = a*cos(kx)*cos(ky);
        gradU[1][1] = -a*sin(kx)*sin(ky);
        return gradU;
    }

    Scalar pressure_(const GlobalPosition& x, Scalar t) const
    {
        using std::cos;
        const Scalar f = decayFactor_(t);
        return -0.25*rho_*u0_*u0_*(cos(2.0*k_*x[0]) + cos(2.0*k_*x[1]))*f*f;
    }

    GlobalPosition pressureGradient_(const GlobalPosition& x, Scalar t) const
    {
        using std::sin;
        const Scalar f = decayFactor_(t);
        const Scalar a = 0.5*rho_*u0_*u0_*k_*f*f;
        return { a*sin(2.0*k_*x[0]), a*sin(2.0*k_*x[1]) };
    }

    Scalar rho_, nu_, mu_;
    Scalar k_, u0_;
    bool isStationary_;
    mutable Scalar time_;
    mutable std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux

#endif
