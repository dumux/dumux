// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Test for the (hybrid) CVFE Navier-Stokes models with analytical solution.
 */
#ifndef DUMUX_SINCOS_TEST_PROBLEM_NEWINTERFACE_HH
#define DUMUX_SINCOS_TEST_PROBLEM_NEWINTERFACE_HH

#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/discretization/dirichletconstraints.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief Test problem for the (hybrid) CVFE schemes.
 *
 * The 2D, incompressible Navier-Stokes equations for zero gravity and a Newtonian
 * flow is solved and compared to an analytical solution (sums/products of trigonometric functions).
 * For the instationary case, the velocities and pressures are periodical in time. The Dirichlet boundary conditions are
 * consistent with the analytical solution and in the instationary case time-dependent.
 */
template <class TypeTag, class BaseProblem>
class SincosTestProblemNewInterface : public BaseProblem
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

    static constexpr auto dimWorld = GridDiscretization::GridView::dimensionworld;
    using Element = typename ElementDiscretization::Element;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    using Indices = typename ModelTraits::Indices;

    SincosTestProblemNewInterface(std::shared_ptr<const GridDiscretization> gridDiscretization,
                                  std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridDiscretization, couplingManager)
    { init_(); }

    /*!
     * \brief Constructor for the momentum problem without coupling to a mass problem
     * \note Pressure, density and viscosity are then taken from the analytical solution
     *       and the fluid parameters, see pressureAtPos() and densityAtPos().
     */
    explicit SincosTestProblemNewInterface(std::shared_ptr<const GridDiscretization> gridDiscretization)
    : ParentType(gridDiscretization)
    { init_(); }

    /*!
     * \brief Return the sources within the domain.
     *
     * \param globalPos The global position
     */
    Sources sourceAtPos(const GlobalPosition& globalPos) const
    {
        Sources source(0.0);
        if constexpr (ParentType::isMomentumProblem())
        {
            const Scalar x = globalPos[0];
            const Scalar y = globalPos[1];
            const Scalar t = time_;

            source[Indices::momentumXBalanceIdx] = rho_*dtU_(x,y,t) - 2.0*mu_*dxxU_(x,y,t) - mu_*dyyU_(x,y,t) - mu_*dxyV_(x,y,t) + dxP_(x,y,t);
            source[Indices::momentumYBalanceIdx] = rho_*dtV_(x,y,t) - 2.0*mu_*dyyV_(x,y,t) - mu_*dxyU_(x,y,t) - mu_*dxxV_(x,y,t) + dyP_(x,y,t);

            if (this->enableInertiaTerms())
            {
                source[Indices::momentumXBalanceIdx] += rho_*dxUU_(x,y,t) + rho_*dyUV_(x,y,t);
                source[Indices::momentumYBalanceIdx] += rho_*dxUV_(x,y,t) + rho_*dyVV_(x,y,t);
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
            const Scalar x = faceIpData.global()[0];
            const Scalar y = faceIpData.global()[1];
            const Scalar t = time_;

            Dune::FieldMatrix<Scalar, dimWorld, dimWorld> momentumFlux(0.0);
            momentumFlux[0][0] = -2.0*mu_*dxU_(x,y,t) + p_(x,y,t);
            momentumFlux[0][1] = -mu_*(dyU_(x,y,t) + dxV_(x,y,t));
            momentumFlux[1][0] = momentumFlux[0][1];
            momentumFlux[1][1] = -2.0*mu_*dyV_(x,y,t) + p_(x,y,t);

            if (this->enableInertiaTerms())
            {
                momentumFlux[0][0] += rho_*u_(x,y,t)*u_(x,y,t);
                momentumFlux[0][1] += rho_*u_(x,y,t)*v_(x,y,t);
                momentumFlux[1][0] += rho_*v_(x,y,t)*u_(x,y,t);
                momentumFlux[1][1] += rho_*v_(x,y,t)*v_(x,y,t);
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
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        DirichletValues values(0.0);

        if constexpr (ParentType::isMomentumProblem())
        {
            values[Indices::velocityXIdx] = u_(x,y,time);
            values[Indices::velocityYIdx] = v_(x,y,time);
        }
        else
            values[Indices::pressureIdx] = p_(x,y,time);

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
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        const Scalar t = time_;
        Dune::FieldVector<GlobalPosition, DirichletValues::dimension> values;

        if constexpr (ParentType::isMomentumProblem())
        {
            values[Indices::velocityXIdx][0] = dxU_(x,y,t);
            values[Indices::velocityXIdx][1] = dyU_(x,y,t);
            values[Indices::velocityYIdx][0] = dxV_(x,y,t);
            values[Indices::velocityYIdx][1] = dyV_(x,y,t);
        }
        else
        {
            values[Indices::pressureIdx][0] = dxP_(x,y,t);
            values[Indices::pressureIdx][1] = dyP_(x,y,t);
        }

        return values;
    }

    /*!
     * \brief Applies the initial solution for all degrees of freedom
     * \note The analytical solution vanishes at t = 0 and the stationary problem starts from zero.
     */
    template<class SolutionVector>
    void applyInitialSolution(SolutionVector& sol) const
    {
        sol.resize(this->gridDiscretization().numDofs());
        sol = 0.0;
    }

    /*!
     * \brief Set the time at which sources and boundary conditions are evaluated
     * \note This is const because multi-stage assemblers set the stage time through a const problem.
     */
    void setTime(const Scalar time) const
    {
        time_ = time;
        updateConstraints_();
    }

    /*!
     * \brief Updates the time and the time-dependent Dirichlet constraints
     */
    void updateTime(const Scalar time)
    { setTime(time); }

    /*!
     * \brief The pressure acting on the momentum balance if not coupled to a mass problem
     */
    Scalar pressureAtPos(const GlobalPosition& globalPos) const
    { return p_(globalPos[0], globalPos[1], time_); }

    /*!
     * \brief The density if not coupled to a mass problem
     */
    Scalar densityAtPos(const GlobalPosition&) const
    { return rho_; }

    /*!
     * \brief The dynamic viscosity if not coupled to a mass problem
     */
    Scalar effectiveViscosityAtPos(const GlobalPosition&) const
    { return mu_; }

private:
    void init_()
    {
        time_ = 0.0;
        isStationary_ = getParam<bool>("Problem.IsStationary");
        rho_ = getParam<Scalar>("Component.LiquidDensity");
        const Scalar nu = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0);
        mu_ = rho_*nu;
        useNeumann_ = getParam<bool>("Problem.UseNeumann", false);

        updateConstraints_();
    }

    bool isMomentumFluxBoundary_(const GlobalPosition& globalPos) const
    {
        if (!useNeumann_)
            return false;

        static constexpr Scalar eps = 1e-6;
        const auto& bBoxMin = this->gridDiscretization().bBoxMin();
        const auto& bBoxMax = this->gridDiscretization().bBoxMax();
        return globalPos[1] > bBoxMax[1] - eps
            || globalPos[0] > bBoxMax[0] - eps
            || globalPos[1] < bBoxMin[1] + eps;
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
                    ConstraintValues values(analyticalSolution(scv.dofPosition(), time_)[Indices::pressureIdx]);
                    constraints_.push_back(DirichletConstraintData{std::move(info), std::move(values), scv.dofIndex()});
                    return;
                }
            }
        }
    }

    Scalar f_(Scalar t) const
    {
        using std::sin;
        if (isStationary_)
            return 1.0;
        else
            return sin(2.0 * t);
    }

    Scalar df_(Scalar t) const
    {
        using std::cos;
        if (isStationary_)
            return 0.0;
        else
            return 2.0 * cos(2.0 * t);
    }

    Scalar f1_(Scalar x) const
    { using std::cos; return -0.25 * cos(2.0 * x); }

    Scalar df1_(Scalar x) const
    { using std::sin; return 0.5 * sin(2.0 * x); }

    Scalar f2_(Scalar x) const
    { using std::cos; return -cos(x); }

    Scalar df2_(Scalar x) const
    { using std::sin; return sin(x); }

    Scalar ddf2_(Scalar x) const
    { using std::cos; return cos(x); }

    Scalar dddf2_(Scalar x) const
    { using std::sin; return -sin(x); }

    Scalar p_(Scalar x, Scalar y, Scalar t) const
    { return (f1_(x) + f1_(y)) * f_(t) * f_(t); }

    Scalar dxP_ (Scalar x, Scalar y, Scalar t) const
    { return df1_(x) * f_(t) * f_(t); }

    Scalar dyP_ (Scalar x, Scalar y, Scalar t) const
    { return df1_(y) * f_(t) * f_(t); }

    Scalar u_(Scalar x, Scalar y, Scalar t) const
    { return f2_(x)*df2_(y) * f_(t); }

    Scalar dtU_ (Scalar x, Scalar y, Scalar t) const
    { return f2_(x)*df2_(y) * df_(t); }

    Scalar dxU_ (Scalar x, Scalar y, Scalar t) const
    { return df2_(x)*df2_(y) * f_(t); }

    Scalar dyU_ (Scalar x, Scalar y, Scalar t) const
    { return f2_(x)*ddf2_(y) * f_(t); }

    Scalar dxxU_ (Scalar x, Scalar y, Scalar t) const
    { return ddf2_(x)*df2_(y) * f_(t); }

    Scalar dxyU_ (Scalar x, Scalar y, Scalar t) const
    { return df2_(x)*ddf2_(y) * f_(t); }

    Scalar dyyU_ (Scalar x, Scalar y, Scalar t) const
    { return f2_(x)*dddf2_(y) * f_(t); }

    Scalar v_(Scalar x, Scalar y, Scalar t) const
    { return -f2_(y)*df2_(x) * f_(t); }

    Scalar dtV_ (Scalar x, Scalar y, Scalar t) const
    { return -f2_(y)*df2_(x) * df_(t); }

    Scalar dxV_ (Scalar x, Scalar y, Scalar t) const
    { return -f2_(y)*ddf2_(x) * f_(t); }

    Scalar dyV_ (Scalar x, Scalar y, Scalar t) const
    { return -df2_(y)*df2_(x) * f_(t); }

    Scalar dyyV_ (Scalar x, Scalar y, Scalar t) const
    { return -ddf2_(y)*df2_(x) * f_(t); }

    Scalar dxyV_ (Scalar x, Scalar y, Scalar t) const
    { return -df2_(y)*ddf2_(x) * f_(t); }

    Scalar dxxV_ (Scalar x, Scalar y, Scalar t) const
    { return -f2_(y)*dddf2_(x) * f_(t); }

    Scalar dxUU_ (Scalar x, Scalar y, Scalar t) const
    { return 2.*u_(x,y,t)*dxU_(x,y,t); }

    Scalar dyVV_ (Scalar x, Scalar y, Scalar t) const
    { return 2.*v_(x,y,t)*dyV_(x,y,t); }

    Scalar dxUV_ (Scalar x, Scalar y, Scalar t) const
    { return v_(x,y,t)*dxU_(x,y,t) + u_(x,y,t)*dxV_(x,y,t); }

    Scalar dyUV_ (Scalar x, Scalar y, Scalar t) const
    { return v_(x,y,t)*dyU_(x,y,t) + u_(x,y,t)*dyV_(x,y,t); }

    Scalar rho_;
    Scalar mu_;
    mutable Scalar time_;
    bool isStationary_;
    bool useNeumann_;
    mutable std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux

#endif
