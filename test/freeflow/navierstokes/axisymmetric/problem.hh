// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Axisymmetric Stokes flow with a manufactured solution with nonzero radial velocity
 *
 * Axial coordinate z = x, radial coordinate r = y (rotation about the x-axis). The velocity
 * derives from the Stokes stream function \f$ \psi = r^2 \cos(\pi r) \sin(\pi z) \f$,
 * \f$ u_r = -\partial_z\psi/r \f$, \f$ u_z = \partial_r\psi/r \f$, and is divergence-free, so the
 * symmetrized and the unsymmetrized viscous stress give the same source term. The pressure
 * \f$ p = \cos(\pi r)\cos(\pi z) \f$ is prescribed (no mass balance). The velocity is
 * prescribed on the whole boundary. The domain may be an annulus (lower radial bound > 0).
 *
 * The problem implements the classic interface (NavierStokesMomentumProblem) and the integral
 * interface (CVFENavierStokesMomentumProblem), where the Dirichlet conditions are constraints.
 */
#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_AXISYMMETRIC_PROBLEM_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_AXISYMMETRIC_PROBLEM_HH

#include <cmath>
#include <memory>
#include <vector>

#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/common/typetraits/griddiscretization.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/discretization/dirichletconstraints.hh>

namespace Dumux {

template<class TypeTag, class BaseProblem>
class AxisymmetricStokesTestProblem : public BaseProblem
{
    using ParentType = BaseProblem;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using GlobalPosition = typename GridGeometry::GlobalCoordinate;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    static constexpr int numEq = ModelTraits::numEq();
    using ConstraintInfo = DirichletConstraintInfo<numEq>;
    using ConstraintValues = Dune::FieldVector<Scalar, numEq>;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using DirichletConstraint = DirichletConstraintData<ConstraintInfo, ConstraintValues, GridIndexType>;

public:
    using Indices = typename ModelTraits::Indices;
    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using DirichletValues = typename ParentType::DirichletValues;
    using Sources = typename ParentType::Sources;

    static constexpr bool integralInterface = requires (BoundaryTypes values) { values.setAllFluxBoundary(); };

    AxisymmetricStokesTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        density_ = getParam<Scalar>("Component.LiquidDensity");
        viscosity_ = getParam<Scalar>("Component.LiquidKinematicViscosity")*density_;

        if constexpr (integralInterface)
            updateConstraints_();
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition&) const
    {
        BoundaryTypes values;
        if constexpr (!integralInterface)
            values.setAllDirichlet();
        return values;
    }

    DirichletValues dirichletAtPos(const GlobalPosition& globalPos) const
    { return analyticalSolution(globalPos); }

    const std::vector<DirichletConstraint>& constraints() const
    requires integralInterface
    { return constraints_; }

    Sources sourceAtPos(const GlobalPosition& globalPos) const
    {
        using std::sin; using std::cos;
        const Scalar z = globalPos[0];
        const Scalar r = globalPos[1];
        const Scalar mu = viscosity_;
        const Scalar pi = M_PI;
        // sin(pi r)/r with its limit pi on the axis
        const Scalar sinOverR = r > 1e-12 ? sin(pi*r)/r : pi;

        Sources source(0.0);
        source[Indices::momentumXBalanceIdx] = -pi*sin(pi*z)*(2.0*pi*pi*mu*r*sin(pi*r) - 7.0*pi*mu*cos(pi*r)
                                                              - 3.0*mu*sinOverR + cos(pi*r));
        source[Indices::momentumYBalanceIdx] = -pi*cos(pi*z)*(2.0*pi*pi*mu*r*cos(pi*r) + 3.0*pi*mu*sin(pi*r) + sin(pi*r));
        return source;
    }

    template<class SolutionVector>
    void applyInitialSolution(SolutionVector& sol) const
    {
        sol.resize(Dumux::gridDiscretization(*this).numDofs());
        sol = 0.0;
    }

    DirichletValues analyticalSolution(const GlobalPosition& globalPos) const
    {
        using std::sin; using std::cos;
        const Scalar z = globalPos[0];
        const Scalar r = globalPos[1];
        const Scalar pi = M_PI;
        DirichletValues values(0.0);
        values[Indices::velocityXIdx] = (2.0*cos(pi*r) - pi*r*sin(pi*r))*sin(pi*z);
        values[Indices::velocityYIdx] = -pi*r*cos(pi*r)*cos(pi*z);
        return values;
    }

    Dune::FieldVector<GlobalPosition, 2> gradAnalyticalSolution(const GlobalPosition& globalPos) const
    {
        using std::sin; using std::cos;
        const Scalar z = globalPos[0];
        const Scalar r = globalPos[1];
        const Scalar pi = M_PI;
        Dune::FieldVector<GlobalPosition, 2> values;
        values[Indices::velocityXIdx][0] = pi*(2.0*cos(pi*r) - pi*r*sin(pi*r))*cos(pi*z);
        values[Indices::velocityXIdx][1] = -pi*(pi*r*cos(pi*r) + 3.0*sin(pi*r))*sin(pi*z);
        values[Indices::velocityYIdx][0] = pi*pi*r*cos(pi*r)*sin(pi*z);
        values[Indices::velocityYIdx][1] = pi*(pi*r*sin(pi*r) - cos(pi*r))*cos(pi*z);
        return values;
    }

    Scalar pressureAtPos(const GlobalPosition& globalPos) const
    { return std::cos(M_PI*globalPos[1])*std::cos(M_PI*globalPos[0]); }

    Scalar densityAtPos(const GlobalPosition&) const
    { return density_; }

    Scalar effectiveViscosityAtPos(const GlobalPosition&) const
    { return viscosity_; }

private:
    void updateConstraints_()
    {
        constraints_.clear();
        const auto& gridGeometry = Dumux::gridDiscretization(*this);
        auto fvGeometry = localView(gridGeometry);
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            fvGeometry.bind(element);
            for (const auto& boundaryFace : boundaryFaces(fvGeometry))
            {
                for (const auto& localDof : localDofs(fvGeometry, boundaryFace))
                {
                    ConstraintInfo info;
                    info.setAll();
                    ConstraintValues values(analyticalSolution(ipData(fvGeometry, localDof).global()));
                    constraints_.push_back(DirichletConstraint{std::move(info), std::move(values), localDof.dofIndex()});
                }
            }
        }
    }

    Scalar density_, viscosity_;
    std::vector<DirichletConstraint> constraints_;
};

} // end namespace Dumux

#endif
