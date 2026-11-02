// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Test problem for the advection of a density interface with uniform velocity.
 */
#ifndef DUMUX_DENSITY_ADVECTION_TEST_PROBLEM_HH
#define DUMUX_DENSITY_ADVECTION_TEST_PROBLEM_HH

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include <dune/common/fvector.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/discretization/dirichletconstraints.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief Momentum balance with a prescribed density field and a constant pressure
 *
 * The density \f$ \rho(x - u_0 t) \f$ is a smooth step between a left and a right value which is
 * translated with the velocity \f$ u_0 \f$. Without viscosity and with a constant pressure, the momentum
 * balance \f$ \partial_t (\rho u) + \nabla\cdot(\rho u \otimes u) = 0 \f$ is solved by the uniform
 * velocity \f$ u = u_0 \f$. The velocity is prescribed on the entire boundary.
 */
template <class TypeTag, class BaseProblem>
class DensityAdvectionTestProblem : public BaseProblem
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
    using Element = typename ElementDiscretization::Element;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using Indices = typename ModelTraits::Indices;

    explicit DensityAdvectionTestProblem(std::shared_ptr<const GridDiscretization> gridDiscretization)
    : ParentType(gridDiscretization)
    {
        velocity_ = getParam<Scalar>("Problem.Velocity");
        interfacePosition_ = getParam<Scalar>("Problem.InterfacePosition");
        interfaceWidth_ = getParam<Scalar>("Problem.InterfaceWidth");
        updateConstraints_();
    }

    Sources sourceAtPos(const GlobalPosition&) const
    { return Sources(0.0); }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition&) const
    { return BoundaryTypes(); }

    DirichletValues dirichletAtPos(const GlobalPosition&) const
    { return exactVelocity(); }

    DirichletValues exactVelocity() const
    {
        DirichletValues values(0.0);
        values[Indices::velocityXIdx] = velocity_;
        return values;
    }

    const auto& constraints() const
    { return constraints_; }

    template<class ElementVariables, class FaceIpData>
    BoundaryFluxes boundaryFlux(const ElementDiscretization&, const ElementVariables&, const FaceIpData&) const
    { return BoundaryFluxes(0.0); }

    template<class ElementVariables, class ElementFluxVariablesCache, class FaceIpData>
    BoundaryFluxes boundaryFlux(const ElementDiscretization&, const ElementVariables&,
                                const ElementFluxVariablesCache&, const FaceIpData&) const
    { return BoundaryFluxes(0.0); }

    template<class SolutionVector>
    void applyInitialSolution(SolutionVector& sol) const
    {
        sol.resize(this->gridDiscretization().numDofs());
        sol = 0.0;
        auto elemDisc = localView(this->gridDiscretization());
        for (const auto& element : elements(this->gridDiscretization().gridView()))
        {
            elemDisc.bind(element);
            for (const auto& scv : scvs(elemDisc))
                sol[scv.dofIndex()][Indices::velocityXIdx] = velocity_;
            // the first bubble dof holds the nodal value at the element center, the remaining ones are corrections
            for (const auto& localDof : nonCVLocalDofs(elemDisc))
            {
                sol[localDof.dofIndex()][Indices::velocityXIdx] = velocity_;
                break;
            }
        }
    }

    //! The prescribed density at a position, translated with the velocity
    Scalar exactDensityAt(const GlobalPosition& globalPos, const Scalar time) const
    {
        using std::tanh;
        const Scalar s = 0.5*(1.0 + tanh((globalPos[0] - interfacePosition_ - velocity_*time)/interfaceWidth_));
        return densityLeft_ + (densityRight_ - densityLeft_)*s;
    }

    /*!
     * \brief Density that is transported by a one-dimensional finite-volume scheme on the nodes of the grid
     *
     * The density is piecewise linear in x and given by implicit Euler steps of the advection equation
     * with the face value taken from the upstream node (upwind) or as the mean of the two nodes (central).
     */
    void computeTransportedDensity(const bool upwind, const int numCells, const Scalar length,
                                   const Scalar dt, const int numSteps)
    {
        transported_ = true;
        nodeSpacing_ = length/numCells;
        dt_ = dt;
        const int n = numCells + 1;
        const Scalar c = velocity_*dt/nodeSpacing_;
        nodalDensity_.assign(numSteps + 1, std::vector<Scalar>(n));
        for (int i = 0; i < n; ++i)
            nodalDensity_[0][i] = exactDensityAt({i*nodeSpacing_, 0.0}, 0.0);

        for (int step = 1; step <= numSteps; ++step)
        {
            // tridiagonal system a_i rho_{i-1} + b_i rho_i + d_i rho_{i+1} = rhs_i with Dirichlet end values
            std::vector<Scalar> a(n, 0.0), b(n, 1.0), d(n, 0.0), rhs = nodalDensity_[step-1];
            for (int i = 1; i < n-1; ++i)
            {
                if (upwind)
                { a[i] = -c; b[i] = 1.0 + c; }
                else
                { a[i] = -0.5*c; d[i] = 0.5*c; }
            }
            for (int i = 1; i < n; ++i)
            {
                const Scalar m = a[i]/b[i-1];
                b[i] -= m*d[i-1];
                rhs[i] -= m*rhs[i-1];
            }
            auto& rho = nodalDensity_[step];
            rho[n-1] = rhs[n-1]/b[n-1];
            for (int i = n-2; i >= 0; --i)
                rho[i] = (rhs[i] - d[i]*rho[i+1])/b[i];
        }
    }

    Scalar densityAt(const GlobalPosition& globalPos, const Scalar time) const
    {
        if (!transported_)
            return exactDensityAt(globalPos, time);

        using std::floor;
        using std::round;
        const auto& nodal = nodalDensity_[static_cast<std::size_t>(round(time/dt_))];
        const Scalar xi = std::clamp(globalPos[0]/nodeSpacing_, Scalar(0.0), Scalar(nodal.size() - 1));
        const std::size_t i = std::min<std::size_t>(static_cast<std::size_t>(floor(xi)), nodal.size() - 2);
        return nodal[i] + (xi - i)*(nodal[i+1] - nodal[i]);
    }

    //! The density at the current or the previous time level
    template<class IpData>
    Scalar density(const Element&, const ElementDiscretization&, const IpData& ipData,
                   const bool isPreviousTimeStep = false) const
    { return densityAt(ipData.global(), isPreviousTimeStep ? timePrevious_ : time_); }

    //! Set the density on the left and on the right of the interface
    void setDensities(const Scalar densityLeft, const Scalar densityRight)
    {
        densityLeft_ = densityLeft;
        densityRight_ = densityRight;
    }

    //! Set the time at the end of the time step and the one at its beginning
    void setTimeStep(const Scalar timePrevious, const Scalar time)
    {
        timePrevious_ = timePrevious;
        time_ = time;
    }

    Scalar pressureAtPos(const GlobalPosition&) const
    { return 0.0; }

    Scalar effectiveViscosityAtPos(const GlobalPosition&) const
    { return 0.0; }

private:
    void updateConstraints_()
    {
        constraints_.clear();
        auto elemDisc = localView(this->gridDiscretization());
        for (const auto& element : elements(this->gridDiscretization().gridView()))
        {
            elemDisc.bind(element);
            for (const auto& boundaryFace : boundaryFaces(elemDisc))
            {
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

    Scalar velocity_, interfacePosition_, interfaceWidth_;
    Scalar densityLeft_ = 1.0, densityRight_ = 1.0;
    Scalar time_ = 0.0, timePrevious_ = 0.0;
    bool transported_ = false;
    Scalar nodeSpacing_ = 1.0, dt_ = 1.0;
    std::vector<std::vector<Scalar>> nodalDensity_;
    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux

#endif
