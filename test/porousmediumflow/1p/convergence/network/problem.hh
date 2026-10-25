// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Single-phase flow on a line or on a junction of straight branches with an
 *        analytical solution
 *
 * The branches leave the junction at the origin with unit directions \f$ \mathbf{d}_i \f$
 * and have unit length. With the arc length \f$ s = \mathbf{x}\cdot\mathbf{d}_i \f$ the
 * pressure on branch \f$ i \f$ is
 * \f[ p_i(s) = 1 + \beta_i s + \sin(\omega_i s), \qquad -p_i'' = \omega_i^2 \sin(\omega_i s), \f]
 * which is continuous at the junction, and \f$ \sum_i (\beta_i + \omega_i) = 0 \f$ makes the
 * fluxes balance there. On a single branch the solution can instead be the quartic
 * \f$ p = s(1-s)(1+s^2) \f$, whose source every quadrature rule of the schemes integrates exactly.
 * The free ends carry the analytical pressure, except that with `Problem.FluxAtEnd` the end
 * of the first branch carries the analytical outflow \f$ -\nabla p \cdot \mathbf{n} \f$ instead.
 */
#ifndef DUMUX_ONEP_NETWORK_CONVERGENCE_TEST_PROBLEM_HH
#define DUMUX_ONEP_NETWORK_CONVERGENCE_TEST_PROBLEM_HH

#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/problemwithspatialparams.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/discretization/dirichletconstraints.hh>

namespace Dumux {

template<class TypeTag>
class OnePNetworkConvergenceProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridDiscretization::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using BoundaryTypes = Dumux::Experimental::BoundaryTypes<ModelTraits::numEq()>;
    using ConstraintInfo = Dumux::DirichletConstraintInfo<ModelTraits::numEq()>;
    using ConstraintValues = Dune::FieldVector<Scalar, ModelTraits::numEq()>;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using DirichletConstraintData = Dumux::DirichletConstraintData<ConstraintInfo, ConstraintValues, GridIndexType>;
    static constexpr int dimWorld = GridView::dimensionworld;

    struct Branch
    {
        GlobalPosition direction;
        Scalar slope;
        Scalar omega;
    };

public:
    OnePNetworkConvergenceProblem(std::shared_ptr<const GridDiscretization> gridDiscretization)
    : ParentType(gridDiscretization)
    , quartic_(getParam<std::string>("Problem.Solution", "sine") == "quartic")
    , fluxAtEnd_(getParam<bool>("Problem.FluxAtEnd", false))
    {
        if constexpr (dimWorld == 1)
            branches_.push_back({GlobalPosition(1.0), 0.5, M_PI});
        else if constexpr (dimWorld == 3)
        {
            const auto unit = [](GlobalPosition d){ d /= d.two_norm(); return d; };
            const Scalar omega0 = M_PI, omega1 = 0.5*M_PI, omega2 = 1.5*M_PI;
            const Scalar slope0 = 0.5, slope1 = -1.0;
            const Scalar slope2 = -(slope0 + slope1 + omega0 + omega1 + omega2);
            branches_.push_back({unit({1.0, 0.0, 0.0}), slope0, omega0});
            branches_.push_back({unit({-0.5, 0.8, 0.1}), slope1, omega1});
            branches_.push_back({unit({-0.3, -0.6, 0.7}), slope2, omega2});
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Branches for dimWorld = " << dimWorld);

        if (quartic_ && branches_.size() != 1)
            DUNE_THROW(Dune::InvalidStateException, "The quartic solution is defined on a single branch only");

        appendDirichletConstraints_();
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        if (isFluxEnd_(globalPos))
            values.setFluxBoundary(0);
        return values;
    }

    NumEqVector boundaryFluxAtPos(const GlobalPosition& globalPos) const
    {
        if (!isFluxEnd_(globalPos))
            return NumEqVector(0.0);
        return NumEqVector(-(analyticalPressureGradient(globalPos)*branches_[0].direction));
    }

    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        const auto [branch, s] = branchAndArcLength_(globalPos);
        if (quartic_)
            return NumEqVector(2.0 - 6.0*s + 12.0*s*s);
        using std::sin;
        return NumEqVector(branch.omega*branch.omega*sin(branch.omega*s));
    }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    //! The analytical pressure
    Scalar analyticalPressure(const GlobalPosition& globalPos) const
    {
        const auto [branch, s] = branchAndArcLength_(globalPos);
        if (quartic_)
            return s*(1.0 - s)*(1.0 + s*s);
        using std::sin;
        return 1.0 + branch.slope*s + sin(branch.omega*s);
    }

    //! The gradient of the analytical pressure, along the branch
    GlobalPosition analyticalPressureGradient(const GlobalPosition& globalPos) const
    {
        const auto [branch, s] = branchAndArcLength_(globalPos);
        using std::cos;
        const Scalar dpds = quartic_ ? 1.0 - 2.0*s + 3.0*s*s - 4.0*s*s*s
                                     : branch.slope + branch.omega*cos(branch.omega*s);
        auto gradient = branch.direction;
        gradient *= dpds;
        return gradient;
    }

    const std::vector<DirichletConstraintData>& constraints() const
    { return constraints_; }

private:
    bool isFluxEnd_(const GlobalPosition& globalPos) const
    { return fluxAtEnd_ && (globalPos - branches_[0].direction).two_norm() < 1e-8; }

    //! the branch whose direction is closest to the position, and the arc length along it
    std::pair<Branch, Scalar> branchAndArcLength_(const GlobalPosition& globalPos) const
    {
        std::size_t best = 0;
        for (std::size_t i = 1; i < branches_.size(); ++i)
            if (globalPos*branches_[i].direction > globalPos*branches_[best].direction)
                best = i;
        return {branches_[best], globalPos*branches_[best].direction};
    }

    void appendDirichletConstraints_()
    {
        auto elemDisc = localView(this->gridDiscretization());
        for (const auto& element : elements(this->gridDiscretization().gridView()))
        {
            elemDisc.bind(element);
            for (const auto& boundaryFace : boundaryFaces(elemDisc))
                for (const auto& localDof : localDofs(elemDisc, boundaryFace))
                {
                    const auto globalPos = ipData(elemDisc, localDof).global();
                    if (isFluxEnd_(globalPos))
                        continue;
                    ConstraintInfo info;
                    info.setAll();
                    ConstraintValues values(analyticalPressure(globalPos));
                    constraints_.push_back(
                        DirichletConstraintData{std::move(info), std::move(values), localDof.dofIndex()}
                    );
                }
        }
    }

    bool quartic_;
    bool fluxAtEnd_;
    std::vector<Branch> branches_;
    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux

#endif
