// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief The convergence test problem stated in terms of the new interfaces
 */
#ifndef DUMUX_INCOMPRESSIBLE_ONEP_CONVERGENCE_TEST_PROBLEM_NEW_INTERFACE_HH
#define DUMUX_INCOMPRESSIBLE_ONEP_CONVERGENCE_TEST_PROBLEM_NEW_INTERFACE_HH

#include <cmath>
#include <vector>
#include <memory>
#include <utility>

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

/*!
 * \ingroup OnePTests
 * \brief The convergence test problem stated in terms of the new interfaces
 *
 * The analytic solution, the source term and the boundary values are those of
 * ConvergenceProblem. Dirichlet values reach the assembler as constraints, which
 * the problem collects once over the boundary faces of the grid.
 */
template<class TypeTag>
class ConvergenceProblemNewInterface : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;

    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridDiscretization::LocalView;
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

    static constexpr auto velocityXIdx = 0;
    static constexpr auto velocityYIdx = 1;
    static constexpr auto pressureIdx = 2;

public:
    /*!
     * \brief The constructor.
     * \param gridDiscretization The grid discretization
     */
    ConvergenceProblemNewInterface(std::shared_ptr<const GridDiscretization> gridDiscretization)
    : ParentType(gridDiscretization)
    , c_(getParam<Scalar>("Problem.C"))
    { appendDirichletConstraints_(); }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     * \param globalPos The position of the boundary segment
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        // the whole boundary is Dirichlet, which the assembler reads off the constraints;
        // a boundary type only marks the equations that receive a flux instead
        return BoundaryTypes{};
    }

    /*!
     * \brief Evaluates Dirichlet boundary conditions.
     * \param globalPos The position at which the value is evaluated
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(analyticalSolution(globalPos)[pressureIdx]); }

    /*!
     * \brief Evaluates the source term.
     * \param globalPos The position at which the source is evaluated
     */
    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        using std::exp;
        using std::sin;
        using std::cos;
        const Scalar cosOmegaX = cos(omega_*x);
        static const Scalar expTwo = exp(2);
        const Scalar expYPlusOne = exp(y+1);

        const Scalar result = ( -(c_*cosOmegaX + 1)*exp(y - 1)
                                + 1.5*c_*expYPlusOne*cosOmegaX
                                + omega_*omega_*(expYPlusOne - expTwo + 2))
                              * sin(omega_*x);

        return NumEqVector(result);
    }

    /*!
     * \brief Evaluates the initial value.
     * \param globalPos The position at which the value is evaluated
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    /*!
     * \brief Returns the analytical solution of the problem at a given position.
     * \param globalPos The global position
     */
    auto analyticalSolution(const GlobalPosition& globalPos) const
    {
        Dune::FieldVector<Scalar, 3> sol(0.0);
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        using std::exp; using std::sin; using std::cos;
        const Scalar sinOmegaX = sin(omega_*x);
        const Scalar cosOmegaX = cos(omega_*x);
        static const Scalar expTwo = exp(2);
        const Scalar expYPlusOne = exp(y+1);

        sol[pressureIdx] = (expYPlusOne + 2 - expTwo)*sinOmegaX + 10.0;
        sol[velocityXIdx] = c_/(2*omega_)*expYPlusOne*sinOmegaX*sinOmegaX
                            -omega_*(expYPlusOne + 2 - expTwo)*cosOmegaX;
        sol[velocityYIdx] = (0.5*c_*(expYPlusOne + 2 - expTwo)*cosOmegaX
                            -(c_*cosOmegaX + 1)*exp(y-1))*sinOmegaX;

        return sol;
    }

    //! The Dirichlet values the assembler constrains the solution to
    const std::vector<DirichletConstraintData>& constraints() const
    { return constraints_; }

private:
    //! Every local dof on the boundary carries the analytical pressure
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

                    ConstraintInfo info;
                    info.setAll();

                    ConstraintValues values(dirichletAtPos(globalPos));
                    constraints_.push_back(
                        DirichletConstraintData{std::move(info), std::move(values), localDof.dofIndex()}
                    );
                }
        }
    }

    static constexpr Scalar omega_ = M_PI;
    Scalar c_;
    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux

#endif
