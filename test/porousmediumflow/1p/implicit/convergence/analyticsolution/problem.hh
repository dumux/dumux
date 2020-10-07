// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup OnePTests
 * \brief The problem setup for the convergence test with analytic solution
 */
#ifndef DUMUX_CONVERGENCE_TEST_ONEP_PROBLEM_HH
#define DUMUX_CONVERGENCE_TEST_ONEP_PROBLEM_HH

#include <cmath>
#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/porousmediumflow/problem.hh>
#include "testcase.hh"

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief The problem setup for the convergence test with analytic solution
 */
template <class TypeTag>
class ConvergenceProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr int dimWorld = GridView::dimensionworld;

public:
    /*!
     * \brief The constructor.
     * \param gridGeometry The finite-volume grid geometry
     */
    ConvergenceProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<typename ParentType::SpatialParams> spatialParams, const TestCase testCase)
    : ParentType(gridGeometry, spatialParams)
    , testCase_(testCase)
    , c_(getParam<Scalar>("Problem.C"))
    , freqFactor_(getParam<Scalar>("Problem.FreqFactor", 0.5))
    {}

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the temperature within the domain in [K].
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10°C
    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        switch (testCase_)
        {
            case TestCase::Schneider:
                values.setAllDirichlet();
                break;
            case TestCase::Sinus:
                values.setAllDirichlet();
                break;
            default:
                DUNE_THROW(Dune::InvalidStateException, "Invalid test case");
        }
        return values;
    }

    /*!
     * \brief Evaluates Dirichlet boundary conditions.
     * \param globalPos The center of the finite volume which ought to be set.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        const auto p = analyticalSolution(globalPos);
        return PrimaryVariables(p);
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the source term for all phases within a given
     *        sub-control volume.
     *
     * For this method, the \a priVars parameter stores the rate mass
     * of a component is generated or annihilated per volume
     * unit. Positive values mean that mass is created, negative ones
     * mean that it vanishes.
     *
     * The units must be according to either using mole or mass fractions (mole/(m^3*s) or kg/(m^3*s)).
     */
    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        switch (testCase_)
        {
            case TestCase::Schneider:
                return sourceSchneiderEtAl_(globalPos);
            case TestCase::Sinus:
                return sourceSinus_(globalPos);
            default:
                DUNE_THROW(Dune::InvalidStateException, "Invalid test case");
        }
    }

    // \}

    /*!
     * \brief Evaluates the initial value for a control volume.
     * \param globalPos The position for which the initial condition should be evaluated
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    /*!
     * \brief Returns the analytical solution of the problem at a given position.
     * \param globalPos The global position
     */
    auto analyticalSolution(const GlobalPosition& globalPos) const
    {
        switch (testCase_)
        {
            case TestCase::Schneider:
                return analyticalSolutionSchneiderEtAl_(globalPos);
            case TestCase::Sinus:
                return analyticalSolutionSinus_(globalPos);
            default:
                DUNE_THROW(Dune::InvalidStateException, "Invalid test case");
        }
    }


private:

    // see Schneider et al., 2019: "Coupling staggered-grid and MPFA finite volume methods for
    // free flow/porous-medium flow problems"
    auto analyticalSolutionSchneiderEtAl_(const GlobalPosition& x) const
    {
        using std::exp; using std::sin; using std::cos;

        return (exp(x[1]+1) + 2 - exp(2))*sin(omega_*x[0]) + 10.0;
    }

    auto analyticalSolutionSinus_(const GlobalPosition& x) const
    {
        Scalar sol(0.0);

        using std::sin; using std::cos;
        if constexpr (dimWorld == 2)
            sol = sin(freqFactor_*omega_*x[0]) * sin(freqFactor_*omega_*x[1]);
        else if constexpr (dimWorld == 3)
            sol = sin(freqFactor_*omega_*x[0]) * sin(freqFactor_*omega_*x[1]) * sin(freqFactor_*omega_*x[2]);

        return sol;
    }

    auto sourceSchneiderEtAl_(const GlobalPosition& x) const
    {
        using std::exp; using std::sin; using std::cos;
        const Scalar cosOmegaX = cos(omega_*x[0]);
        static const Scalar expTwo = exp(2);
        const Scalar expYPlusOne = exp(x[1]+1);

        const Scalar source = ( -(c_*cosOmegaX + 1)*exp(x[1] - 1)
                                + 1.5*c_*expYPlusOne*cosOmegaX
                                + omega_*omega_*(expYPlusOne - expTwo + 2))
                              * sin(omega_*x[0]);

        return NumEqVector(source);
    }

    auto sourceSinus_(const GlobalPosition& x) const
    {
        using std::sin; using std::cos;
        const auto K = this->spatialParams().permeabilityAtPos(x);

        Scalar source(0.0);
        if constexpr (dimWorld == 2)
            source =  freqFactor_*omega_ * freqFactor_*omega_
                    * sin(freqFactor_*omega_*x[0]) * sin(freqFactor_*omega_*x[1])
                    * (K[0][0] + K[1][1]);
        else if constexpr (dimWorld == 3)
            source =  freqFactor_*omega_ * freqFactor_*omega_
                    * sin(freqFactor_*omega_*x[0]) * sin(freqFactor_*omega_*x[1]) * sin(freqFactor_*omega_*x[2])
                    * (K[0][0] + K[1][1] + K[2][2]);

        return NumEqVector(source);
    }

    static constexpr Scalar eps_ = 1e-7;
    static constexpr Scalar omega_ = M_PI;
    TestCase testCase_;
    Scalar c_;
    Scalar freqFactor_;
};

} // end namespace Dumux

#endif
