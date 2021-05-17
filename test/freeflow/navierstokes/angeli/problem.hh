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
 * \ingroup NavierStokesTests
 * \brief Test for the instationary staggered grid Navier-Stokes model
 *        with analytical solution (Angeli et al. 2017, \cite Angeli2017).
 */

#ifndef DUMUX_ANGELI_TEST_PROBLEM_HH
#define DUMUX_ANGELI_TEST_PROBLEM_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/problem.hh>

#include "../l2error.hh"

namespace Dumux {
/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the staggered grid (Angeli et al. 2017, \cite Angeli2017).
 *
 * The unsteady, 2D, incompressible Navier-Stokes equations for a zero source and a Newtonian
 * flow is solved and compared to an analytical solution (sums/products of trigonometric functions).
 * The velocities and pressures decay exponentially. The Dirichlet boundary conditions are
 * time-dependent and consistent with the analytical solution.
 */
template <class TypeTag>
class AngeliTestProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using BoundaryTypes = Dumux::NavierStokesBoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

public:
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    AngeliTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        kinematicViscosity_ = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0);
        rho_ = getParam<Scalar>("Component.LiquidDensity", 1.0);
    }

   /*!
     * \brief Returns the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 298.0; }

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
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        // set Dirichlet values for the velocity everywhere
        values.setDirichlet(Indices::velocityXIdx);
        values.setDirichlet(Indices::velocityYIdx);

        return values;
    }

    /*!
     * \brief Returns whether a fixed Dirichlet value shall be used at a given cell.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param scv The sub control volume
     * \param pvIdx The primary variable index in the solution vector
     */
    bool isDirichletCell(const Element& element,
                         const typename GridGeometry::LocalView& fvGeometry,
                         const typename GridGeometry::SubControlVolume& scv,
                         int pvIdx) const
    {
        // set a fixed pressure in all cells at the boundary
        for (const auto& scvf : scvfs(fvGeometry))
            if (scvf.boundary())
                return true;

        return false;
    }

   /*!
     * \brief Returns Dirichlet boundary values at a given position.
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition & globalPos) const
    {
        // use the values of the analytical solution
        return this->instationaryAnalyticalSolution(globalPos, time_+timeStepSize_);
    }

    /*!
     * \brief Returns the analytical solution of the problem at a given position.
     *
     * \param globalPos The global position
     */
    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos) const
    {
        return this->instationaryAnalyticalSolution(globalPos, time_+timeStepSize_);
    }

    PrimaryVariables instationaryAnalyticalSolutionAtPos(const GlobalPosition& globalPos, const Scalar time) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        const Scalar t = time;

        PrimaryVariables values = {};

        for (unsigned int j = 0; j < 2; ++j)
        {
            for (unsigned int i = 0; i < values.size(); ++i)
            {
                values[i] += xFactorAnalyticalSolutionAtPos(x,t)[j][i]*yFactorAnalyticalSolutionAtPos(y,t)[j][i];
            }
        }

        return values;
    }

    /*!
     * \brief Returns the analytical solution of the problem at a given position.
     *
     * \param globalPos The global position
     */
    PrimaryVariables analyticalSolutionAtPos(const GlobalPosition& globalPos) const
    {
        return this->instationaryAnalyticalSolutionAtPos(globalPos, time_+timeStepSize_);
    }


    std::array<PrimaryVariables, 2> xFactorAnalyticalSolutionAtPos(Scalar x, Scalar t) const
    {
        PrimaryVariables values1;
        values1[Indices::pressureIdx] = - 0.25 * f_(t) * f_(t)  * 4.0 * M_PI * M_PI * std::cos(2.0 * M_PI * x)*rho_;
        values1[Indices::velocityXIdx] = f_(t) * f2_(x) ;
        values1[Indices::velocityYIdx] = - f_(t) * df2_(x);

        PrimaryVariables values2 = {};
        values2[Indices::pressureIdx] = 1.;

        std::array<PrimaryVariables, 2> retPair;
        retPair[0] = values1;
        retPair[1] = values2;

        return retPair;
    }

    std::array<PrimaryVariables, 2> yFactorAnalyticalSolutionAtPos(Scalar y, Scalar t) const
    {
        PrimaryVariables values1;
        values1[Indices::pressureIdx] = 1.;
        values1[Indices::velocityXIdx] = df3_(y);
        values1[Indices::velocityYIdx] = f3_(y);

        PrimaryVariables values2 = {};
        values2[Indices::pressureIdx] = - 0.25 * f_(t) * f_(t)  * M_PI * M_PI * std::cos(4.0 * M_PI * y)*rho_;

        std::array<PrimaryVariables, 2> retPair;
        retPair[0] = values1;
        retPair[1] = values2;

        return retPair;
    }

    std::array<PrimaryVariables, 2> xFactorAnalyticalSolutionAntiderivativeAtPos(Scalar x, Scalar t) const
    {
        PrimaryVariables values1;
        values1[Indices::pressureIdx] = - 0.25 * f_(t) * f_(t)  * 2.0 * M_PI * std::sin(2.0 * M_PI * x)*rho_;
        values1[Indices::velocityXIdx] = f_(t) * intf2_(x) ;
        values1[Indices::velocityYIdx] = - f_(t) * f2_(x);

        PrimaryVariables values2 = {};
        values2[Indices::pressureIdx] = x;

        std::array<PrimaryVariables, 2> retPair;
        retPair[0] = values1;
        retPair[1] = values2;

        return retPair;
    }

    std::array<PrimaryVariables, 2> yFactorAnalyticalSolutionAntiderivativeAtPos(Scalar y, Scalar t) const
    {
        PrimaryVariables values1;
        values1[Indices::pressureIdx] = y;
        values1[Indices::velocityXIdx] = f3_(y);
        values1[Indices::velocityYIdx] = intf3_(y);

        PrimaryVariables values2 = {};
        values2[Indices::pressureIdx] = - 1./16 * f_(t) * f_(t)  * M_PI * std::sin(4.0 * M_PI * y)*rho_;

        std::array<PrimaryVariables, 2> retPair;
        retPair[0] = values1;
        retPair[1] = values2;

        return retPair;
    }

    // \}

   /*!
     * \name Volume terms
     */
    // \{

   /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        return this->instationaryAnalyticalSolution(globalPos, 0);
    }

    /*!
     * \brief Updates the time
     */
    void updateTime(const Scalar time)
    {
        time_ = time;
    }

    /*!
     * \brief Updates the time step size
     */
    void updateTimeStepSize(const Scalar timeStepSize)
    {
        timeStepSize_ = timeStepSize;
    }

private:
    Scalar f_(Scalar t) const
    { return std::exp(- 5.0 * kinematicViscosity_ * M_PI * M_PI * t); }

    Scalar df_(Scalar t) const
    { return - 5.0 * kinematicViscosity_ * M_PI * M_PI *  f_(t); }

    Scalar intf2_ (Scalar x) const
    { return std::sin(M_PI * x)/M_PI; }

    Scalar f2_ (Scalar x) const
    { return std::cos(M_PI * x); }

    Scalar df2_ (Scalar x) const
    { return -M_PI * std::sin(M_PI * x); }

    Scalar intf3_ ( Scalar x) const
    { return std::sin(2.0 * M_PI * x)/(2.0 * M_PI); }

    Scalar f3_ ( Scalar x) const
    { return std::cos(2.0 * M_PI * x); }

    Scalar df3_ ( Scalar x) const
    { return - 2.0 * M_PI * std::sin(2.0 * M_PI * x); }

    Scalar kinematicViscosity_;
    Scalar rho_;
    Scalar time_ = 0;
    Scalar timeStepSize_ = 0;
};
} // end namespace Dumux

#endif
