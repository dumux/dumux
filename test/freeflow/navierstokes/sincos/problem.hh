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
 * \brief Test for the staggered grid Navier-Stokes model with analytical solution.
 */
#ifndef DUMUX_SINCOS_TEST_PROBLEM_HH
#define DUMUX_SINCOS_TEST_PROBLEM_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/problem.hh>

#include "../l2error.hh"

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the staggered grid.
 *
 * The 2D, incompressible Navier-Stokes equations for zero gravity and a Newtonian
 * flow is solved and compared to an analytical solution (sums/products of trigonometric functions).
 * For the instationary case, the velocities and pressures are periodical in time. The Dirichlet boundary conditions are
 * consistent with the analytical solution and in the instationary case time-dependent.
 */
template <class TypeTag>
class SincosTestProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using BoundaryTypes = Dumux::NavierStokesBoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using FVElementGeometry = typename GridGeometry::LocalView;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

public:
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    SincosTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry), time_(0.0), timeStepSize_(0.0)
    {
        isStationary_ = getParam<bool>("Problem.IsStationary");
        enableInertiaTerms_ = getParam<bool>("Problem.EnableInertiaTerms");
        rho_ = getParam<Scalar>("Component.LiquidDensity");
        //kinematic
        Scalar nu = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0);
        //dynamic
        mu_ = rho_*nu;
    }

   /*!
     * \brief Returns the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 298.0; }

   /*!
     * \brief Return the sources within the domain.
     *
     * \param globalPos The global position
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        NumEqVector source(0.0);
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        const Scalar t = time_ + timeStepSize_;

        using std::cos;
        using std::sin;

        source[Indices::momentumXBalanceIdx] = rho_*dtU_(x,y,t) - 2.*mu_*dxxU_(x,y,t) - mu_*dyyU_(x,y,t) - mu_*dxyV_(x,y,t) + dxP_(x,y,t);
        source[Indices::momentumYBalanceIdx] = rho_*dtV_(x,y,t) - 2.*mu_*dyyV_(x,y,t) - mu_*dxyU_(x,y,t) - mu_*dxxV_(x,y,t) + dyP_(x,y,t);

        if (enableInertiaTerms_)
        {
            source[Indices::momentumXBalanceIdx] += rho_*dxUU_(x,y,t) + rho_*dyUV_(x,y,t);
            source[Indices::momentumYBalanceIdx] += rho_*dxUV_(x,y,t) + rho_*dyVV_(x,y,t);
        }

        return source;
    }

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
                         const FVElementGeometry& fvGeometry,
                         const SubControlVolume& scv,
                         int pvIdx) const
    {
        // set fixed pressure in one cell
        return (scv.dofIndex() == 0) && pvIdx == Indices::pressureIdx;
    }

   /*!
     * \brief Returns Dirichlet boundary values at a given position.
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition & globalPos) const
    {
        // use the values of the analytical solution
        return analyticalSolution(globalPos, time_+timeStepSize_);
    }

    /*!
     * \brief Returns the analytical solution of the problem at a given time and position.
     *
     * \param globalPos The global position
     * \param time The current simulation time
     */
    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos, const Scalar time) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        const Scalar t = time;

        PrimaryVariables values;

        values[Indices::pressureIdx] = (f1_(x) + f1_(y)) * f_(t) * f_(t);
        values[Indices::velocityXIdx] = u_(x,y,t);
        values[Indices::velocityYIdx] = v_(x,y,t);

        return values;
    }

    /*!
     * \brief Returns the analytical solution of the problem at a given position.
     *
     * \param globalPos The global position
     */
    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos) const
    {
        return analyticalSolution(globalPos, time_+timeStepSize_);
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
        if (isStationary_)
        {
            PrimaryVariables values;
            values[Indices::pressureIdx] = 0.0;
            values[Indices::velocityXIdx] = 0.0;
            values[Indices::velocityYIdx] = 0.0;

            return values;
        }
        else
        {
            return analyticalSolution(globalPos, 0.0);
        }
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
    {
        if (isStationary_)
            return 1.0;
        else
            return std::sin(2.0 * t);
    }

    Scalar df_(Scalar t) const
    {
        if (isStationary_)
            return 0.0;
        else
            return 2.0 * std::cos(2.0 * t);
    }

    Scalar f1_(Scalar x) const
    { return -0.25 * std::cos(2.0 * x); }

    Scalar df1_(Scalar x) const
    { return 0.5 * std::sin(2.0 * x); }

    Scalar f2_(Scalar x) const
    { return - std::cos(x); }

    Scalar df2_(Scalar x) const
    { return std::sin(x); }

    Scalar ddf2_(Scalar x) const
    { return std::cos(x); }

    Scalar dddf2_(Scalar x) const
    { return -std::sin(x); }

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
    bool enableInertiaTerms_;
    Scalar time_;
    Scalar timeStepSize_;

    bool isStationary_;
};

} // end namespace Dumux

#endif
