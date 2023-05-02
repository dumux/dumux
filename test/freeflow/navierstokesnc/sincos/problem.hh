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
    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    static constexpr auto compIdx = 1;
public:
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    SincosTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        enableInertiaTerms_ = getParam<bool>("Problem.EnableInertiaTerms");
        rho_ = getParam<Scalar>("Component.LiquidDensity");
        //kinematic
        Scalar nu = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0);
        //dynamic
        mu_ = rho_*nu;
        diffusionCoeff_ = getParam<Scalar>("Component.DiffusionCoefficient");
        useNeumann_ = getParam<bool>("Problem.UseNeumann", false);
        constrainCell_ = getParam<bool>("Problem.ConstrainCell",true);
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

        source[Indices::momentumXBalanceIdx] = -2.*mu_*dxxU_(x,y) - mu_*dyyU_(x,y) - mu_*dxyV_(x,y) + dxP_(x,y);
        source[Indices::momentumYBalanceIdx] = -2.*mu_*dyyV_(x,y) - mu_*dxyU_(x,y) - mu_*dxxV_(x,y) + dyP_(x,y);

        if (enableInertiaTerms_)
        {
            source[Indices::momentumXBalanceIdx] += rho_*dxUU_(x,y) + rho_*dyUV_(x,y);
            source[Indices::momentumYBalanceIdx] += rho_*dxUV_(x,y) + rho_*dyVV_(x,y);
        }
            source[Indices::conti0EqIdx] = 0.0;
            source[Indices::conti0EqIdx + compIdx] = rho_*(dxUC_(x,y) + dyVC_(x,y)) - diffusionCoeff_*rho_*(dxxC_(x,y)+dyyC_(x,y));
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
        if (scv.dofIndex() == 0)
        {
            if (pvIdx == Indices::pressureIdx)
                return true;
            if (pvIdx == Indices::conti0EqIdx + compIdx)
                return true;
        }
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
        return analyticalSolution(globalPos);
    }

    /*!
     * \brief Returns the analytical solution of the problem at a given time and position.
     *
     * \param globalPos The global position
     * \param time The current simulation time
     */
    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        PrimaryVariables values;

        values[Indices::pressureIdx] = p_(x,y);
        values[Indices::conti0EqIdx + compIdx] = c_(x,y);
        values[Indices::velocityXIdx] = u_(x,y);
        values[Indices::velocityYIdx] = v_(x,y);

        return values;
    }

    // \}

   /*!
     * \name Volume terms
     */
    // \{


    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return PrimaryVariables(0.0); }


private:
    Scalar f_() const
    { return 1.0; }

    Scalar df_() const
    { return 0.0; }

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

    Scalar f3_(Scalar x) const
    { using std::cos; return 0.25*(cos(x)+1.0); }

    Scalar df3_(Scalar x) const
    { using std::sin; return -0.25*sin(x); }

    Scalar ddf3_(Scalar x) const
    { using std::cos; return -0.25*cos(x); }

    Scalar p_(Scalar x, Scalar y) const
    { return (f1_(x) + f1_(y)); }

    Scalar dxP_ (Scalar x, Scalar y) const
    { return df1_(x); }

    Scalar dyP_ (Scalar x, Scalar y) const
    { return df1_(y); }

    Scalar u_(Scalar x, Scalar y) const
    { return f2_(x)*df2_(y); }

    Scalar dxU_ (Scalar x, Scalar y) const
    { return df2_(x)*df2_(y); }

    Scalar dyU_ (Scalar x, Scalar y) const
    { return f2_(x)*ddf2_(y); }

    Scalar dxxU_ (Scalar x, Scalar y) const
    { return ddf2_(x)*df2_(y); }

    Scalar dxyU_ (Scalar x, Scalar y) const
    { return df2_(x)*ddf2_(y); }

    Scalar dyyU_ (Scalar x, Scalar y) const
    { return f2_(x)*dddf2_(y); }

    Scalar v_(Scalar x, Scalar y) const
    { return -f2_(y)*df2_(x); }

    Scalar dxV_ (Scalar x, Scalar y) const
    { return -f2_(y)*ddf2_(x); }

    Scalar dyV_ (Scalar x, Scalar y) const
    { return -df2_(y)*df2_(x); }

    Scalar dyyV_ (Scalar x, Scalar y) const
    { return -ddf2_(y)*df2_(x); }

    Scalar dxyV_ (Scalar x, Scalar y) const
    { return -df2_(y)*ddf2_(x); }

    Scalar dxxV_ (Scalar x, Scalar y) const
    { return -f2_(y)*dddf2_(x); }

    Scalar dxUU_ (Scalar x, Scalar y) const
    { return 2.*u_(x,y)*dxU_(x,y); }

    Scalar dyVV_ (Scalar x, Scalar y) const
    { return 2.*v_(x,y)*dyV_(x,y); }

    Scalar dxUV_ (Scalar x, Scalar y) const
    { return (v_(x,y) * dxU_(x,y)) + (u_(x,y) * dxV_(x,y)); }

    Scalar dyUV_ (Scalar x, Scalar y) const
    { return (v_(x,y) * dyU_(x,y)) + (u_(x,y) * dyV_(x,y)); }

    Scalar c_(Scalar x, Scalar y) const
    { return (f3_(x) + f3_(y)); }

    Scalar dxC_(Scalar x, Scalar) const
    { return df3_(x); }

    Scalar dyC_(Scalar, Scalar y) const
    { return df3_(y); }

    Scalar dxxC_(Scalar x, Scalar) const
    { return ddf3_(x); }

    Scalar dyyC_(Scalar, Scalar y) const
    { return ddf3_(y); }

    Scalar dxUC_(Scalar x, Scalar y) const
    { return (dxU_(x,y) * c_(x,y)) + (dxC_(x,y) * u_(x,y)); }

    Scalar dyVC_(Scalar x, Scalar y) const
    { return (dyV_(x,y) * c_(x,y)) + (dyC_(x,y) * v_(x,y)); }
    Scalar rho_;
    Scalar mu_;
    bool enableInertiaTerms_;
    Scalar time_;
    Scalar timeStepSize_;
    Scalar diffusionCoeff_;
    bool useNeumann_;
    bool constrainCell_;
};

} // end namespace Dumux

#endif
