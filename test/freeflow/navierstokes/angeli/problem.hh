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

#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/problem.hh>

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
    using GridView = typename GridGeometry::GridView;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    AngeliTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        kinematicViscosity_ = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0);
        rho_ = getParam<Scalar>("Component.LiquidDensity", 1.0);
        useVelocityAveragingForDirichlet_ = getParam<bool>("Problem.UseVelocityAveragingForDirichlet", false);
        useVelocityAveragingForInitial_ = getParam<bool>("Problem.UseVelocityAveragingForInitial", false);
    }

    /*!
     * \brief Returns the temperature within the domain in [K].
     * This problem assumes a temperature of 20 degrees Celsius (unused)
     */
    Scalar temperature() const
    { return 293.15; }

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
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume face (velocities)
     *
     * \param element The finite element
     * \param scvf the sub control volume face
     *
     * Concerning the usage of averagedVelocity_, see the explanation of the initial function.
     */
    PrimaryVariables dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    {
        if (useVelocityAveragingForDirichlet_)
            return averagedVelocity_(scvf, time_);
        else
            return analyticalSolution(scvf.center(), time_);
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume face (pressure)
     *
     * \param element The finite element
     * \param scv the sub control volume
     */
    PrimaryVariables dirichlet(const Element& element, const SubControlVolume& scv) const
    {
        PrimaryVariables priVars(0.0);
        priVars[Indices::pressureIdx] = analyticalSolution(scv.center(), time_)[Indices::pressureIdx];
        return priVars;
    }

    /*!
     * \brief Returns the analytical solution of the problem at a given time and position.
     * \param globalPos The global position
     * \param time The current simulation time
     */
    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos, Scalar time) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        PrimaryVariables values;

        values[Indices::pressureIdx] = - 0.25 * std::exp(-10.0 * kinematicViscosity_ * M_PI * M_PI * time) * M_PI * M_PI * (4.0 * std::cos(2.0 * M_PI * x) + std::cos(4.0 * M_PI * y))*rho_;
        values[Indices::velocityXIdx] = - 2.0 * M_PI * std::exp(- 5.0 * kinematicViscosity_ * M_PI * M_PI * time) * std::cos(M_PI * x) * std::sin(2.0 * M_PI * y);
        values[Indices::velocityYIdx] = M_PI * std::exp(- 5.0 * kinematicViscosity_ * M_PI * M_PI * time) * std::sin(M_PI * x) * std::cos(2.0 * M_PI * y);

        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the initial value for a control volume (pressure)
     */
    PrimaryVariables initial(const SubControlVolume& scv) const
    {
        PrimaryVariables priVars(0.0);
        priVars[Indices::pressureIdx] = analyticalSolution(scv.center(), time_)[Indices::pressureIdx];
        return priVars;
    }

    /*!
     * \brief Evaluates the initial value for a sub control volume face (velocities)
     *
     * Simply assigning the value of the analytical solution at the face center
     * gives a discrete solution that is not divergence-free. For small initial
     * time steps, this has a negative impact on the pressure solution
     * after the first time step. The flag UseVelocityAveragingForInitial triggers the
     * function averagedVelocity_ which uses a higher order quadrature formula to
     * bring the discrete solution sufficiently close to being divergence-free.
     */
    PrimaryVariables initial(const SubControlVolumeFace& scvf) const
    {
        if (useVelocityAveragingForInitial_)
            return averagedVelocity_(scvf, time_);
        else
            return analyticalSolution(scvf.center(), time_);
    }

    // \}

    /*!
     * \brief Updates the time information
     */
    void setTime(Scalar t)
    {
        time_ = t;
    }

private:
    PrimaryVariables averagedVelocity_(const SubControlVolumeFace& scvf, Scalar t) const
    {
        PrimaryVariables priVars(0.0);
        const auto geo = scvf.geometry();
        const auto& quad = Dune::QuadratureRules<Scalar, GridView::dimension-1>::rule(geo.type(), 3);
        for (auto&& qp : quad)
        {
            const auto w = qp.weight()*geo.integrationElement(qp.position());
            const auto globalPos = geo.global(qp.position());
            const auto sol = analyticalSolution(globalPos, t);
            priVars[Indices::velocityXIdx] += sol[Indices::velocityXIdx]*w;
            priVars[Indices::velocityYIdx] += sol[Indices::velocityYIdx]*w;
        }
        priVars[Indices::velocityXIdx] /= scvf.area();
        priVars[Indices::velocityYIdx] /= scvf.area();
        return priVars;
    }

    Scalar kinematicViscosity_, rho_;
    Scalar time_ = 0;
    bool useVelocityAveragingForDirichlet_;
    bool useVelocityAveragingForInitial_;
};
} // end namespace Dumux

#endif
