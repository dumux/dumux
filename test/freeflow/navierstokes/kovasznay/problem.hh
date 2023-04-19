// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Test for the staggered grid Navier-Stokes model with analytical solution (Kovasznay 1948, \cite Kovasznay1948)
 */

#ifndef DUMUX_KOVASZNAY_TEST_PROBLEM_HH
#define DUMUX_KOVASZNAY_TEST_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/staggered/problem.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the staggered grid (Kovasznay 1948, \cite Kovasznay1948)
 *
 * A two-dimensional Navier-Stokes flow with a periodicity in one direction
 * is considered. The set-up represents a wake behind a two-dimensional grid
 * and is chosen in a way such that an exact solution is available.
 */
template <class TypeTag>
class KovasznayTestProblem : public NavierStokesStaggeredProblem<TypeTag>
{
    using ParentType = NavierStokesStaggeredProblem<TypeTag>;

    using BoundaryTypes = Dumux::NavierStokesBoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr auto upwindSchemeOrder = getPropValue<TypeTag, Properties::UpwindSchemeOrder>();

public:
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    KovasznayTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        std::cout<< "upwindSchemeOrder is: " << GridGeometry::upwindStencilOrder() << "\n";
        rho_ = getParam<Scalar>("Component.LiquidDensity", 1.0);
        kinematicViscosity_ = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0);
        Scalar reynoldsNumber = 1.0 / kinematicViscosity_;
        lambda_ = 0.5 * reynoldsNumber
                        - std::sqrt(reynoldsNumber * reynoldsNumber * 0.25 + 4.0 * M_PI * M_PI);
    }

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
        // set fixed pressure in all cells at the left boundary
        auto isAtLeftBoundary = [&](const FVElementGeometry& fvGeometry)
        {
            if (fvGeometry.hasBoundaryScvf())
            {
                for (const auto& scvf : scvfs(fvGeometry))
                    if (scvf.boundary() && scvf.center()[0] < this->gridGeometry().bBoxMin()[0] + eps_)
                        return true;
            }
            return false;
        };
        return (isAtLeftBoundary(fvGeometry) && pvIdx == Indices::pressureIdx);
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
     * \brief Returns the analytical solution of the problem at a given position.
     *
     * \param globalPos The global position
     * \param time A parameter for consistent signatures. It is ignored here as this is a stationary test
     */
    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos, Scalar time = 0.0) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];

        PrimaryVariables values;
        values[Indices::pressureIdx] = rho_ * 0.5 * (1.0 - std::exp(2.0 * lambda_ * x));
        values[Indices::velocityXIdx] = 1.0 - std::exp(lambda_ * x) * std::cos(2.0 * M_PI * y);
        values[Indices::velocityYIdx] = 0.5 * lambda_ / M_PI * std::exp(lambda_ * x) * std::sin(2.0 * M_PI * y);

        return values;
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
        PrimaryVariables values;
        values[Indices::pressureIdx] = 0.0;
        values[Indices::velocityXIdx] = 0.0;
        values[Indices::velocityYIdx] = 0.0;

        return values;
    }

private:
    static constexpr Scalar eps_=1e-6;

    Scalar rho_;
    Scalar kinematicViscosity_;
    Scalar lambda_;
};
} // end namespace Dumux

#endif
