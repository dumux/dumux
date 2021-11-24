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
 * \brief The properties & problem setup for the convergence test
 */
#ifndef DUMUX_INCOMPRESSIBLE_ONEP_CONVERGENCETEST_PROBLEM_HH
#define DUMUX_INCOMPRESSIBLE_ONEP_CONVERGENCETEST_PROBLEM_HH

#include <cmath>

#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief problem setup for the convergence test
 */
template <class TypeTag>
class OnePTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    /*!
     * \brief The constructor.
     * \param gridGeometry The finite-volume grid geometry
     */
    OnePTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     * \param globalPos The center of the finite volume for which it is to be set.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return exact(globalPos);
    }

    /*!
     * \brief Evaluates the source term within a sub-control volume.
     * \param element The finite element
     * \param fvGeometry The element finite-volume geometry
     * \param elemVolVars The element volume variables
     * \param scv The sub-control volume for which the source term is evaluated
     */
    template <class ElementVolumeVariables>
    NumEqVector source(const Element &element,
                        const FVElementGeometry &fvGeometry,
                        const ElementVolumeVariables &elemVolVars,
                        const SubControlVolume &scv) const
    {
        static const auto order = getParam<Scalar>("Problem.SourceIntegrationOrder");
        static const auto periodLength = getParam<Scalar>("Problem.ExactSolPeriodLength");
        const auto &k = this->spatialParams().permeabilityAtPos(scv.center());

        using std::cos;
        using std::sin;

        const auto eg = element.geometry();
        const auto rule = Dune::QuadratureRules<Scalar, GridView::dimension>::rule(eg.type(), order);

        Scalar source = 0.0;
        for (auto qp : rule)
        {
            const auto p = eg.global(qp.position());
            const auto x = p[0];
            const auto y = p[1];

            const auto preFactor = -1.0 * periodLength * periodLength * M_PI * M_PI;
            const auto preFactorArg = periodLength * M_PI;
            const auto sineTerm = sin(preFactorArg * x);
            const auto cosTerm = cos(preFactorArg * y);
            const auto secondDeriv = preFactor * sineTerm * cosTerm;

            // derivative in x and y are identical
            source -= 2.0 * k * secondDeriv * qp.weight() * eg.integrationElement(qp.position());
        }

        source /= eg.volume();
        return NumEqVector(source);
    }

    /*!
     * \brief Returns the exact solution at a position.
     * \param globalPos The center of the finite volume for which it is to be set.
     */
    static PrimaryVariables exact(const GlobalPosition& globalPos)
    {
        const auto x = globalPos[0];
        const auto y = globalPos[1];

        using std::cos;
        using std::sin;

        static const auto periodLength = getParam<Scalar>("Problem.ExactSolPeriodLength");
        const auto preFactorArg = periodLength * M_PI;
        const auto u = sin(preFactorArg * x) * cos(preFactorArg * y);

        return PrimaryVariables(u);
    }
};

} // end namespace Dumux

#endif
