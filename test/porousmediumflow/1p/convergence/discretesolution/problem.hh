// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
    {
        periodLength_ = getParam<Scalar>("Problem.ExactSolPeriodLength");
        sourceIntegrationOrder_ = getParam<Scalar>("Problem.SourceIntegrationOrder");
    }

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
        return exact(globalPos, periodLength_);
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
        const auto &k = this->spatialParams().permeabilityAtPos(scv.center());

        using std::cos;
        using std::sin;

        const auto eg = element.geometry();
        const auto rule = Dune::QuadratureRules<Scalar, GridView::dimension>::rule(eg.type(), sourceIntegrationOrder_);

        Scalar source = 0.0;
        for (auto qp : rule)
        {
            const auto p = eg.global(qp.position());
            const auto x = p[0];
            const auto y = p[1];

            const auto sineTerm = sin(periodLength_*x);
            const auto cosTerm = cos(periodLength_*y);
            const auto secondDeriv = -1.0*periodLength_*periodLength_*sineTerm*cosTerm;

            // derivative in x and y are identical
            source -= 2.0*k*secondDeriv*qp.weight()*eg.integrationElement(qp.position());
        }

        source /= eg.volume();
        return NumEqVector(source);
    }

    /*!
     * \brief Returns the exact solution at a position.
     * \param globalPos The center of the finite volume for which it is to be set.
     */
    static PrimaryVariables exact(const GlobalPosition& globalPos, const Scalar periodLength)
    {
        const auto x = globalPos[0];
        const auto y = globalPos[1];

        using std::cos;
        using std::sin;

        const auto u = sin(periodLength*x)*cos(periodLength*y);
        return PrimaryVariables(u);
    }

private:
    Scalar periodLength_;
    Scalar sourceIntegrationOrder_;
};

} // end namespace Dumux

#endif
