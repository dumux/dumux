// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Static single-phase flow on a 1d network embedded in 3d, used to validate the
 *        parallel (distributed, overlapping) FoamGrid against the sequential solution.
 */
#ifndef DUMUX_TEST_1P_NETWORK_PARALLEL_PROBLEM_HH
#define DUMUX_TEST_1P_NETWORK_PARALLEL_PROBLEM_HH

#include <limits>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

template<class TypeTag>
class OnePNetworkProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    OnePNetworkProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        // Neumann on the boundary part with x < NeumannXMax, Dirichlet on the rest. Default is a
        // large negative value, i.e. Dirichlet everywhere (well-posed on the multi-component
        // network). For connected grids a positive value enables a mixed BC. Position-based, so it
        // is identical on the full (reference) and the distributed grid.
        neumannXMax_ = getParam<Scalar>("Problem.NeumannXMax", std::numeric_limits<Scalar>::lowest());
    }

    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes bcTypes;
        if (globalPos[0] < neumannXMax_)
            bcTypes.setAllNeumann();
        else
            bcTypes.setAllDirichlet();
        return bcTypes;
    }

    //! a position-dependent pressure so the solution has a non-trivial gradient
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(1.0e5);
        values[0] = 1.0e5 + 1.0e9 * globalPos[0];
        return values;
    }

    //! prescribed mass flux on the Neumann ends [kg/(m^2 s)]
    NumEqVector neumannAtPos(const GlobalPosition& /*globalPos*/) const
    { return NumEqVector(1.0e-2); }

    //! no source
    PrimaryVariables sourceAtPos(const GlobalPosition& /*globalPos*/) const
    { return PrimaryVariables(0.0); }

    PrimaryVariables initialAtPos(const GlobalPosition& /*globalPos*/) const
    { return PrimaryVariables(1.0e5); }

private:
    Scalar neumannXMax_;
};

} // end namespace Dumux

#endif
