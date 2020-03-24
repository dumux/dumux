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

#ifndef DUMUX_ONEP_TEST_PROBLEM_HH
#define DUMUX_ONEP_TEST_PROBLEM_HH

// ## The file `problem_1p.hh`
//
// This file contains the __problem class__ which defines the initial and boundary
// conditions for the single-phase flow simulation.
//
// ### Include files
// This header contains the porous medium problem class that this class is derived from:
#include <dumux/porousmediumflow/problem.hh>
// This header contains the class that specifies all spatially variable parameters
// related to this problem.
#include "spatialparams_1p.hh"

// ### The problem class
// We enter the problem class where all necessary boundary conditions and initial conditions are set for our simulation.
// As this is a porous medium flow problem, we inherit from the base class `PorousMediumFlowProblem`.
namespace Dumux {

template<class TypeTag>
class OnePTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    // We use convenient declarations that we derive from the property system.
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;

    static constexpr int dimWorld = GridView::dimensionworld;

public:
    // This is the constructor of our problem class:
    OnePTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry) {}

    // First, we define the type of boundary conditions depending on the location. Two types of boundary conditions
    // can be specified: Dirichlet or Neumann boundary condition. On a Dirichlet boundary, the values of the
    // primary variables need to be fixed. On a Neumann boundary condition, values for derivatives need to be fixed.
    // Mixed boundary conditions (different types for different equations on the same boundary) are not accepted for
    // cell-centered finite volume schemes.
    BoundaryTypes boundaryTypes(const Element &element,
                                const SubControlVolumeFace &scvf) const
    {
        BoundaryTypes values;
        // we retrieve the global position, i.e. the vector with the global coordinates,
        // of the integration point on the boundary sub-control volume face `scvf`
        const auto globalPos = scvf.ipGlobal();
        // we define a small epsilon value
        Scalar eps = 1.0e-6;
        // We specify Dirichlet boundaries on the top and bottom of our domain:
        if (globalPos[dimWorld-1] < eps || globalPos[dimWorld-1] > this->gridGeometry().bBoxMax()[dimWorld-1] - eps)
            values.setAllDirichlet();
        // The top and bottom of our domain are Neumann boundaries:
        else
            values.setAllNeumann();

        return values;
    }

    // Second, we specify the values for the Dirichlet boundaries. We need to fix values of our primary variable
    PrimaryVariables dirichlet(const Element &element,
                               const SubControlVolumeFace &scvf) const
    {
        // we retreive again the global position
        const auto& pos = scvf.ipGlobal();
        PrimaryVariables values(0);
        // and assign pressure values in [Pa] according to a pressure gradient to 1e5 Pa at the top and 1.1e5 Pa at the bottom.
        values[0] = 1.0e+5*(1.1 - pos[dimWorld-1]*0.1);
        return values;
    }

    // We need to specify a constant temperature for our isothermal problem.
    // Fluid properties that depend on temperature will be calculated with this value.
    Scalar temperature() const
    {
        return 283.15; // 10°C
    }

    // This is everything the one phase problem class contains.
};

// We leave the namespace Dumux.
} // end namespace Dumux
#endif
