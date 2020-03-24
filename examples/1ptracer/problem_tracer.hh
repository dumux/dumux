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

#ifndef DUMUX_TRACER_TEST_PROBLEM_HH
#define DUMUX_TRACER_TEST_PROBLEM_HH

// ## The file `problem_tracer.hh`
//
// This file contains the __problem class__ which defines the initial and boundary
// conditions for the tracer transport simulation.
//
// ### Include files
// This header contains the porous medium problem class that this class is derived from:
#include <dumux/porousmediumflow/problem.hh>
// This header contains the class that specifies all spatially variable parameters
// related to this problem.
#include "spatialparams_tracer.hh"

// ### The problem class
// We enter the problem class where all necessary boundary conditions and initial conditions are set for our simulation.
// As this is a porous medium flow problem, we inherit from the base class `PorousMediumFlowProblem`.
namespace Dumux {

template <class TypeTag>
class TracerTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    // We use convenient declarations that we derive from the property system.
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    // We create a convenience bool stating whether mole or mass fractions are used
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();
    // We create additional convenience integers to make dimWorld and numComponents available in the problem
    static constexpr int dimWorld = GridView::dimensionworld;
    static const int numComponents = FluidSystem::numComponents;

public:
    // This is the constructor of our problem class:
    TracerTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        // We print to the terminal whether mole or mass fractions are used
        if(useMoles)
            std::cout<<"problem uses mole fractions" << '\n';
        else
            std::cout<<"problem uses mass fractions" << '\n';
    }

    // We define the type of boundary conditions depending on the location.
    // All boundaries are set to a neumann-type flow boundary condition.
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    // We specify the initial conditions for the primary variable (tracer concentration) depending on the location.
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables initialValues(0.0);

        // The initial contamination is located at the bottom of the domain:
        if (globalPos[1] < 0.1 + eps_)
        {
            // We chose a mole fraction of $`1e-9`$, but in case the mass fractions
            // are used by the model, we have to convert this value:
            if (useMoles)
                initialValues = 1e-9;
            else
                initialValues = 1e-9*FluidSystem::molarMass(0)
                                    /this->spatialParams().fluidMolarMassAtPos(globalPos);
        }
        return initialValues;
    }

    // We implement an outflow boundary on the top of the domain and prescribe zero-flux Neumann boundary conditions on all other boundaries.
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);
        const auto& volVars = elemVolVars[scvf.insideScvIdx()];
        const auto& globalPos = scvf.center();

        // This is the outflow boundary, where tracer is transported by advection with the given flux field.
        if (globalPos[dimWorld-1] > this->gridGeometry().bBoxMax()[dimWorld-1] - eps_)
        {
            values = this->spatialParams().volumeFlux(element, fvGeometry, elemVolVars, scvf)
                     * volVars.massFraction(0, 0) * volVars.density(0)
                     / scvf.area();
            assert(values>=0.0 && "Volume flux at outflow boundary is expected to have a positive sign");
        }

        // Prescribe zero-flux Neumann boundary conditions elsewhere
        else
            values = 0.0;

        return values;
    }

private:
    // We assign a private global variable for the epsilon:
    static constexpr Scalar eps_ = 1e-6;

// This is everything the tracer problem class contains.
};

// We leave the namespace Dumux here.
} // end namespace Dumux

#endif
