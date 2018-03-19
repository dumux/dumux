// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
/**
 * \file
 * \ingroup TracerTests
 * \brief Definition of a problem, for the tracer problem:
 * A rotating velocity field mixes a tracer band in a porous groundwater reservoir.
 */
#ifndef DUMUX_TRACER_TEST_PROBLEM_HH
#define DUMUX_TRACER_TEST_PROBLEM_HH

#include <fstream>

#include <dumux/discretization/cellcentered/mpfa/properties.hh>
#include <dumux/porousmediumflow/tracer/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/fluidsystems/base.hh>

#include "tracertestspatialparams.hh"

namespace Dumux
{
/**
 * \ingroup TracerTests
 * \brief Definition of a problem, for the tracer problem:
 * A rotating velocity field mixes a tracer band in a porous groundwater reservoir.
 */
template <class TypeTag>
class TracerTestProblem;

namespace Properties
{
NEW_TYPE_TAG(TracerTestProblem, INHERITS_FROM(Tracer));
NEW_TYPE_TAG(TracerTestCCProblem, INHERITS_FROM(CCMpfaModel, TracerTestProblem));

// enable caching
SET_BOOL_PROP(TracerTestProblem, EnableGridVolumeVariablesCache, true);
SET_BOOL_PROP(TracerTestProblem, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(TracerTestProblem, EnableFVGridGeometryCache, true);

// Set the grid type
SET_TYPE_PROP(TracerTestProblem, Grid, Dune::UGGrid<3>);

// Set the problem property
SET_TYPE_PROP(TracerTestProblem, Problem, TracerTestProblem<TypeTag>);

// Set the spatial parameters
SET_TYPE_PROP(TracerTestProblem, SpatialParams, TracerTestSpatialParams<TypeTag>);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(TracerTestProblem, UseMoles, false);
SET_BOOL_PROP(TracerTestCCProblem, EnableMolecularDiffusion, false);
SET_BOOL_PROP(TracerTestCCProblem, SolutionDependentAdvection, false);

//! A simple fluid system with one tracer component
template<class TypeTag>
class TracerFluidSystem : public FluidSystems::BaseFluidSystem<typename GET_PROP_TYPE(TypeTag, Scalar),
                                                               TracerFluidSystem<TypeTag>>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

public:
    //! If the fluid system only contains tracer components
    static constexpr bool isTracerFluidSystem() { return true; }
    //! No component is the main component
    static constexpr int getMainComponent(int phaseIdx) { return -1; }
    //! The number of components
    static constexpr int numComponents = 1;
    //! Human readable component name (index compIdx) (for vtk output)
    static std::string componentName(int compIdx) { return "tracer_" + std::to_string(compIdx); }
    //! Human readable phase name (index phaseIdx) (for velocity vtk output)
    static std::string phaseName(int phaseIdx = 0) { return "Groundwater"; }
    //! Molar mass in kg/mol of the component with index compIdx
    static constexpr Scalar molarMass(unsigned int compIdx) { return 1.; }
    //! binary diffusion coefficient
    //! (might depend on spatial parameters like pressure / temperature)
    static constexpr Scalar binaryDiffusionCoefficient(unsigned int compIdx,
                                                       const Problem& problem,
                                                       const Element& element,
                                                       const SubControlVolume& scv)
    { return 0.0; }
};

SET_TYPE_PROP(TracerTestProblem, FluidSystem, TracerFluidSystem<TypeTag>);

} // end namespace Properties


/*!
 * \ingroup TracerModel
 * \ingroup ImplicitTestProblems
 *
 * \brief Definition of a problem, for the tracer problem:
 * A lens of contaminant tracer is diluted by diffusion and a base groundwater flow
 *
 * This problem uses the \ref TracerModel model.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_boxtracer -ParameterFile ./test_boxtracer.input</tt> or
 * <tt>./test_cctracer -ParameterFile ./test_cctracer.input</tt>
 */
template <class TypeTag>
class TracerTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);

    //! property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    TracerTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeom) : ParentType(fvGridGeom)
    {
        boundaryMoleFrac_ = getParam<Scalar>("Problem.BoundaryMoleFrac");
        outputFileName_ = getParam<std::string>("OutputFile.Name");

        // state in the console whether mole or mass fractions are used
        if(useMoles) std::cout<<"problem uses mole fractions" << '\n';
        else std::cout<<"problem uses mass fractions" << '\n';
    }

    //! writes the mass distribution in the layers to output file
    void writeMassDistribution(Scalar t,
                               const SolutionVector& sol,
                               const FVGridGeometry& fvGridGeometry,
                               const GridVariables& gridVars) const
    {
        std::ofstream file(outputFileName_, std::ofstream::out | std::ofstream::app);

        Scalar massUpperLayer = 0.0;
        Scalar massLowerLayer = 0.0;
        Scalar massFracture = 0.0;
        Scalar massInflow = 0.0;
        Scalar massOutflow = 0.0;

        using GridCreator = typename GET_PROP_TYPE(TypeTag, GridCreator);
        using SpatialParams = typename GET_PROP_TYPE(TypeTag, SpatialParams);
        for (const auto& element : elements(fvGridGeometry.gridView()))
        {
            const auto eIdx = fvGridGeometry.elementMapper().index(element);
            const auto marker = GridCreator::getElementDomainMarker(eIdx);

            using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
            ElementSolutionVector elemSol(element, sol, fvGridGeometry);
            auto fvGeometry = localView(fvGridGeometry);
            auto elemVolVars = localView(gridVars.curGridVolVars());

            fvGeometry.bindElement(element);
            elemVolVars.bindElement(element, fvGeometry, sol);

            // obtain reference to this elements variables
            const auto& vv = elemVolVars[eIdx];
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto mass = this->spatialParams().porosity(element, scv, elemSol)
                                  *FluidSystem::molarMass(0)
                                  *vv.moleFraction(0, 0)
                                  /this->spatialParams().fluidMolarMass()
                                  *this->spatialParams().fluidDensity()
                                  *scv.volume();

                if (marker == SpatialParams::LayerMarkers::lowerLayer)
                    massLowerLayer += mass;
                else if (marker == SpatialParams::LayerMarkers::upperLayer)
                    massUpperLayer += mass;
                else if (marker == SpatialParams::LayerMarkers::fracture)
                    massFracture += mass;
                else
                    DUNE_THROW(Dune::InvalidStateException, "Unknown layer marker");
            }

            // evaluate influx/outflux of tracer
            const auto c = element.geometry().center();
            if (c[2] > 90.0 || c[2] < 10.0)
            {
                for (const auto& scvf : scvfs(fvGeometry))
                {
                    const auto& ip = scvf.ipGlobal();
                    if (isOnInflowBoundary(ip))
                        massInflow += scvf.area()*neumann(element, fvGeometry, elemVolVars, scvf)[0];
                    else if (isOnOutflowBoundary(ip))
                        massOutflow += scvf.area()*neumann(element, fvGeometry, elemVolVars, scvf)[0];
                }
            }
        }

        file << t << "\t\t\t"
             << massInflow << "\t\t\t"
             << massOutflow << "\t\t\t"
             << massUpperLayer << "\t\t\t"
             << massFracture << "\t\t\t"
             << massLowerLayer << "\n";
        file.close();
    }

    //! Specifies which kind of boundary condition should be used
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    //! Evaluate a Neumann segment
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const SubControlVolumeFace& scvf) const
    {
        //! influx
        if (isOnInflowBoundary(scvf.ipGlobal()))
        {
            const auto influx = this->spatialParams().volumeFlux(scvf)
                                *this->spatialParams().fluidDensity()
                                /this->spatialParams().fluidMolarMass()
                                *boundaryMoleFrac_
                                *FluidSystem::molarMass(0)
                                /scvf.area();
            return NumEqVector(influx);
        }
        //! outflux
        else if (isOnOutflowBoundary(scvf.ipGlobal()))
        {
            const auto outFlux = this->spatialParams().volumeFlux(scvf)
                                 *this->spatialParams().fluidDensity()
                                 /this->spatialParams().fluidMolarMass()
                                 *elemVolVars[scvf.insideScvIdx()].moleFraction(0, 0)
                                 *FluidSystem::molarMass(0)/
                                 scvf.area();
            return NumEqVector(outFlux);
        }
        //! remaining boundaries
        return NumEqVector(0.0);
    }

    //! Evaluate the initial value for a control volume.
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const { return PrimaryVariables(0.0); }

    //! Returns if a given position is on inflow boundary
    bool isOnInflowBoundary(const GlobalPosition& globalPos) const { return globalPos[0] < 1.0e-6 && globalPos[2] > 90.0; }

    //! Returns if a given position is on outflow boundary
    bool isOnOutflowBoundary(const GlobalPosition& globalPos) const { return globalPos[1] < 1.0e-6 && globalPos[2] < 10.0; }

private:
    Scalar boundaryMoleFrac_;
    std::string outputFileName_;
};

} //end namespace Dumux

#endif
