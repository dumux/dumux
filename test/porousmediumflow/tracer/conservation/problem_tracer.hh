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
 * \ingroup TracerTests
 * \brief  problem for a tracer conservation test
 */
#ifndef DUMUX_TEST_PROBLEM_TRACER_HH
#define DUMUX_TEST_PROBLEM_TRACER_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/discretization/cctpfa.hh>

#include <dumux/porousmediumflow/tracer/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/material/fluidsystems/base.hh>

#include "spatialparams_tracer.hh"

namespace Dumux {
/*!
 * \ingroup TracerTests
 * \brief problem for a tracer conservation test
 */
template <class TypeTag>
class TracerConservationTestProblem;

namespace Properties {
//Create new type tags
namespace TTag {
struct TracerConservationTest { using InheritsFrom = std::tuple<Tracer>; };
struct TracerConservationTestTpfa { using InheritsFrom = std::tuple<TracerConservationTest, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TracerConservationTest> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TracerConservationTest> { using type = TracerConservationTestProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TracerConservationTest>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = TracerConservationTestSpatialParams<GridGeometry, Scalar>;
};

//! A simple fluid system with one tracer component
template<class TypeTag>
class TracerFluidSystem : public FluidSystems::Base<GetPropType<TypeTag, Properties::Scalar>,
                                                                TracerFluidSystem<TypeTag>>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

public:
    //! If the fluid system only contains tracer components
    static constexpr bool isTracerFluidSystem()
    { return true; }

    //! No component is the main component
    static constexpr int getMainComponent(int phaseIdx)
    { return -1; }

    //! The number of components
    static constexpr int numComponents = 1;

    //! Human readable component name (for vtk output)
    static std::string componentName(int compIdx)
    { return "Cs137"; }

    //! Molar mass in kg/mol
    static Scalar molarMass(int compIdx)
    { return 0.1; }

    //! Binary diffusion coefficient
    //! (might depend on spatial parameters like pressure / temperature)
    static Scalar binaryDiffusionCoefficient(unsigned int compIdx,
                                             const Problem& problem,
                                             const Element& element,
                                             const SubControlVolume& scv)
    {
        static const Scalar D = getParam<Scalar>("Problem.BinaryDiffusionCoefficient");
        return D;
    }
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TracerConservationTest> { using type = TracerFluidSystem<TypeTag>; };

} // end namespace Properties

template <class TypeTag>
class TracerConservationTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();

public:
    TracerConservationTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        // stating in the console whether mole or mass fractions are used
        if(useMoles)
            std::cout<<"problem uses mole fractions" << '\n';
        else
            std::cout<<"problem uses mass fractions" << '\n';

        useDirichlet_ = getParam<bool>("Problem.UseDirichletForTracer", false);
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        if (useDirichlet_ && onUpperBoundary_(globalPos))
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        return values;
    }

    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector values(0.0);

        if (onUpperBoundary_(globalPos))
        {
            if (useMoles)
                values = -1e-10;
            else
                values = -1e-10*FluidSystem::molarMass(0)/this->spatialParams().fluidMolarMass(globalPos);
        }

        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;

        if (onUpperBoundary_(globalPos))
            values = 1e-9;
        else
            values = initialAtPos(globalPos);

        return values;
    }

    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables initialValues(0.0);

        return initialValues;
    }

    //! Helper function for mass balance evaluations
    template <class TracerSolutionVector, class GridVariables>
    Scalar printTracerMass(const TracerSolutionVector& sol,
                           const GridVariables& gridVars,
                           Scalar dt,
                           Scalar inflowMass)

    {
        // compute the mass in the entire domain to make sure the tracer is conserved
        Scalar tracerMass = 0.0;

        const auto& gg = this->gridGeometry();
        for (const auto& element : elements(gg.gridView()))
        {
            auto fvGeometry = localView(gg);
            fvGeometry.bind(element);

            auto elemVolVars = localView(gridVars.curGridVolVars());
            elemVolVars.bind(element, fvGeometry, sol);

            auto elemFluxVarsCache = localView(gridVars.gridFluxVarsCache());
            elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

            for (auto&& scv : scvs(fvGeometry))
            {
                const auto& volVars = elemVolVars[scv];
                if(useMoles)
                {
                    tracerMass += volVars.moleFraction(/*phaseIdx*/0, /*compIdx*/0) * volVars.molarDensity(/*phaseIdx*/0)
                                * scv.volume() * volVars.saturation(/*phaseIdx*/0) * volVars.porosity() * volVars.extrusionFactor();
                }
                else
                    tracerMass += volVars.massFraction(/*phaseIdx*/0, /*compIdx*/0) * volVars.density(/*phaseIdx*/0)
                                * scv.volume() * volVars.saturation(/*phaseIdx*/0) * volVars.porosity() * volVars.extrusionFactor();
            }
            for (auto&& scvf : scvfs(fvGeometry))
            {
                if (scvf.boundary())
                {
                    const auto& volVars = elemVolVars[scvf.insideScvIdx()];
                    if (useDirichlet_)
                    {
                        const auto moleFraction = dirichletAtPos(scvf.ipGlobal());
                        inflowMass -= this->spatialParams().volumeFlux(element, fvGeometry, elemVolVars, scvf)
                                    * moleFraction
                                    * volVars.molarDensity(/*phaseIdx*/0)
                                    * scvf.area() * volVars.extrusionFactor()
                                    * dt;
                    }
                    else
                        inflowMass -= neumannAtPos(scvf.ipGlobal())[0]
                                    * scvf.area() * volVars.extrusionFactor()
                                    * dt;
                }
            }
        }

        std::cout << "\033[1;31m" << tracerMass << " moles tracer are inside the domain. \033[0m" << '\n';
        std::cout << "\033[1;31m" << inflowMass << " moles should have flown into the domain. \033[0m" << '\n';
        std::cout << "\033[1;31m" << inflowMass - tracerMass << " moles are missing. \033[0m" << '\n';

        return inflowMass;
    }

private:
    static constexpr Scalar eps_ = 1e-6;

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_;
    }

    bool useDirichlet_;
};

} //end namespace Dumux

#endif
