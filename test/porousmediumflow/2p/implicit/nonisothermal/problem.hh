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
 * \ingroup TwoPTests
 * \brief Non-isothermal gas injection problem where a gas (e.g. air) is injected into a fully
 *        water saturated medium.
 *
 * During buoyancy driven upward migration the gas passes a high temperature area.
 */

#ifndef DUMUX_INJECTION_PROBLEM_2PNI_HH
#define DUMUX_INJECTION_PROBLEM_2PNI_HH

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/boundarytypes.hh>

#include <dumux/porousmediumflow/2p/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/components/n2.hh>

// use the spatial parameters as the injection problem of the 2p2c test program
#include <test/porousmediumflow/2p2c/implicit/injection/spatialparams.hh>

#ifndef GRIDTYPE // default to yasp grid if not provided by CMake
#define GRIDTYPE Dune::YaspGrid<2>
#endif

namespace Dumux {

//! Forward declaration of the problem class
template <class TypeTag> class InjectionProblem2PNI;

namespace Properties {
// Create new type tags
namespace TTag {
struct Injection2PNITypeTag { using InheritsFrom = std::tuple<TwoPNI>; };
struct InjectionBox2PNITypeTag { using InheritsFrom = std::tuple<Injection2PNITypeTag, BoxModel>; };
struct InjectionCC2PNITypeTag { using InheritsFrom = std::tuple<Injection2PNITypeTag, CCTpfaModel>; };
} // end namespace TTag

// Obtain grid type from COMPILE_DEFINITIONS
template<class TypeTag>
struct Grid<TypeTag, TTag::Injection2PNITypeTag> { using type = GRIDTYPE; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Injection2PNITypeTag> { using type = InjectionProblem2PNI<TypeTag>; };

// Use the same fluid system as the 2p2c injection problem
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Injection2PNITypeTag> { using type = FluidSystems::H2ON2<GetPropType<TypeTag, Properties::Scalar>, FluidSystems::H2ON2DefaultPolicy</*fastButSimplifiedRelations=*/true>>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Injection2PNITypeTag>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = InjectionSpatialParams<GridGeometry, Scalar>;
};

} // namespace Properties

/*!
 * \ingroup TwoPTests
 * \brief Non-isothermal gas injection problem where a gas (e.g. air) is injected into a fully
 *        water saturated medium.
 *
 * During buoyancy driven upward migration the gas passes a high temperature area.
 *
 * The domain is sized 60m times 40m. The rectangular area with the increased
 * temperature (380K) starts at (20m, 5m) and ends at (30m, 35m).
 *
 * For the mass conservation equation Neumann boundary conditions are used on
 * the top, on the bottom and on the right of the domain, while Dirichlet conditions
 * are applied on the left boundary.
 * For the energy conservation equation Dirichlet boundary conditions are applied
 * on all boundaries.
 *
 * Gas is injected at the right boundary from 5m to 15m at a rate of
 * 0.001kg/(s m), the remaining Neumann boundaries are no-flow boundaries.
 *
 * At the Dirichlet boundaries a hydrostatic pressure, a gas saturation of zero and
 * a geothermal temperature gradient of 0.03K/m are applied.
 *
 * This problem uses the \ref TwoPModel and \ref NIModel model.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box2pni -parameterFile test_box2pni.input</tt> or
 * <tt>./test_cc2pni -parameterFile test_cc2pni.input</tt>
 */
template<class TypeTag>
class InjectionProblem2PNI : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    enum
    {
        //! Primary variable indices
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        temperatureIdx = Indices::temperatureIdx,

        //! Equation indices
        contiN2EqIdx = Indices::conti0EqIdx + FluidSystem::N2Idx,
        energyEqIdx = Indices::energyEqIdx,

        //! Phase indices
        wPhaseIdx = FluidSystem::H2OIdx,
        nPhaseIdx = FluidSystem::N2Idx,

        // world dimension
        dimWorld = GridView::dimensionworld
    };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;

public:
    InjectionProblem2PNI(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        maxDepth_ = 2700.0; // [m]

        // initialize the tables of the fluid system
        FluidSystem::init(/*tempMin=*/273.15,
                          /*tempMax=*/423.15,
                          /*numTemp=*/50,
                          /*pMin=*/0.0,
                          /*pMax=*/30e6,
                          /*numP=*/300);

        name_ = getParam<std::string>("Problem.Name");
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string name() const
    { return name_; }

    /*!
     * \brief Returns the source term
     *
     * \param globalPos The global position
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        NumEqVector values(0.0);
        values = 0;
        return values;
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        if (globalPos[0] < eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     *
     * \param globalPos The global position
     */
     PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        Scalar densityW = 1000.0;
        values[pressureIdx] = 1e5 + (maxDepth_ - globalPos[dimWorld-1])*densityW*9.81;
        values[saturationIdx] = 0.0;
        values[temperatureIdx] = 283.0 + (maxDepth_ - globalPos[dimWorld-1])*0.03;
        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     *
     * \param globalPos The global position
     *
     * The \a values store the mass flux of each phase normal to the boundary.
     * Negative values indicate an inflow.
     */
    NumEqVector neumannAtPos(const GlobalPosition &globalPos) const
    {
        NumEqVector values(0.0);

        if (globalPos[1] < 13.75 + eps_ && globalPos[1] > 6.875 - eps_)
        {
            // inject air. negative values mean injection
            values[contiN2EqIdx] = -1e-3; // kg/(s*m^2)

            // compute enthalpy flux associated with this injection [(J/(kg*s)]
            using FluidState = GetPropType<TypeTag, Properties::FluidState>;
            FluidState fs;

            const auto initialValues = initialAtPos(globalPos);
            fs.setPressure(wPhaseIdx, initialValues[pressureIdx]);
            fs.setPressure(nPhaseIdx, initialValues[pressureIdx]); // assume pressure equality here
            fs.setTemperature(wPhaseIdx,initialValues[temperatureIdx]);
            fs.setTemperature(nPhaseIdx,initialValues[temperatureIdx]);

            // energy flux is mass flux times specific enthalpy
            values[energyEqIdx] = values[contiN2EqIdx]*FluidSystem::enthalpy(fs, nPhaseIdx);
        }

        return values;
    }

    // \}


    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the initial values for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        Scalar densityW = 1000.0;
        values[pressureIdx] = 1e5 + (maxDepth_ - globalPos[1])*densityW*9.81;
        values[saturationIdx] = 0.0;
        values[temperatureIdx] = 283.0 + (maxDepth_ - globalPos[1])*0.03;

        if (globalPos[0] > 21.25 - eps_ && globalPos[0] < 28.75 + eps_ && globalPos[1] > 6.25 - eps_ && globalPos[1] < 33.75 + eps_)
            values[temperatureIdx] = 380;

        return values;
    }
    // \}

private:
    Scalar maxDepth_;
    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;
};

} // end namespace Dumux

#endif
