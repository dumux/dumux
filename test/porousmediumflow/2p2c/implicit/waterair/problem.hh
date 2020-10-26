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
 * \ingroup TwoPTwoCTests
 * \brief Non-isothermal gas injection problem where a gas (e.g. air)
 *        is injected into a fully water saturated medium.
 */

#ifndef DUMUX_WATER_AIR_PROBLEM_HH
#define DUMUX_WATER_AIR_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/boundarytypes.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/material/components/n2.hh>
#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/porousmediumflow/2p2c/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include "spatialparams.hh"

namespace Dumux {

template <class TypeTag>
class WaterAirProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct WaterAir { using InheritsFrom = std::tuple<TwoPTwoCNI>; };
struct WaterAirBox { using InheritsFrom = std::tuple<WaterAir, BoxModel>; };
struct WaterAirCCTpfa { using InheritsFrom = std::tuple<WaterAir, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::WaterAir> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::WaterAir> { using type = WaterAirProblem<TypeTag>; };

// Set the wetting phase
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::WaterAir> { using type = FluidSystems::H2ON2<GetPropType<TypeTag, Properties::Scalar>>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::WaterAir>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = WaterAirSpatialParams<GridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::WaterAir> { static constexpr bool value = true; };

// Enable caching
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::WaterAir> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::WaterAir> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::WaterAir> { static constexpr bool value = true; };
} // end namespace Dumux

/*!
 * \ingroup TwoPTwoCModel
 * \brief Non-isothermal gas injection problem where a gas (e.g. air)
 *        is injected into a fully water saturated medium.
 *
 * During buoyancy driven upward migration the gas passes a high temperature area.
 *
 * The domain is sized 40m times 40m in a depth of 1000m. The rectangular area
 * with the increased temperature (380K) starts at (20m, 1m) and ends at
 * (30m, 30m).
 *
 * For the mass conservation equation Neumann boundary conditions are used on
 * the top and on the bottom of the domain, while Dirichlet conditions
 * apply on the left and the right boundary.
 * For the energy conservation equation Dirichlet boundary conditions are applied
 * on all boundaries.
 *
 * Gas is injected at the bottom boundary from 15m to 25m at a rate of
 * 0.001kg/(s m), the remaining Neumann boundaries are no-flow boundaries.
 *
 * At the Dirichlet boundaries a hydrostatic pressure, a gas saturation of zero and
 * a geothermal temperature gradient of 0.03K/m are applied.
 *
 * The model is able to use either mole or mass fractions. The property useMoles
 * can be set to either true or false in the problem file. Make sure that the
 * according units are used in the problem set-up.
 * The default setting for useMoles is true.
 *
 * This problem uses the \ref TwoPTwoCModel and \ref NIModel model.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box2p2cni</tt> or
 * <tt>./test_cc2p2cni</tt>
 *  */
template <class TypeTag >
class WaterAirProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    // primary variable indices
    enum
    {
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,
        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx
    };

    // equation indices
    enum
    {
        contiH2OEqIdx = Indices::conti0EqIdx + FluidSystem::H2OIdx,
        contiN2EqIdx = Indices::conti0EqIdx + FluidSystem::N2Idx
    };

    // phase presence
    enum { wPhaseOnly = Indices::firstPhaseOnly };
    // component index
    enum { N2Idx = FluidSystem::N2Idx };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    //! Property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = ModelTraits::useMoles();

public:
    WaterAirProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        maxDepth_ = 1000.0; // [m]

        FluidSystem::init();

        name_ = getParam<std::string>("Problem.Name");
        useDirichlet_ = name_.find("buoyancy") != std::string::npos;

        //stating in the console whether mole or mass fractions are used
        if(useMoles)
            std::cout << "The problem uses mole-fractions" << std::endl;
        else
            std::cout << "The problem uses mass-fractions" << std::endl;

        this->spatialParams().plotMaterialLaw();
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        if(globalPos[0] > 40 - eps_ || globalPos[0] < eps_)
            bcTypes.setAllDirichlet();
        else
            bcTypes.setAllNeumann();

        if (useDirichlet_)
        {
            if (isInjectionArea_(globalPos))
                bcTypes.setAllDirichlet();
        }

        return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     *
     * \param globalPos The position for which the Dirichlet condition should be evaluated
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars = initial_(globalPos);

        if (useDirichlet_)
        {
            if (isInjectionArea_(globalPos))
            {
                priVars.setState(Indices::bothPhases);
                priVars[switchIdx] = 0.2;
                return priVars;
            }
        }

        return priVars;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     *
     * \param globalPos The position for which the Neumann condition should be evaluated
     * \return the mole/mass flux of each component. Negative values mean influx.
     *
     * The units must be according to using either mole or mass fractions (mole/(m^2*s) or kg/(m^2*s)).
     */
    NumEqVector neumannAtPos(const GlobalPosition &globalPos) const
    {
        NumEqVector values(0.0);

        // we inject pure gasious nitrogen at the initial condition temperature and pressure  from the bottom (negative values mean injection)
        if (isInjectionArea_(globalPos))
        {
            values[contiN2EqIdx] = useMoles ? -1e-3/FluidSystem::molarMass(N2Idx) : -1e-3; // kg/(m^2*s) or mole/(m^2*s)

            const auto initialValues = initial_(globalPos);
            const auto& fluidMatrixInteraction = this->spatialParams().fluidMatrixInteractionAtPos(globalPos);
            const auto pn = initialValues[pressureIdx] + fluidMatrixInteraction.endPointPc();
            const auto t = initialValues[temperatureIdx];

            // note: energy equation is always formulated in terms of mass specific quantities, not per mole
            values[energyEqIdx] = -1e-3*Components::N2<Scalar>::gasEnthalpy(t, pn);
        }

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
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        auto priVars = initial_(globalPos);
        // initially there is a heated lens in the domain
        if (globalPos[0] > 20 - eps_ && globalPos[0] < 30 + eps_ && globalPos[1] < 30 + eps_)
            priVars[temperatureIdx] = 380.0;

        return priVars;
    }

private:
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars.setState(wPhaseOnly);
        Scalar densityW = 1000.0;
        priVars[pressureIdx] = 1e5 + (maxDepth_ - globalPos[1])*densityW*9.81;
        priVars[switchIdx] = 0.0;
        priVars[temperatureIdx] = initialTemperatureProfile_(globalPos);
        return priVars;
    }

    Scalar initialTemperatureProfile_(const GlobalPosition &globalPos) const
    { return 283.0 + (maxDepth_ - globalPos[1])*0.03; }

    bool isInjectionArea_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] > 14.8 - eps_ && globalPos[0] < 25.2 + eps_ && globalPos[1] < eps_;
    }

    Scalar maxDepth_;
    static constexpr Scalar eps_ = 1e-2;
    std::string name_;
    bool useDirichlet_;
};

} // end namespace Dumux

#endif
