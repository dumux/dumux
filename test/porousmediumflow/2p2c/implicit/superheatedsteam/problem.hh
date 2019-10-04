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

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/method.hh>

#include <dumux/material/components/air.hh>
#include <dumux/material/fluidsystems/h2oair.hh>
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

// Set the primary variable combination for the model
template<class TypeTag>
struct Formulation<TypeTag, TTag::WaterAir>
{
    static constexpr auto value = TwoPFormulation::p1s0;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::WaterAir> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::WaterAir> { using type = WaterAirProblem<TypeTag>; };

// Set the wetting phase
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::WaterAir> { using type = FluidSystems::H2OAir<GetPropType<TypeTag, Properties::Scalar>>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::WaterAir>
{
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = WaterAirSpatialParams<FVGridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::WaterAir> { static constexpr bool value = true; };
} // end namespace Dumux

/*!
 * \ingroup TwoPTwoCModel
 * \brief Non-isothermal gas injection problem where superheated steam
 *        is injected into a partially water saturated medium.
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
    using GridView = GetPropType<TypeTag, Properties::GridView>;
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
        contiAirEqIdx = Indices::conti0EqIdx + FluidSystem::AirIdx
    };

    // phase presence
    enum { bothPhases = Indices::bothPhases };
    // component index
    enum { AirIdx = FluidSystem::AirIdx, H2OIdx = FluidSystem::H2OIdx };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;

    //! Property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = ModelTraits::useMoles();

public:
    WaterAirProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        FluidSystem::init();

        name_ = getParam<std::string>("Problem.Name");

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
        if(globalPos[0] > 1. - eps_)
            bcTypes.setAllDirichlet();
        else
            bcTypes.setAllNeumann();

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

        // we inject superheated steam from the left boundary
        if (isInjectionArea_(globalPos))
        {
            values[contiAirEqIdx] = useMoles ? -1e-5/FluidSystem::molarMass(AirIdx) : -1e-5; // kg/(m^2*s) or mole/(m^2*s)
            values[contiH2OEqIdx] = useMoles ? -1e-1/FluidSystem::molarMass(H2OIdx) : -1e-1; // kg/(m^2*s) or mole/(m^2*s)

            // note: energy equation is always formulated in terms of mass specific quantities, not per mole
            values[energyEqIdx] = -1.8e5;
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

        return priVars;
    }

private:
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars.setState(bothPhases);
        priVars[pressureIdx] = 1e5;
        priVars[switchIdx] = 0.1;
        priVars[temperatureIdx] = 293.15;
        return priVars;
    }

    bool isInjectionArea_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] < eps_;
    }

    static constexpr Scalar eps_ = 1e-2;
    std::string name_;
};

} // end namespace Dumux

#endif
