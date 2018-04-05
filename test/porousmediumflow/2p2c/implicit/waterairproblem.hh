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
/*!
 * \file
 * \ingroup TwoPTwoCTests
 * \brief Non-isothermal gas injection problem where a gas (e.g. air)
 *        is injected into a fully water saturated medium.
 */
#ifndef DUMUX_WATER_AIR_PROBLEM_HH
#define DUMUX_WATER_AIR_PROBLEM_HH

#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/box/properties.hh>
#include <dumux/discretization/methods.hh>

#include <dumux/material/components/n2.hh>
#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/porousmediumflow/2p2c/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include "waterairspatialparams.hh"

namespace Dumux
{
/*!
 * \ingroup TwoPTwoCTests
 * \brief Non-isothermal gas injection problem where a gas (e.g. air)
 *        is injected into a fully water saturated medium.
 */
template <class TypeTag>
class WaterAirProblem;

namespace Properties
{
NEW_TYPE_TAG(WaterAirTypeTag, INHERITS_FROM(TwoPTwoCNI, WaterAirSpatialParams));
NEW_TYPE_TAG(WaterAirBoxTypeTag, INHERITS_FROM(BoxModel, WaterAirTypeTag));
NEW_TYPE_TAG(WaterAirCCTpfaTypeTag, INHERITS_FROM(CCTpfaModel, WaterAirTypeTag));

// Set the grid type
SET_TYPE_PROP(WaterAirTypeTag, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(WaterAirTypeTag, Problem, WaterAirProblem<TypeTag>);

// Set the wetting phase
SET_TYPE_PROP(WaterAirTypeTag, FluidSystem, FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), true>);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(WaterAirTypeTag, UseMoles, true);
}


/*!
 * \ingroup TwoPTwoCModel
 * \ingroup ImplicitTestProblems
 * \brief Non-isothermal gas injection problem where a gas (e.g. air)
 *        is injected into a fully water saturated medium. During
 *        buoyancy driven upward migration the gas passes a high
 *        temperature area.
 *
 * The domain is sized 40 m times 40 m in a depth of 1000 m. The rectangular area
 * with the increased temperature (380 K) starts at (20 m, 1 m) and ends at
 * (30 m, 30 m).
 *
 * For the mass conservation equation neumann boundary conditions are used on
 * the top and on the bottom of the domain, while dirichlet conditions
 * apply on the left and the right boundary.
 * For the energy conservation equation dirichlet boundary conditions are applied
 * on all boundaries.
 *
 * Gas is injected at the bottom boundary from 15 m to 25 m at a rate of
 * 0.001 kg/(s m), the remaining neumann boundaries are no-flow
 * boundaries.
 *
 * At the dirichlet boundaries a hydrostatic pressure, a gas saturation of zero and
 * a geothermal temperature gradient of 0.03 K/m are applied.
 *
 * The model is able to use either mole or mass fractions. The property useMoles can be set to either true or false in the
 * problem file. Make sure that the according units are used in the problem setup. The default setting for useMoles is true.
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

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
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
        conti0EqIdx = Indices::conti0EqIdx,
        contiNEqIdx = Indices::contiNEqIdx
    };

    // phase presence
    enum { wPhaseOnly = Indices::wPhaseOnly };
    // component index
    enum { nCompIdx = FluidSystem::nCompIdx };

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;

    //! property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = ModelTraits::useMoles();

public:
    /*!
     * \brief The constructor.
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    WaterAirProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        maxDepth_ = 1000.0; // [m]

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
     * \param values The boundary types for the conservation equations
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        if(globalPos[0] > 40 - eps_ || globalPos[0] < eps_)
            bcTypes.setAllDirichlet();
        else
            bcTypes.setAllNeumann();

        return bcTypes;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Dirichlet boundary segment.
     *
     * \param globalPos The position for which the Dirichlet condition should be evaluated
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return initial_(globalPos);
    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann boundary segment.
     *
     * \param globalPos The position for which the Neumann condition should be evaluated
     *
     * \return the mole/mass flux of each component. Negative values mean influx.
     *
     * The units must be according to using either mole or mass fractions. (mole/(m^2*s) or kg/(m^2*s))
     */
    NumEqVector neumannAtPos(const GlobalPosition &globalPos) const
    {
        NumEqVector values(0.0);

        // we inject pure gasious nitrogen at the initial condition temperature and pressure  from the bottom (negative values mean injection)
        if (globalPos[0] > 14.8 - eps_ && globalPos[0] < 25.2 + eps_ && globalPos[1] < eps_)
        {
            values[contiNEqIdx] = useMoles ? -1e-3/FluidSystem::molarMass(nCompIdx) : -1e-3; // kg/(m^2*s) or mole/(m^2*s)

            const auto initialValues = initial_(globalPos);
            const auto& mParams = this->spatialParams().materialLawParamsAtPos(globalPos);
            using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
            const auto pn = initialValues[pressureIdx] + MaterialLaw::endPointPc(mParams);
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
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The initial values for the primary variables
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

    Scalar maxDepth_;
    static constexpr Scalar eps_ = 1e-2;
    std::string name_;
};

} //end namespace Dumux

#endif
