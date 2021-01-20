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
 * \ingroup ThreePThreeCTests
 * \brief Isothermal NAPL infiltration problem: LNAPL contaminates
 *        the unsaturated and the saturated groundwater zone.
 */

#ifndef DUMUX_INFILTRATION_THREEPTHREEC_PROBLEM_HH
#define DUMUX_INFILTRATION_THREEPTHREEC_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/boundarytypes.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/3p3c/model.hh>
#include <dumux/material/fluidsystems/h2oairmesitylene.hh>

#include "spatialparams.hh"

namespace Dumux {

/*!
 * \ingroup ThreePThreeCTests
 * \brief Isothermal NAPL infiltration problem: LNAPL contaminates
 *        the unsaturated and the saturated groundwater zone.
 */
template <class TypeTag>
class InfiltrationThreePThreeCProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct InfiltrationThreePThreeC { using InheritsFrom = std::tuple<ThreePThreeC>; };
struct InfiltrationThreePThreeCBox { using InheritsFrom = std::tuple<InfiltrationThreePThreeC, BoxModel>; };
struct InfiltrationThreePThreeCCCTpfa { using InheritsFrom = std::tuple<InfiltrationThreePThreeC, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::InfiltrationThreePThreeC> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::InfiltrationThreePThreeC> { using type = InfiltrationThreePThreeCProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::InfiltrationThreePThreeC>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = InfiltrationThreePThreeCSpatialParams<GridGeometry, Scalar>;
};

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::InfiltrationThreePThreeC>
{ using type = FluidSystems::H2OAirMesitylene<GetPropType<TypeTag, Properties::Scalar>>; };
}

/*!
 * \ingroup ThreePThreeCTests
 * \brief Isothermal NAPL infiltration problem: LNAPL contaminates
 *        the unsaturated and the saturated groundwater zone.
 *
 * The 2D domain of this test problem is 500 m long and 10 m deep, where
 * the lower part represents a slightly inclined groundwater table, and the
 * upper part is the vadose zone.
 * A LNAPL (Non-Aqueous Phase Liquid which is lighter than water) infiltrates
 * (modelled with a Neumann boundary condition) into the vadose zone. Upon
 * reaching the water table, it spreads (since lighter than water) and migrates
 * on top of the water table in the direction of the slope.
 * On its way through the vadose zone, it leaves a trace of residually trapped
 * immobile NAPL, which can in the following dissolve and evaporate slowly,
 * and eventually be transported by advection and diffusion.
 *
 * Left and right boundaries are constant head boundaries (Dirichlet),
 * Top and bottom are Neumann boundaries, all no-flow except for the small
 * infiltration zone in the upper left part.
 *
 * This problem uses the \ref ThreePThreeCModel.
 *
 * This problem should typically be simulated for 30 days.
 * A good choice for the initial time step size is 60 s.
 * To adjust the simulation time it is necessary to edit the file test_box3p3c.input
 * or test_cc3p3c.input.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box3p3c test_box3p3c.input</tt> or
 * <tt>./test_cc3p3c test_cc3p3c.input</tt>
 *  */
template <class TypeTag >
class InfiltrationThreePThreeCProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    // copy some indices for convenience
    enum {
        pressureIdx = Indices::pressureIdx,
        switch1Idx = Indices::switch1Idx,
        switch2Idx = Indices::switch2Idx,

        // phase state
        wgPhaseOnly = Indices::wgPhaseOnly,

        contiWEqIdx = Indices::conti0EqIdx + FluidSystem::wCompIdx, //!< Index of the mass conservation equation for the water component
        contiNEqIdx = Indices::conti0EqIdx + FluidSystem::nCompIdx,//!< Index of the mass conservation equation for the contaminant component
        contiAEqIdx = Indices::conti0EqIdx + FluidSystem::gCompIdx,//!< Index of the mass conservation equation for the gas component

        // world dimension
        dimWorld = GridView::dimensionworld
    };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    InfiltrationThreePThreeCProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        temperature_ = 273.15 + 10.0; // -> 10 degrees Celsius
        FluidSystem::init(/*tempMin=*/temperature_ - 1,
                          /*tempMax=*/temperature_ + 1,
                          /*nTemp=*/3,
                          /*pressMin=*/0.8*1e5,
                          /*pressMax=*/3*1e5,
                          /*nPress=*/200);

        name_ = getParam<std::string>("Problem.Name");
    }

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * \param globalPos The global position
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    {
        return temperature_;
    }

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
        BoundaryTypes values;
        if(globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_)
            values.setAllDirichlet();
        else if(globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();
        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    { return initial_(globalPos); }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    NumEqVector neumannAtPos(const GlobalPosition &globalPos) const
    {
        NumEqVector values(0.0);

        // negative values for injection
        if ((globalPos[0] < 80.0 + eps_) && (globalPos[0] > 55.0 - eps_) && (globalPos[1] > 10.0 - eps_))
        {
            values[contiWEqIdx] = -0.0;
            //mole flow conversion to mass flow with molar mass M(Mesit.)=0,120 kg/mol --> 1.2e-4 kg/(sm)
            //the 3p3c model uses mole fractions
            values[contiNEqIdx] = -0.001;
            values[contiAEqIdx] = -0.0;
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
    { return initial_(globalPos); }

private:
    // internal method for the initial condition
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values.setState(wgPhaseOnly);

        Scalar y = globalPos[1];
        Scalar x = globalPos[0];
        Scalar sw, swr=0.12, sgr=0.03;

        if(y >(-1.E-3*x+5) - eps_)
        {
            Scalar pc = 9.81 * 1000.0 * (y - (-5E-4*x+5));
            if (pc < 0.0) pc = 0.0;

            sw = invertPcgw_(pc, this->spatialParams().fluidMatrixInteractionAtPos(globalPos));
            if (sw < swr) sw = swr;
            if (sw > 1.-sgr) sw = 1.-sgr;

            values[pressureIdx] = 1e5 ;
            values[switch1Idx] = sw;
            values[switch2Idx] = 1.e-6;
        }else {
            values[pressureIdx] = 1e5 + 9.81 * 1000.0 * ((-5E-4*x+5) - y);
            values[switch1Idx] = 1.-sgr;
            values[switch2Idx] = 1.e-6;
        }
        return values;
    }

    template<class FluidMatrixInteraction>
    static Scalar invertPcgw_(Scalar pcIn, const FluidMatrixInteraction& fluidmatrixinteraction)
    {
        Scalar lower,upper;
        int k;
        int maxIt = 50;
        Scalar bisLimit = 1.;
        Scalar sw, pcgw;
        lower=0.0; upper=1.0;
        for (k=1; k<=25; k++)
        {
            sw = 0.5*(upper+lower);
            pcgw = fluidmatrixinteraction.pcgw(sw, 0.0/*sn*/);
            Scalar delta = pcgw-pcIn;
            if (delta<0.) delta*=-1.;
            if (delta<bisLimit)
            {
                return(sw);
            }
            if (k==maxIt) {
                return(sw);
            }
            if (pcgw>pcIn) lower=sw;
            else upper=sw;
        }
        return(sw);
    }

    Scalar temperature_;
    static constexpr Scalar eps_ = 1e-6;
    std::string name_;
};
} //end namespace Dumux

#endif
