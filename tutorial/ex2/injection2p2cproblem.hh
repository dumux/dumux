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
 *
 * \brief Problem where air is injected under a low permeable layer in a depth of 2700m.
 */
#ifndef DUMUX_INJECTION_2P2C_PROBLEM_HH
#define DUMUX_INJECTION_2P2C_PROBLEM_HH

#include <dumux/porousmediumflow/2p2c/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>

#include "injection2p2cspatialparams.hh"

// TODO: dumux-course-task
// Include the local residual header

namespace Dumux
{

// foward declaration
template <class TypeTag>
class Injection2p2cProblem;

// setup property TypeTag
namespace Properties
{
NEW_TYPE_TAG(Injection2p2cTypeTag, INHERITS_FROM(TwoPTwoC, InjectionSpatialParams));
NEW_TYPE_TAG(Injection2p2pcCCTypeTag, INHERITS_FROM(CCTpfaModel, Injection2p2cTypeTag));

// Set the grid type
SET_TYPE_PROP(Injection2p2cTypeTag, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(Injection2p2cTypeTag, Problem, Injection2p2cProblem<TypeTag>);

// TODO: dumux-course-task
// change the local residual type to MyTwoPTwoCLocalResidual<TypeTag>

// Set fluid configuration
SET_TYPE_PROP(Injection2p2cTypeTag,
              FluidSystem,
              FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), true /*useComplexRelations*/>);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(Injection2p2cTypeTag, UseMoles, true);

} // end namespace Properties

/*!
 * \ingroup TwoPTwoCModel
 * \ingroup ImplicitTestProblems
 * \brief Problem where air is injected under a low permeable layer in a depth of 2700m.
 *
 * The domain is sized 60m times 40m and consists of two layers, a moderately
 * permeable one for \f$ y<22m\f$ and one with a lower permeablility
 * in the rest of the domain.
 *
 * Nitrogen is injected into a water-filled aquifer through a well. First, we inject for one month.
 * Then, we continue simulating the development of the nitrogen plume for 10 years.
 * This is realized with a Neumann boundary condition at the right boundary
 * (\f$ 7m<y<15m\f$). The aquifer is situated 2700m below sea level (the depth can be changed through the input file).
 * The injected fluid phase migrates upwards due to buoyancy.
 * It accumulates and partially enters the top layer lower permeable aquitard.
 *
 * The default setting for useMoles is true, i.e. each component is balaced in units of mole.
 * The property useMoles can be set to either true or false in the
 * problem file. If you change this, make sure that the according units are used in the problem setup.
 *
 * This problem uses the \ref TwoPTwoCModel.
 */
template <class TypeTag>
class Injection2p2cProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;

    enum { dimWorld = GridView::dimensionworld };
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    Injection2p2cProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
         // initialize the tables of the fluid system
        FluidSystem::init(/*tempMin=*/273.15,
                          /*tempMax=*/423.15,
                          /*numTemp=*/50,
                          /*pMin=*/0.0,
                          /*pMax=*/30e6,
                          /*numP=*/300);

        // name of the problem and output file
        name_ = getParam<std::string>("Problem.Name");
        // depth of the aquifer, units: m
        aquiferDepth_ = getParam<Scalar>("Problem.AquiferDepth");
        // inflow rate of nitrogen water vapor mixture, units: kg/(s m^2)
        totalAreaSpecificInflow_ = getParam<Scalar>("Problem.TotalAreaSpecificInflow");
        // the duration of the injection, units: second
        injectionDuration_ = getParam<Scalar>("Problem.InjectionDuration");
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
    const std::string& name() const
    { return name_; }

    /*!
     * \brief Returns the temperature in \f$ K \f$
     */
    Scalar temperature() const
    { return 273.15 + 30; }

    // \}

     /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param bcTypes The boundary types for the conservation equations
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
         BoundaryTypes bcTypes;
        if (globalPos[dimWorld-1] < eps_)
            bcTypes.setAllDirichlet();

        // and Neuman boundary conditions everywhere else
        // note that we don't differentiate between Neumann and Robin boundary types
        else
            bcTypes.setAllNeumann();

        return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return initialAtPos(globalPos);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param globalPos The position of the integration point of the boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    PrimaryVariables neumannAtPos(const GlobalPosition &globalPos) const
    {
        // initialize values to zero, i.e. no-flow Neumann boundary conditions
        PrimaryVariables values(0.0);

        //if we are inside the injection zone set inflow Neumann boundary conditions
        if (time_ < injectionDuration_
            && globalPos[1] < 15 + eps_ && globalPos[1] > 7 - eps_ && globalPos[0] > 0.9*this->fvGridGeometry().bBoxMax()[0])
        {
            // set the Neumann values for the Nitrogen component balance
            // convert from units kg/(s*m^2) to mole/(s*m^2)
            values[Indices::conti0EqIdx + FluidSystem::H2OIdx] = -totalAreaSpecificInflow_/FluidSystem::molarMass(FluidSystem::N2Idx);
            values[Indices::conti0EqIdx + FluidSystem::N2Idx] = 0.0;
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
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values.setState(Indices::firstPhaseOnly);

        // get the water density at atmospheric conditions
        const Scalar densityW = FluidSystem::H2O::liquidDensity(temperature(), 1.0e5);

        // assume an intially hydrostatic liquid pressure profile
        // note: we subtract rho_w*g*h because g is defined negative
        const Scalar pw = 1.0e5 - densityW*this->gravity()[dimWorld-1]*(aquiferDepth_ - globalPos[dimWorld-1]);

        // initially we have some nitrogen dissolved
        // saturation mole fraction would be
        // moleFracLiquidN2 = (pw + pc + p_vap^sat)/henry;
        const Scalar moleFracLiquidN2 = pw*0.95/BinaryCoeff::H2O_N2::henry(temperature());

        // note that because we start with a single phase system the primary variables
        // are pl and x^w_N2. This will switch as soon after we start injecting to a two
        // phase system so the primary variables will be pl and Sn (non-wetting saturation).
        values[Indices::switchIdx] = moleFracLiquidN2;
        values[Indices::pressureIdx] = pw;

        return values;
    }

    // \}

    void setTime(Scalar time)
    {
        time_ = time;
    }

private:
    static constexpr Scalar eps_ = 1e-6;
    std::string name_; //! Problem name
    Scalar aquiferDepth_; //! Depth of the aquifer in m
    Scalar totalAreaSpecificInflow_; //! Area specific inflow rate in mole/(s*m^2)
    Scalar time_;
    Scalar injectionDuration_; //! Duration of the injection in seconds
};

} // end namespace Dumux

#endif
