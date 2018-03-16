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
 * \brief Problem where air is injected under a low permeable layer in a depth of 2700m.
 */
#ifndef DUMUX_INJECTION_PROBLEM_HH
#define DUMUX_INJECTION_PROBLEM_HH

#include <dumux/discretization/cellcentered/mpfa/properties.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/box/properties.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/2p2c/model.hh>
#include <dumux/material/fluidsystems/h2on2.hh>

#include "injectionspatialparams.hh"

namespace Dumux
{

#ifndef ENABLECACHING
#define ENABLECACHING 0
#endif

template <class TypeTag>
class InjectionProblem;

namespace Properties
{
NEW_TYPE_TAG(InjectionTypeTag, INHERITS_FROM(TwoPTwoC, InjectionSpatialParams));
NEW_TYPE_TAG(InjectionBoxTypeTag, INHERITS_FROM(BoxModel, InjectionTypeTag));
NEW_TYPE_TAG(InjectionCCTpfaTypeTag, INHERITS_FROM(CCTpfaModel, InjectionTypeTag));
NEW_TYPE_TAG(InjectionCCMpfaTypeTag, INHERITS_FROM(CCMpfaModel, InjectionTypeTag));

// Set the grid type
SET_TYPE_PROP(InjectionTypeTag, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(InjectionTypeTag, Problem, InjectionProblem<TypeTag>);

// Set fluid configuration
SET_TYPE_PROP(InjectionTypeTag,
              FluidSystem,
              FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), false /*useComplexRelations*/>);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(InjectionTypeTag, UseMoles, true);

// Enable caching or not (reference solutions created without caching)
SET_BOOL_PROP(InjectionTypeTag, EnableFVGridGeometryCache, ENABLECACHING);
SET_BOOL_PROP(InjectionTypeTag, EnableGridVolumeVariablesCache, ENABLECACHING);
SET_BOOL_PROP(InjectionTypeTag, EnableGridFluxVariablesCache, ENABLECACHING);
}


/*!
 * \ingroup TwoPTwoCTests
 * \brief Problem where air is injected under a low permeable layer in a depth of 2700m.
 *
 * The domain is sized 60m times 40m and consists of two layers, a moderately
 * permeable one (\f$ K=10e-12\f$) for \f$ y<22m\f$ and one with a lower permeablility (\f$ K=10e-13\f$)
 * in the rest of the domain.
 *
 * A mixture of Nitrogen and Water vapor, which is composed according to the prevailing conditions (temperature, pressure)
 * enters a water-filled aquifer. This is realized with a solution-dependent Neumann boundary condition at the right boundary
 * (\f$ 5m<y<15m\f$). The aquifer is situated 2700m below sea level. The injected fluid phase migrates upwards due to buoyancy.
 * It accumulates and partially enters the lower permeable aquitard.
 *
 * The model is able to use either mole or mass fractions. The property useMoles can be set to either true or false in the
 * problem file. Make sure that the according units are used in the problem setup. The default setting for useMoles is true.
 *
 * This problem uses the \ref TwoPTwoCModel.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box2p2c</tt> or
 * <tt>./test_cc2p2c</tt>
 */
template <class TypeTag>
class InjectionProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum {
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,

        // world dimension
        dimWorld = GridView::dimensionworld
    };

    enum {
        nPhaseIdx = Indices::nPhaseIdx,

        wPhaseOnly = Indices::wPhaseOnly,

        wCompIdx = FluidSystem::wCompIdx,
        nCompIdx = FluidSystem::nCompIdx,

        contiH2OEqIdx = Indices::contiWEqIdx,
        contiN2EqIdx = Indices::contiNEqIdx
    };

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    //! property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    InjectionProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        nTemperature_       = getParam<int>("Problem.NTemperature");
        nPressure_          = getParam<int>("Problem.NPressure");
        pressureLow_        = getParam<Scalar>("Problem.PressureLow");
        pressureHigh_       = getParam<Scalar>("Problem.PressureHigh");
        temperatureLow_     = getParam<Scalar>("Problem.TemperatureLow");
        temperatureHigh_    = getParam<Scalar>("Problem.TemperatureHigh");
        temperature_        = getParam<Scalar>("Problem.InitialTemperature");
        depthBOR_           = getParam<Scalar>("Problem.DepthBOR");
        name_               = getParam<std::string>("Problem.Name");

        // initialize the tables of the fluid system
        FluidSystem::init(/*Tmin=*/temperatureLow_,
                          /*Tmax=*/temperatureHigh_,
                          /*nT=*/nTemperature_,
                          /*pmin=*/pressureLow_,
                          /*pmax=*/pressureHigh_,
                          /*np=*/nPressure_);

        //stateing in the console whether mole or mass fractions are used
        if(useMoles)
            std::cout<<"problem uses mole-fractions"<<std::endl;
        else
            std::cout<<"problem uses mass-fractions"<<std::endl;
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
     * \brief Returns the temperature [K]
     */
    Scalar temperature() const
    { return temperature_; }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param values Stores the value of the boundary type
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        if (globalPos[0] < eps_)
            bcTypes.setAllDirichlet();
        else
            bcTypes.setAllNeumann();
        return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param values Stores the Dirichlet values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} ] \f$
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return initial_(globalPos);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann
     *        boundary segment in dependency on the current solution.
     *
     * \param values Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param intersection The intersection between element and boundary
     * \param scvIdx The local index of the sub-control volume
     * \param boundaryFaceIdx The index of the boundary face
     * \param elemVolVars All volume variables for the element
     *
     * This method is used for cases, when the Neumann condition depends on the
     * solution and requires some quantities that are specific to the fully-implicit method.
     * The \a values store the mass flux of each phase normal to the boundary.
     * Negative values indicate an inflow.
     */
    NumEqVector neumann(const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const ElementVolumeVariables& elemVolVars,
                          const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        const auto& globalPos = scvf.ipGlobal();

        Scalar injectedPhaseMass = 1e-3;
        Scalar moleFracW = elemVolVars[scvf.insideScvIdx()].moleFraction(nPhaseIdx, wCompIdx);
        if (globalPos[1] < 14 - eps_ && globalPos[1] > 6.5 - eps_)
        {
            values[contiN2EqIdx] = -(1-moleFracW)*injectedPhaseMass/FluidSystem::molarMass(nCompIdx); //mole/(m^2*s) -> kg/(s*m^2)
            values[contiH2OEqIdx] = -moleFracW*injectedPhaseMass/FluidSystem::molarMass(wCompIdx); //mole/(m^2*s) -> kg/(s*m^2)
        }
        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * \param values Stores the initial values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        return initial_(globalPos);
    }

    // \}

private:
    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * The internal method for the initial condition
     *
     * \param values Stores the initial values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars.setState(wPhaseOnly);

        Scalar densityW = FluidSystem::H2O::liquidDensity(temperature_, 1e5);

        Scalar pl = 1e5 - densityW*this->gravity()[1]*(depthBOR_ - globalPos[1]);
        Scalar moleFracLiquidN2 = pl*0.95/BinaryCoeff::H2O_N2::henry(temperature_);
        Scalar moleFracLiquidH2O = 1.0 - moleFracLiquidN2;

        Scalar meanM =
            FluidSystem::molarMass(wCompIdx)*moleFracLiquidH2O +
            FluidSystem::molarMass(nCompIdx)*moleFracLiquidN2;
        if(useMoles)
        {
            //mole-fraction formulation
            priVars[switchIdx] = moleFracLiquidN2;
        }
        else
        {
            //mass fraction formulation
            Scalar massFracLiquidN2 = moleFracLiquidN2*FluidSystem::molarMass(nCompIdx)/meanM;
            priVars[switchIdx] = massFracLiquidN2;
        }
        priVars[pressureIdx] = pl;
        return priVars;
    }

    Scalar temperature_;
    Scalar depthBOR_;
    static constexpr Scalar eps_ = 1e-6;

    int nTemperature_;
    int nPressure_;
    Scalar pressureLow_, pressureHigh_;
    Scalar temperatureLow_, temperatureHigh_;

    std::string name_;
};

} // end namespace Dumux

#endif
