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
 * \brief Definition of a problem, for the 1p2cadsorption problem:
 * Component transport of CO2 in the CH4 phase.
 */
#ifndef DUMUX_1P2C_ADSORPTION_OUTFLOW_PROBLEM_HH
#define DUMUX_1P2C_ADSORPTION_OUTFLOW_PROBLEM_HH

#include <dumux/porousmediumflow/1p2cadsorption/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>

#include <dumux/material/fluidsystems/ch4co2.hh>
#include "1p2cadsorptionoutflowspatialparams.hh"

#define NONISOTHERMAL 0

namespace Dumux
{

template <class TypeTag>
class OnePTwoCAdsorptionOutflowProblem;

namespace Properties
{
#if NONISOTHERMAL
NEW_TYPE_TAG(OnePTwoCAdsorptionOutflowProblem, INHERITS_FROM(OnePTwoCAdsorptionNI));
NEW_TYPE_TAG(OnePTwoCAdsorptionOutflowBoxProblem, INHERITS_FROM(BoxModel, OnePTwoCAdsorptionOutflowProblem));
NEW_TYPE_TAG(OnePTwoCAdsorptionOutflowCCProblem, INHERITS_FROM(CCModel, OnePTwoCAdsorptionOutflowProblem));
#else
NEW_TYPE_TAG(OnePTwoCAdsorptionOutflowProblem, INHERITS_FROM(OnePTwoCAdsorption));
NEW_TYPE_TAG(OnePTwoCAdsorptionOutflowBoxProblem, INHERITS_FROM(BoxModel, OnePTwoCAdsorptionOutflowProblem));
NEW_TYPE_TAG(OnePTwoCAdsorptionOutflowCCProblem, INHERITS_FROM(CCModel, OnePTwoCAdsorptionOutflowProblem));
#endif

// Set the grid type
#if HAVE_UG
SET_TYPE_PROP(OnePTwoCOutflowProblem, Grid, Dune::UGGrid<2>);
#else
SET_TYPE_PROP(OnePTwoCAdsorptionOutflowProblem, Grid, Dune::YaspGrid<2>);
#endif

// Set the problem property
SET_TYPE_PROP(OnePTwoCAdsorptionOutflowProblem, Problem, OnePTwoCAdsorptionOutflowProblem<TypeTag>);

// Set fluid configuration
SET_PROP(OnePTwoCAdsorptionOutflowProblem, FluidSystem)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dumux::CO2Tables CO2Tables;
    static const bool useComplexRelations = false;
    typedef Dumux::FluidSystems::CH4CO2<TypeTag, Scalar, CO2Tables, useComplexRelations> type;
};

// Set the spatial parameters
SET_TYPE_PROP(OnePTwoCAdsorptionOutflowProblem,
              SpatialParams,
              OnePTwoCAdsorptionOutflowSpatialParams<TypeTag>);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(OnePTwoCAdsorptionOutflowProblem, UseMoles, true);

// Enable velocity output
SET_BOOL_PROP(OnePTwoCAdsorptionOutflowProblem, VtkAddVelocity, true);

// Disable gravity
SET_BOOL_PROP(OnePTwoCAdsorptionOutflowProblem, ProblemEnableGravity, false);
}


/*!
 * \ingroup OnePTwoCAdsorptionModel
 * \ingroup ImplicitTestProblems
 *
 * \brief Definition of a problem, for the 1p2c problem:
 * CH4 is present and CO2 is injected and transportet from the bottom to the top.
 *
 * The model domain is a column with
 * and homogeneous coal properties (\f$ \mathrm{K=17e-7, \Phi=0.4851, \tau=0.28}\f$).
 * Initially the domain is filled with pure CH4.
 *
 * The CO2 phase flows from the bottom to top.
 *
 * The model is able to use either mole or mass fractions. The property useMoles can be set to either true or false in the
 * problem file. Make sure that the according units are used in the problem setup. The default setting for useMoles is true.
 *
 * This problem uses the \ref OnePTwoCModel model.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box1p2cadsorption -parameterFile ./test_1p2cadsorption.input</tt> or
 * <tt>./test_cc1p2cadsorption -parameterFile ./test_1p2cadsorption.input</tt>
 */
template <class TypeTag>
class OnePTwoCAdsorptionOutflowProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        // world dimension
        dimWorld        = GridView::dimensionworld,
        nPhaseIdx       = FluidSystem::nPhaseIdx,
        nCompIdx        = FluidSystem::nCompIdx,
        numComponents   = FluidSystem::numComponents,
        CH4Idx          = FluidSystem::CH4Idx,
        TCIdx           = FluidSystem::CO2Idx,
    };
    enum {
        // indices of the primary variables
        pressureIdx         = Indices::pressureIdx,
        massOrMoleFracIdx   = Indices::massOrMoleFracIdx,
#if NONISOTHERMAL
        temperatureIdx = Indices::temperatureIdx
#endif
    };
    enum {
        // index of the transport equation
        conti0EqIdx     = Indices::conti0EqIdx,
        transportEqIdx  = Indices::transportEqIdx,
#if NONISOTHERMAL
        energyEqIdx = Indices::energyEqIdx
#endif
    };


    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    //! property that defines whether mole or mass fractions are used
        static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

public:
    OnePTwoCAdsorptionOutflowProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        //initialize fluid system
        FluidSystem::init();

        name_       = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);
        CO2inj_     = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, CO2Inj);
        pressure_   = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Initial, ReservoirPressure);
        initCO2_    = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Initial, initCO2);
        VCH4_       = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SorptionCoefficients, VCH4);
        bCH4_       = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SorptionCoefficients, bCH4);
        VCO2_       = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SorptionCoefficients, VCO2);
        bCO2_       = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SorptionCoefficients, bCO2);

        //stating in the console whether mole or mass fractions are used
        if(useMoles)
        {
            std::cout<<"problem uses mole fractions"<<std::endl;
        }
        else
        {
            std::cout<<"problem uses mass fractions"<<std::endl;
        }
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
    {
        return name_;
    }

#if !NONISOTHERMAL
    /*!
     * \brief Returns the temperature within the domain [K].
     *
     * This problem assumes a temperature of 20 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 40; } // in [K]

#endif

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
    void boundaryTypesAtPos(BoundaryTypes &values,
                            const GlobalPosition &globalPos) const
    {
        if(globalPos[1] > this->bBoxMax()[1] / 2 - eps_ && globalPos[1] < this->bBoxMax()[1] /2 + eps_)
        {
            values.setDirichlet(pressureIdx);
        }
        else
        {
            values.setAllNeumann();
        }
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The position for which the bc type should be evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        values[pressureIdx] = pressure_;
//        priVars[massOrMoleFracIdx] = initCO2_
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * This is the method for the case where the Neumann condition is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param values The neumann values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^2 \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param intersection The intersection between element and boundary
     * \param scvIdx The local subcontrolvolume index
     * \param boundaryFaceIdx The index of the boundary face
     * \param elemVolVars All volume variables for the element
     *
     * For this method, the \a values parameter stores the flux
     * in normal direction of each phase. Negative values mean influx.
     * E.g. for the mass balance that would the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    void solDependentNeumann(PrimaryVariables &values,
                      const Element &element,
                      const FVElementGeometry &fvGeometry,
                      const Intersection &intersection,
                      const int scvIdx,
                      const int boundaryFaceIdx,
                      const ElementVolumeVariables &elemVolVars) const
    {
        values = 0;
        const GlobalPosition globalPos = element.geometry().corner(scvIdx);

        //Set no flow boundary conditions for the symmetry axes of the generic models
        Scalar diameter = this->bBoxMax()[0];
        Scalar flux = CO2inj_/(3.14*diameter*diameter/4.); //[m/s] v=1 m/d = 1.1547e-5 m/s mol/(m^2*s)

        if(globalPos[1] <= eps_)
        {
            values[TCIdx] = -flux;
        }
        else if(globalPos[1] > this->bBoxMax()[1] - eps_)
        {
            Scalar moleFracCH4 = elemVolVars[scvIdx].moleFraction(CH4Idx);
            Scalar moleFracTC = elemVolVars[scvIdx].moleFraction(TCIdx);
            values[CH4Idx] = flux*moleFracCH4; //mol/(m^2*s) TODO:Check units again
            values[TCIdx]  = flux*moleFracTC; //mol/(m^2*s)
        }
        else
            values = 0.0;
    }

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * For this method, the \a priVars parameter stores the rate mass
     * of a component is generated or annihilate per volume
     * unit. Positive values mean that mass is created, negative ones
     * mean that it vanishes.
     *
     * The units must be according to either using mole or mass fractions. (mole/(m^3*s) or kg/(m^3*s))
     */
    void sourceAtPos(PrimaryVariables &priVars,
                     const GlobalPosition &globalPos) const
    {
        priVars = Scalar(0.0);
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The initial values for the primary variables
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    void initialAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
    }

    /*!
     * \brief Returns the V constant for ad- and desorption for each component in \f$\mathrm{[mol/m^3]}\f$.
     *
     * \param compIdx The component index
     */
    const Scalar V(int compIdx) const
    {
        Scalar V_[numComponents] = {};

        V_[nCompIdx] = VCH4_*1000; // *1000 for conversion from mol/l to mol/m3
        V_[TCIdx]    = VCO2_*1000;
        return V_[compIdx];
    }

    /*!
     * \brief Returns the b constant for ad- and desorption for each component in \f$\mathrm{[1/Pa]}\f$.
     *
     * \param compIdx The component index
     */
    const Scalar b(int compIdx) const
    {
        Scalar b_[numComponents] = {};
        b_[nCompIdx] = bCH4_;
        b_[TCIdx]   = bCO2_;
        return b_[compIdx];
    }
    // \}

private:
    // the internal method for the initial condition
    void initial_(PrimaryVariables &priVars,
                  const GlobalPosition &globalPos) const
    {
        priVars[pressureIdx] = pressure_; // initial condition for the pressure
        priVars[massOrMoleFracIdx] = initCO2_;  // initial condition for the CO2 molefraction
    }

    static constexpr Scalar eps_ = 1e-6;
    std::string name_;

    Scalar pressure_;
    Scalar CO2inj_;
    Scalar initCO2_;
    Scalar VCH4_;
    Scalar bCH4_;
    Scalar VCO2_;
    Scalar bCO2_;
};

} //end namespace
#endif
