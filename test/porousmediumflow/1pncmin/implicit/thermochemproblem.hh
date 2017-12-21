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
 * \ingroup OnePNCMinTests
 * \brief Definition of a problem for thermochemical heat storage using \f$ \textnormal{CaO},   \textnormal{Ca} \left( \textnormal{OH} \right)_2\f$.
 */
#ifndef DUMUX_THERMOCHEM_PROBLEM_HH
#define DUMUX_THERMOCHEM_PROBLEM_HH

#include <dumux/porousmediumflow/1pncmin/model.hh>
#include <dumux/discretization/box/properties.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/cellcentered/mpfa/properties.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>

#include "thermochemspatialparams.hh"
#include "thermochemreaction.hh"
#include "modifiedsteamn2cao2h2.hh"


namespace Dumux
{

template <class TypeTag>
class ThermoChemProblem;

namespace Properties
{
NEW_TYPE_TAG(ThermoChemProblem, INHERITS_FROM(OnePNCMinNI, ThermoChemSpatialParams));
NEW_TYPE_TAG(ThermoChemBoxProblem, INHERITS_FROM(BoxModel, ThermoChemProblem));

// Set the grid type
SET_TYPE_PROP(ThermoChemProblem, Grid, Dune::YaspGrid<2>);
// Set the problem property
SET_TYPE_PROP(ThermoChemProblem, Problem, ThermoChemProblem<TypeTag>);
// Set fluid configuration
SET_PROP(ThermoChemProblem, FluidSystem)
{ /*private:*/
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::ModifiedSteamN2CaO2H2<Scalar>;
};

// // Enable velocity output
// SET_BOOL_PROP(ThermoChemProblem, VtkAddVelocity, false);

// Set the spatial parameters
SET_TYPE_PROP(ThermoChemProblem, SpatialParams, ThermoChemSpatialParams<TypeTag>);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(ThermoChemProblem, UseMoles, true);
}

/*!
 * \ingroup OnePNCMinTests
 *
 * \brief Test for the 1pncmin model in combination with the NI model for a quasi batch
 * reaction of Calciumoxyde to Calciumhydroxide.
 *
 * The boundary conditions of the batch test are such, that there are no gradients for temperature, pressure and gas water concentration within the reactor.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_1pncminni_box -ParameterFile </tt>
 * The test only runs for the box discretization.
 */
template <class TypeTag>
class ThermoChemProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);;
    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using ResidualVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using ReactionRate =ThermoChemReaction<TypeTag>;


    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    enum
    {
        numComponents = FluidSystem::numComponents,

        // Indices of the primary variables
        pressureIdx = Indices::pressureIdx, //gas-phase pressure
        firstMoleFracIdx = Indices::firstMoleFracIdx, // mole fraction water

        CaOIdx = FluidSystem::numComponents,
        CaO2H2Idx = FluidSystem::numComponents+1,

        //Equation Indices
        conti0EqIdx = Indices::conti0EqIdx,
        firstTransportEqIdx = Indices::firstTransportEqIdx,

        // Phase Indices
        phaseIdx = FluidSystem::gPhaseIdx,
        cPhaseIdx = FluidSystem::cPhaseIdx,
        hPhaseIdx = FluidSystem::hPhaseIdx,

        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx
    };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;


public:
    /*!
     * \brief The constructor
     *
     * \param fvGridGeometry The finite volume grid geometry
     */
    ThermoChemProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
        : ParentType(fvGridGeometry)
    {
        name_      = getParam<std::string>("Problem.Name");

        // obtain BCs
        boundaryPressure_ = getParam<Scalar>("Problem.BoundaryPressure");
        boundaryVaporMoleFrac_ = getParam<Scalar>("Problem.BoundaryMoleFraction");
        boundaryTemperature_ = getParam<Scalar>("Problem.BoundaryTemperature");

        unsigned int codim = GET_PROP_VALUE(TypeTag, DiscretizationMethod) == DiscretizationMethods::Box ? dim : 0;
        permeability_.resize(fvGridGeometry->gridView().size(codim));
        porosity_.resize(fvGridGeometry->gridView().size(codim));
        reactionRate_.resize(fvGridGeometry->gridView().size(codim));
    }

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string name() const
    { return name_; }

    /*!
     * \brief Sets the currently used time step size.
     *
     * This is necessary to limit the source terms to the maximum possible rate.
     */
    void setTimeStepSize( Scalar timeStepSize )
     {
        timeStepSize_ = timeStepSize;
     }

    /*!
     * \name Boundary conditions
     *
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos( const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        // we don't set any BCs for the solid phases
        values.setDirichlet(pressureIdx);
        values.setDirichlet(firstMoleFracIdx);
        values.setDirichlet(temperatureIdx);

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);

        priVars[pressureIdx] = boundaryPressure_;
        priVars[firstMoleFracIdx] = boundaryVaporMoleFrac_;
        priVars[temperatureIdx] = boundaryTemperature_;
        priVars[CaO2H2Idx] = 0.0;
        priVars[CaOIdx] = 0.2;

        return priVars;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann
     *        boundary segment in dependency on the current solution.
     *
     * \param element The element
     * \param fvGeometry The finite volume geometry
     * \param elemVolVars The element volume variables
     * \param scvf The subcontrolvolume face
     *
     * \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * Negative values indicate an inflow.
     */

    ResidualVector neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const SubControlVolumeFace& scvf) const
    {
        ResidualVector flux(0.0);
        return flux;
    }

    /*!
     * \brief Evaluates the initial values for a control volume in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables priVars(0.0);

        Scalar pInit;
        Scalar tInit;
        Scalar h2oInit;
        Scalar CaOInit;
        Scalar CaO2H2Init;

        pInit = getParam<Scalar>("Problem.PressureInitial");
        tInit = getParam<Scalar>("Problem.TemperatureInitial");
        h2oInit = getParam<Scalar>("Problem.VaporInitial");
        CaOInit = getParam<Scalar>("Problem.CaOInitial");
        CaO2H2Init = getParam<Scalar>("Problem.CaO2H2Initial");

        priVars[pressureIdx] = pInit;
        priVars[firstMoleFracIdx]   = h2oInit;
        priVars[temperatureIdx] = tInit;

        // these values are not used, as we didn't set BCs
        // for the solid phases. For cell-centered models it is
        // important to set the values to fully define Dirichlet BCs
        priVars[CaOIdx] = CaOInit;
        priVars[CaO2H2Idx]   = CaO2H2Init;

        return priVars;
    }

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume in units of \f$ [ \textnormal{unit of conserved quantity} / (m^3 \cdot s )] \f$.
     *
     * This is the method for the case where the source term is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The subcontrolvolume
     *
     * For this method, the \a values parameter stores the conserved quantity rate
     * generated or annihilate per volume unit. Positive values mean
     * that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
     */
    PrimaryVariables source(const Element &element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume &scv) const
    {

        PrimaryVariables source(0.0);
        const auto& volVars = elemVolVars[scv];

        Scalar qMass = rrate_.thermoChemReaction(volVars);

        ElementSolutionVector elemSol(element, elemVolVars, fvGeometry);
        Scalar qMole = qMass/FluidSystem::molarMass(firstMoleFracIdx)*(1-this->spatialParams().porosity(element, scv, elemSol));

        // make sure not more solid reacts than present
        // In this test, we only consider discharge. Therefore, we use the cPhaseIdx for CaO.
        if (-qMole*timeStepSize_ + volVars.precipitateVolumeFraction(cPhaseIdx)* volVars.molarDensity(cPhaseIdx) < 0 + eps_)
        {
            qMole = -volVars.precipitateVolumeFraction(cPhaseIdx)* volVars.molarDensity(cPhaseIdx)/timeStepSize_;
        }

        source[conti0EqIdx+CaO2H2Idx] = qMole;
        source[conti0EqIdx+CaOIdx] = - qMole;
        source[conti0EqIdx+firstMoleFracIdx] = - qMole;

        Scalar deltaH = 108e3; // J/mol
        source[energyEqIdx] = qMole * deltaH;

        return source;
    }


   /*!
     * \brief Return the permeability
     */
    const std::vector<Scalar>& getPerm()
    {
        return permeability_;
    }

   /*!
     * \brief Return the porosity
     */
    const std::vector<Scalar>& getPoro()
    {
        return porosity_;
    }

     /*!
     * \brief Return the reaction rate
     */
    const std::vector<Scalar>& getRRate()
    {
        return reactionRate_;
    }

    /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    void updateVtkOutput(const SolutionVector& curSol)
    {
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            ElementSolutionVector elemSol(element, curSol, this->fvGridGeometry());

            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                VolumeVariables volVars;
                volVars.update(elemSol, *this, element, scv);
                const auto dofIdxGlobal = scv.dofIndex();
                permeability_[dofIdxGlobal] = this->spatialParams().permeability(element, scv, elemSol);
                porosity_[dofIdxGlobal] = this->spatialParams().porosity(element, scv, elemSol);
                PrimaryVariables reactionRate;
                reactionRate_[dofIdxGlobal] = rrate_.thermoChemReaction(volVars);
            }
        }
    }

private:
    std::string name_;

    static constexpr Scalar eps_ = 1e-6;

    // boundary conditions
    Scalar boundaryPressure_;
    Scalar boundaryVaporMoleFrac_;
    Scalar boundaryTemperature_;

    std::vector<double> permeability_;
    std::vector<double> porosity_;
    std::vector<double> reactionRate_;

    ReactionRate rrate_;
    Scalar timeStepSize_;
};

} //end namespace

#endif
