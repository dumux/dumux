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

#include <dune/grid/yaspgrid.hh>

#include <dumux/porousmediumflow/1pncmin/model.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/box/properties.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/cellcentered/mpfa/properties.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/fluidsystems/1padapter.hh>
#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>
#include <dumux/material/components/cao2h2.hh>
#include <dumux/material/solidsystems/compositionalsolidphase.hh>

#include "spatialparams.hh"
#include "reaction.hh"
#include "modifiedcao.hh"

namespace Dumux {

template <class TypeTag>
class ThermoChemProblem;

namespace Properties {
NEW_TYPE_TAG(ThermoChemTypeTag, INHERITS_FROM(OnePNCMinNI));
NEW_TYPE_TAG(ThermoChemBoxTypeTag, INHERITS_FROM(BoxModel, ThermoChemTypeTag));

// Set the grid type
SET_TYPE_PROP(ThermoChemTypeTag, Grid, Dune::YaspGrid<2>);
// Set the problem property
SET_TYPE_PROP(ThermoChemTypeTag, Problem, ThermoChemProblem<TypeTag>);

// The fluid system
SET_PROP(ThermoChemTypeTag, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using H2ON2 = FluidSystems::H2ON2<Scalar>;
    static constexpr auto phaseIdx = H2ON2::gasPhaseIdx; // simulate the air phase
    using type = FluidSystems::OnePAdapter<H2ON2, phaseIdx>;
};

SET_PROP(ThermoChemTypeTag, SolidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ComponentOne = Components::ModifiedCaO<Scalar>;
    using ComponentTwo = Components::CaO2H2<Scalar>;
    using type = SolidSystems::CompositionalSolidPhase<Scalar, ComponentOne, ComponentTwo>;
};

// // Enable velocity output
// SET_BOOL_PROP(ThermoChemTypeTag, VtkAddVelocity, false);

// Set the spatial parameters
SET_PROP(ThermoChemTypeTag, SpatialParams)
{
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = ThermoChemSpatialParams<FVGridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(ThermoChemTypeTag, UseMoles, true);
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
    using SolidSystem = typename GET_PROP_TYPE(TypeTag, SolidSystem);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ReactionRate = ThermoChemReaction;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    enum
    {
        // Indices of the primary variables
        pressureIdx = Indices::pressureIdx, //gas-phase pressure
        H2OIdx = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::H2OIdx), // mole fraction water

        CaOIdx = FluidSystem::numComponents,
        CaO2H2Idx = FluidSystem::numComponents+1,

        // Equation Indices
        conti0EqIdx = Indices::conti0EqIdx,

        // Phase Indices
        cPhaseIdx = SolidSystem::comp0Idx,

        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx
    };

    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;

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
        FluidSystem::init(/*tempMin=*/473.15,
                          /*tempMax=*/623.0,
                          /*numTemptempSteps=*/25,
                          /*startPressure=*/0,
                          /*endPressure=*/9e6,
                          /*pressureSteps=*/200);

        // obtain BCs
        boundaryPressure_ = getParam<Scalar>("Problem.BoundaryPressure");
        boundaryVaporMoleFrac_ = getParam<Scalar>("Problem.BoundaryMoleFraction");
        boundaryTemperature_ = getParam<Scalar>("Problem.BoundaryTemperature");

        unsigned int codim = GET_PROP_TYPE(TypeTag, FVGridGeometry)::discMethod == DiscretizationMethod::box ? dim : 0;
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
        values.setDirichlet(H2OIdx);
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
        priVars[H2OIdx] = boundaryVaporMoleFrac_;
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

    NumEqVector neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const SubControlVolumeFace& scvf) const
    {
        NumEqVector flux(0.0);
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
        priVars[H2OIdx]   = h2oInit;
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
    NumEqVector source(const Element &element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume &scv) const
    {

        NumEqVector source(0.0);
        const auto& volVars = elemVolVars[scv];

        Scalar qMass = rrate_.thermoChemReaction(volVars);
        const auto elemSol = elementSolution(element, elemVolVars, fvGeometry);
        Scalar qMole = qMass/FluidSystem::molarMass(H2OIdx)*(1-volVars.porosity());

        // make sure not more solid reacts than present
        // In this test, we only consider discharge. Therefore, we use the cPhaseIdx for CaO.
        if (-qMole*timeStepSize_ + volVars.solidVolumeFraction(cPhaseIdx)* volVars.solidComponentMolarDensity(cPhaseIdx) < 0 + eps_)
        {
            qMole = -volVars.solidVolumeFraction(cPhaseIdx)* volVars.solidComponentMolarDensity(cPhaseIdx)/timeStepSize_;
        }
        source[conti0EqIdx+CaO2H2Idx] = qMole;
        source[conti0EqIdx+CaOIdx] = - qMole;
        source[conti0EqIdx+H2OIdx] = - qMole;

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
            const auto elemSol = elementSolution(element, curSol, this->fvGridGeometry());

            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                VolumeVariables volVars;
                volVars.update(elemSol, *this, element, scv);
                const auto dofIdxGlobal = scv.dofIndex();
                permeability_[dofIdxGlobal] = this->spatialParams().permeability(element, scv, elemSol);
                porosity_[dofIdxGlobal] = volVars.porosity();
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
