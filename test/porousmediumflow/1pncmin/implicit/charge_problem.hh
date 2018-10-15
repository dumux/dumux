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
#ifndef DUMUX_CHARGE_PROBLEM_HH
#define DUMUX_CHARGE_PROBLEM_HH

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
//changes made
#include <dumux/material/fluidsystems/1padapter.hh>
// box solution dependent Neumann, outflow
#include <dumux/discretization/evalgradients.hh>
#include <dumux/discretization/evalsolution.hh>
// #include "steamn2.hh"
#include <dumux/material/fluidsystems/h2on2.hh>



#include "thermochemspatialparams.hh"
#include "thermochemreaction.hh"
#include "modifiedcao.hh"

namespace Dumux {

template <class TypeTag>
class ChargeProblem;

namespace Properties {
NEW_TYPE_TAG(ChargeTypeTag, INHERITS_FROM(OnePNCMinNI));
NEW_TYPE_TAG(ChargeBoxTypeTag, INHERITS_FROM(BoxModel, ChargeTypeTag));

// Set the grid type
SET_TYPE_PROP(ChargeTypeTag, Grid, Dune::YaspGrid<2>); // 2 dimensional grid
// Set the problem property
SET_TYPE_PROP(ChargeTypeTag, Problem, ChargeProblem<TypeTag>);

// The fluid system
SET_PROP(ChargeTypeTag, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using H2ON2 = FluidSystems::H2ON2<Scalar>;
    static constexpr auto phaseIdx = H2ON2::gasPhaseIdx; // simulate the air phase
    using type = FluidSystems::OnePAdapter<H2ON2, phaseIdx>;
};

SET_PROP(ChargeTypeTag, SolidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ComponentOne = Components::ModifiedCaO<Scalar>;
    using ComponentTwo = Components::CaO2H2<Scalar>;
    using type = SolidSystems::CompositionalSolidPhase<Scalar, ComponentOne, ComponentTwo>;
};

// // Enable velocity output
// SET_BOOL_PROP(ChargeTypeTag, VtkAddVelocity, false);

// Set the spatial parameters
SET_TYPE_PROP(ChargeTypeTag, SpatialParams, ThermoChemSpatialParams<TypeTag>);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(ChargeTypeTag, UseMoles, true);
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
class ChargeProblem : public PorousMediumFlowProblem<TypeTag>
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

    //change made
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    enum
    {
        // Indices of the primary variables
        pressureIdx = Indices::pressureIdx, //gas-phase pressure


        //change made
        firstMoleFracIdx = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::H2OIdx), // mole fraction water
        //H2OIdx = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::H2OIdx), // mole fraction water


        CaOIdx = FluidSystem::numComponents,
        CaO2H2Idx = FluidSystem::numComponents+1,

        // Equation Indices
        conti0EqIdx = Indices::conti0EqIdx,

        // Phase Indices
        phaseIdx = FluidSystem::phase0Idx,
        cPhaseIdx = SolidSystem::comp0Idx,
        //change made
        hPhaseIdx = SolidSystem::comp0Idx+1,

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
    ChargeProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
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

        //change made
        // create and initialize file for flux and storage calculations
        outputFile_.open("balance.out", std::ios::out);
        outputFile_ << "inFluxH20|outFluxH20|inFluxN2|outFluxN2|inFluxEnthaply|outFluxEnthaply" << std::endl;

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
        /*values.setDirichlet(pressureIdx);
        values.setDirichlet(H2OIdx);
        values.setDirichlet(temperatureIdx);*/

        //changes made
 //       if (globalPos[0] < eps_)
//         {
// //             values.setNeumann(pressureIdx);
// //             values.setNeumann(firstMoleFracIdx);
// //             values.setNeumann(temperatureIdx);
// //             values.setNeumann(CaO2H2Idx);
// //             values.setNeumann(CaOIdx);
//             values.setDirichlet(pressureIdx);
//             values.setDirichlet(firstMoleFracIdx);
//             values.setDirichlet(temperatureIdx);
//             values.setDirichlet(CaO2H2Idx);
//             values.setDirichlet(CaOIdx);
//         }
//
//         if( globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_ )
//         {
//             values.setDirichlet(pressureIdx);
//             values.setDirichlet(firstMoleFracIdx);
//             values.setDirichlet(temperatureIdx);
            values.setNeumann(pressureIdx);
            values.setNeumann(firstMoleFracIdx);
            values.setNeumann(temperatureIdx);
            values.setNeumann(CaO2H2Idx);
            values.setNeumann(CaOIdx);
//         }


        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param globalPos The global position
     */

    //changes made
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);

        if (globalPos[0] < eps_ )
        {
            priVars[pressureIdx] = 2e5;
            priVars[firstMoleFracIdx] = 0.464;
            priVars[temperatureIdx] = 573.15;
            priVars[CaO2H2Idx] = 0.0;
            priVars[CaOIdx] = 0.2;
        }

    /*PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);

        priVars[pressureIdx] = boundaryPressure_;
        priVars[H2OIdx] = boundaryVaporMoleFrac_;
        priVars[temperatureIdx] = boundaryTemperature_;
        priVars[CaO2H2Idx] = 0.0;
        priVars[CaOIdx] = 0.2;*/

        return priVars;
    }

    /*
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

    //changes made
    /*NumEqVector neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const SubControlVolumeFace& scvf) const
    {
        NumEqVector flux(0.0);
        return flux;
    }
    */
    NumEqVector neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const SubControlVolumeFace& scvf) const
    {
        NumEqVector flux(0.0);

        const auto& globalPos = scvf.ipGlobal();
        const auto& scv = fvGeometry.scv(scvf.insideScvIdx());

        if(globalPos[0] < eps_)
        {
           Scalar InFlowAir = 4.64;
           Scalar InFlowH2O = 0.1;
           Scalar tIn = 773.15;

           FluidState fluidstateBorder;

           fluidstateBorder.setTemperature(tIn);
           fluidstateBorder.setPressure(phaseIdx, elemVolVars[scv].pressure(phaseIdx));

           Scalar T = elemVolVars[scv].temperature();

           Scalar deltaH = 0.0; // if temperature at the right border > 573.15 K, cool down
                                //temperature of injected fluid via the enthalpyflux deltaH

           deltaH = (T-tIn)* 1e20;

           Scalar hInAir = InFlowAir*FluidSystem::molarMass(firstMoleFracIdx-1)
                           *FluidSystem::componentEnthalpy(fluidstateBorder, phaseIdx, firstMoleFracIdx-1);

           Scalar hInH2O = InFlowH2O*FluidSystem::molarMass(firstMoleFracIdx)
                           *FluidSystem::componentEnthalpy(fluidstateBorder, phaseIdx, firstMoleFracIdx);
           flux[pressureIdx] = - InFlowAir; //[mol/s] gas inflow of the air component
           flux[firstMoleFracIdx] = - InFlowH2O;//[mol/s] gas inflow of the water component
           flux[temperatureIdx] = - (hInAir + hInH2O -deltaH); //[J/s] enthalpy inflow
        }

        // outflow BC
        if(globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_ )
        {
            // set a fixed pressure on the right side of the domain
            const Scalar dirichletPressure = 1.0e5;
            const Scalar K = elemVolVars[scv].permeability();
            const Scalar molarDensity = elemVolVars[scv].molarDensity(phaseIdx);
            const Scalar density = elemVolVars[scv].density(phaseIdx);

/// BOX //////////////////////
            // construct the element solution
            const auto elemSol = [&]()
            {
                auto sol = elementSolution(element, elemVolVars, fvGeometry);

                for(auto&& scvf : scvfs(fvGeometry))
                    if(scvf.center()[0] > this->fvGridGeometry().bBoxMax()[0] - eps_)
                        sol[fvGeometry.scv(scvf.insideScvIdx()).localDofIndex()][pressureIdx] = dirichletPressure;

                return sol;
            }();
            // evaluate the gradient
            const auto gradient = [&]()->GlobalPosition
            {
                const auto grads = evalGradients(element, element.geometry(), fvGeometry.fvGridGeometry(), elemSol, globalPos);
                return grads[pressureIdx];

            }();

            // calculate the flux
            Scalar tpfaFlux = gradient * scvf.unitOuterNormal();
            tpfaFlux *= -1.0  * K;

            if(tpfaFlux < 0) tpfaFlux = 0.0; //make sure that there is no influx from the right

            Scalar tpfaFluxMole = tpfaFlux * molarDensity * elemVolVars[scv].mobility(phaseIdx);
            Scalar tpfaFluxMass = tpfaFlux * density * elemVolVars[scv].mobility(phaseIdx);

            // emulate an outflow condition for the component transport on the right side
            flux[pressureIdx] = tpfaFluxMole * (elemVolVars[scv].moleFraction(phaseIdx, (firstMoleFracIdx-1)));
            flux[firstMoleFracIdx] = tpfaFluxMole * elemVolVars[scv].moleFraction(phaseIdx, firstMoleFracIdx);
            flux[temperatureIdx] = tpfaFluxMass * (FluidSystem::enthalpy(elemVolVars[scv].fluidState(), phaseIdx));
        }
/////////////////////////////////////////////////////////////////////////////////////////////

        return flux;
    }





    /*
     *
     *
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
        //change made
        priVars[firstMoleFracIdx]   = h2oInit;
        //priVars[H2OIdx]   = h2oInit;
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


    //change made

     PrimaryVariables source(const Element &element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume &scv) const
    {

        PrimaryVariables source(0.0);
        const auto& volVars = elemVolVars[scv];

        Scalar qMass = rrate_.thermoChemReactionSimple(volVars);

        const auto elemSol = elementSolution(element, elemVolVars, fvGeometry);
        Scalar qMole = qMass/FluidSystem::molarMass(firstMoleFracIdx)*(1-volVars.porosity());

        // make sure not more solid reacts than present
        // In this test, we only consider discharge. Therefore, we use the cPhaseIdx for CaO.
        if (-qMole*timeStepSize_ + volVars.solidVolumeFraction(cPhaseIdx)* volVars.solidComponentMolarDensity(cPhaseIdx) < 0 + eps_)
        {
            qMole = -volVars.solidVolumeFraction(cPhaseIdx)* volVars.solidComponentMolarDensity(cPhaseIdx)/timeStepSize_;
        }

        source[conti0EqIdx+CaO2H2Idx] = qMole;
        source[conti0EqIdx+CaOIdx] = - qMole;
        source[conti0EqIdx+firstMoleFracIdx] = - qMole;

        Scalar deltaH = 108.3e3; // J/mol
        source[energyEqIdx] = qMole * (deltaH - 4*(volVars.pressure(phaseIdx)/volVars.molarDensity(phaseIdx))) ;

        return source;
    }



    /*NumEqVector source(const Element &element,
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
        source[energyEqIdx] = qMole * (deltaH - 4*(volVars.pressure(phaseIdx)/volVars.molarDensity(phaseIdx)));

        return source;
    }
   */

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




    //change made
    void postTimeStep(NumEqVector& neumannInFlux, NumEqVector& neumannOutFlux, const SolutionVector& curSol)
    {
       Scalar inNeumannH20 = 0.0;
        Scalar inNeumannEnthalpy = 0.0;
        Scalar outNeumannH2O = 0.0;
        Scalar outNeumannEnthaply = 0.0;

        // Loop over all elements
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            const auto elemSol = elementSolution(element, curSol, this->fvGridGeometry());
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            auto eIdx = this->fvGridGeometry().elementMapper().index(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                // loop over all subcontrolvolumefaces
                for (auto&& scvf : scvfs(fvGeometry))
                {
                    const auto idx = scvf.index();
                    const auto& ipGlobal = scvf.ipGlobal();
                    VolumeVariables volVars;
                    volVars.update(elemSol, *this, element, scv);
//                     std::cout<<"eIdx = "<<eIdx <<" idx = "<< idx << "\n";

//                     neumannFlux = this->neumann(element, fvGeometry, elemVolVars, scvf);

//                  //   advective fluxes for h2o and enthalpy
//                     if(ipGlobal[0] < /*2e-4*/0.0 + eps_)
//                     {
// //                         std::cout<<"test \n";
//                         inFluxH20 += volumeFlux[eIdx/*idx*/]
//                                     *volVars.molarDensity(phaseIdx)
//                                     *volVars.moleFraction(phaseIdx, firstMoleFracIdx);
//                         inFluxN2 += volumeFlux[eIdx/*idx*/]
//                                     *volVars.molarDensity(phaseIdx)
//                                     *volVars.moleFraction(phaseIdx, firstMoleFracIdx-1);
//                         inFluxEnthaply += volumeFlux[eIdx/*idx*/]
//                                         *volVars.density(phaseIdx)
//                                         *(FluidSystem::enthalpy(volVars.fluidState(), phaseIdx));
//                     }

                    if(ipGlobal[0] < 0.0 + eps_)
                    {
                        inNeumannH20 = neumannInFlux[firstMoleFracIdx];
                        inNeumannEnthalpy = neumannInFlux[temperatureIdx];
                    }

//                     if(ipGlobal[0] > 0.08 - /*10*/*eps_)
//                     {
//                         outFluxH20 += volumeFlux[eIdx/*idx*/]
//                                     *volVars.molarDensity(phaseIdx)
//                                     *volVars.moleFraction(phaseIdx, firstMoleFracIdx);
//                         outFluxN2 += volumeFlux[eIdx/*idx*/]
//                                     *volVars.molarDensity(phaseIdx)
//                                     *volVars.moleFraction(phaseIdx, firstMoleFracIdx-1);
//                         outFluxEnthaply += volumeFlux[eIdx/*idx*/]
//                                         *volVars.density(phaseIdx)
//                                         *(FluidSystem::enthalpy(volVars.fluidState(), phaseIdx));
//                     }

                    if(ipGlobal[0] > 0.08 - eps_)
                    {
                        outNeumannH2O = neumannOutFlux[firstMoleFracIdx];
                        outNeumannEnthaply = neumannOutFlux[temperatureIdx];
                    }
                }
            }
        }

        outputFile_ << /*inFluxH20 <<" | "<< outFluxH20 <<" | "<<*/ inNeumannH20 <<" | "<< outNeumannH2O <<" | "<< /*inFluxN2 <<" | "<< outFluxN2 <<" | "<< inFluxEnthaply <<" | "<< outFluxEnthaply <<" | "<<*/ inNeumannEnthalpy <<" | "<< outNeumannEnthaply  << std::endl;
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

    //change made
    std::ofstream outputFile_;

    ReactionRate rrate_;
    Scalar timeStepSize_;
};

} //end namespace

#endif
