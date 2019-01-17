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
 * \ingroup TwoPNCMinTests
 * \brief Problem where water is injected in a for flushing precipitated salt clogging a gas reservoir.
 */
#ifndef DUMUX_DISSOLUTION_PROBLEM_HH
#define DUMUX_DISSOLUTION_PROBLEM_HH

#include<stdio.h>
// #include<conio.h>
#include <dune/grid/yaspgrid.hh>
#include <iostream>

#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/box/properties.hh>
#include <dumux/porousmediumflow/2pncmin/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/fluidsystems/brineair.hh>

#include <dumux/material/components/nacl.hh>
#include <dumux/material/components/granite.hh>
#include <dumux/material/solidsystems/compositionalsolidphase.hh>

#include "salinizationspatialparams.hh"

namespace Dumux {
/*!
 * \ingroup TwoPNCMinTests
 * \brief Problem where water is injected in a for flushing precipitated salt clogging a gas reservoir.
 */
template <class TypeTag>
class DissolutionProblem;

namespace Properties {
NEW_TYPE_TAG(Dissolution, INHERITS_FROM(TwoPNCMinNI));
NEW_TYPE_TAG(DissolutionBox, INHERITS_FROM(BoxModel, Dissolution));
NEW_TYPE_TAG(DissolutionCCTpfa, INHERITS_FROM(CCTpfaModel, Dissolution));

// Set the grid type
// SET_TYPE_PROP(Dissolution, Grid, Dune::YaspGrid<2>);
SET_TYPE_PROP(Dissolution, Grid, Dune::YaspGrid<2, Dune::TensorProductCoordinates<double, 2> >);

// Set the problem property
SET_TYPE_PROP(Dissolution, Problem, DissolutionProblem<TypeTag>);

// Set fluid configuration
SET_PROP(Dissolution, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::BrineAir<Scalar, Components::H2O<Scalar>>;
};

SET_PROP(Dissolution, SolidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ComponentOne = Components::NaCl<Scalar>;
    using ComponentTwo = Components::Granite<Scalar>;
    static constexpr int numInertComponents = 1;
    using type = SolidSystems::CompositionalSolidPhase<Scalar, ComponentOne, ComponentTwo, numInertComponents>;
};

// Set the spatial parameters
SET_PROP(Dissolution, SpatialParams)
{
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = DissolutionSpatialParams<FVGridGeometry, Scalar>;
};

//Set properties here to override the default property settings
SET_INT_PROP(Dissolution, ReplaceCompEqIdx, 5); //!< Replace gas balance by total mass balance
SET_PROP(Dissolution, Formulation)
{ static constexpr auto value = TwoPFormulation::p0s1; };

} // end namespace Properties

/*!
 * \ingroup TwoPNCMinModel
 * \ingroup ImplicitTestProblems
 * \brief Problem where water is injected to flush precipitated salt in a gas reservoir clogged due to precipitated salt.
 *
 * The domain is sized 10m times 20m and contains a vertical low-permeable layer of precipitated salt near an extraction well.
 *
 * To flush this precipitated salt, water is injected through the gas extraction well in order to dissolve the precipitated salt increasing the permeability and thereby achieving high gas extraction rates later. Here, the system is assumed to be isothermal.
 * Neumann no-flow boundary condition is applied at the top and bottom boundary and Dirichlet boundary condition is used on the right and left sides.
 * The injected water phase migrates downwards due to increase in density as the precipitated salt dissolves.
 *
 * The model uses mole fractions of dissolved components and volume fractions of precipitated salt as primary variables. Make sure that the according units are used in the problem setup.
 *
 * This problem uses the \ref TwoPNCMinModel.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box2pncmin</tt>
 */
template <class TypeTag>
class DissolutionProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using SolidSystem = typename GET_PROP_TYPE(TypeTag, SolidSystem);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);

    enum
    {
        // primary variable indices
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,

        // component indices
        // TODO: using xwNaClIdx as privaridx works here, but
        //       looks like magic. Can this be done differently??
        xwNaClIdx = FluidSystem::NaClIdx,
        precipNaClIdx = FluidSystem::numComponents,

        // Indices of the components
        H2OIdx = FluidSystem::H2OIdx,
        NaClIdx = FluidSystem::NaClIdx,
        AirIdx = FluidSystem::AirIdx,

        // Indices of the phases
        liquidPhaseIdx = FluidSystem::liquidPhaseIdx,
        gasPhaseIdx = FluidSystem::gasPhaseIdx,

        // index of the solid phase
        sPhaseIdx = SolidSystem::comp0Idx,


        // Index of the primary component of G and L phase
        conti0EqIdx = Indices::conti0EqIdx, //water component
        conti1EqIdx = Indices::conti0EqIdx + 1, //air component
        precipNaClEqIdx = Indices::conti0EqIdx + FluidSystem::numComponents,
        energyEqIdx = Indices::energyEqIdx,

        // Phase State
        bothPhases = Indices::bothPhases,
        liquidPhaseOnly = Indices::firstPhaseOnly,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
    };

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GlobalPosition = typename SubControlVolume::GlobalPosition;
    using NeumannFluxes = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);

public:
    DissolutionProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {

        temperature_            = getParam<Scalar>("Problem.Temperature");
        initPressure_      = getParam<Scalar>("Problem.InitialPressure");
        initLiqSaturation_      = getParam<Scalar>("Problem.InitialLiquidSaturation");
        initSalinity_          = getParam<Scalar>("Problem.InitialSalinity");
        bottomLiqSaturation_     = getParam<Scalar>("Problem.BottomLiqSaturation");
        bottomSalinity_          = getParam<Scalar>("Problem.BottomSalinity");
        bottomPressure_          = getParam<Scalar>("Problem.BottomPressure");
        bottomTemperature_       = getParam<Scalar>("Problem.BottomTemperature");
        evaporationStartTime_       = getParam<Scalar>("Problem.EvaporationStartTime");
        waterDepth_       = getParam<Scalar>("Problem.WaterDepth");

        nTemperature_           = getParam<int>("FluidSystem.NTemperature");
        nPressure_              = getParam<int>("FluidSystem.NPressure");
        pressureLow_            = getParam<Scalar>("FluidSystem.PressureLow");
        pressureHigh_           = getParam<Scalar>("FluidSystem.PressureHigh");
        temperatureLow_         = getParam<Scalar>("FluidSystem.TemperatureLow");
        temperatureHigh_        = getParam<Scalar>("FluidSystem.TemperatureHigh");
        name_                   = getParam<std::string>("Problem.Name");

        unsigned int codim = GET_PROP_TYPE(TypeTag, FVGridGeometry)::discMethod == DiscretizationMethod::box ? dim : 0;
        permeability_.resize(fvGridGeometry->gridView().size(codim));

        FluidSystem::init(/*Tmin=*/temperatureLow_,
                          /*Tmax=*/temperatureHigh_,
                          /*nT=*/nTemperature_,
                          /*pmin=*/pressureLow_,
                          /*pmax=*/pressureHigh_,
                          /*np=*/nPressure_);

        this->spatialParams().plotMaterialLaw();

    }

    void setTime( Scalar time )
    {
        time_ = time;
    }

    void setTimeStepSize( Scalar timeStepSize )
     {
        timeStepSize_ = timeStepSize;
     }

    /*!
     * \name Problem parameters
     */


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
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return temperature_; }


           //! Called after every time step
    //! Output the total evaporation rate
    void postTimeStep(const SolutionVector& curSol,
                      const GridVolumeVariables &curGridVolVars)
    {
        PrimaryVariables source(0.0);
        Scalar evaporation = 0.0;
        Scalar saturation = 0.0;
        Scalar i = 0.0;

       for (const auto& element : elements(this->fvGridGeometry().gridView()))
       {
           auto fvGeometry = localView(this->fvGridGeometry());
           fvGeometry.bindElement(element);

           auto elemVolVars = localView(curGridVolVars);
           elemVolVars.bindElement(element, fvGeometry, curSol);

           for (auto&& scvf : scvfs(fvGeometry))
           {
               if (scvf.boundary())
               {
                    evaporation += neumann(element, fvGeometry, elemVolVars, scvf)[conti0EqIdx]
                                        * scvf.area() * elemVolVars[scvf.insideScvIdx()].extrusionFactor();


               }
            saturation += elemVolVars[scvf.insideScvIdx()].saturation(liquidPhaseIdx);
            i = i+1;
           }
        }

        // convert to kg/s if using mole fractions
        evaporation = evaporation * FluidSystem::molarMass(H2OIdx);
        saturation /= i;
        std::cout << "Soil evaporation rate: " << evaporation << " kg/s." << '\n';

        //do a gnuplot
        x_.push_back(time_); // in seconds
        y_.push_back(evaporation);

        gnuplot_.resetPlot();
        gnuplot_.setXRange(0, time_);
        gnuplot_.setYRange(0, 2e-7);
        gnuplot_.setXlabel("time [s]");
        gnuplot_.setYlabel("kg/s");
        gnuplot_.addDataSetToPlot(x_, y_, "evaporation");


        gnuplot_.plot("");


        // compute the mass in the entire domain
        Scalar massNaCl = 0.0;

        // bulk elements
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(curGridVolVars);
            elemVolVars.bindElement(element, fvGeometry, curSol);

            for (auto&& scv : scvs(fvGeometry))
            {
                const auto& volVars = elemVolVars[scv];
                for(int phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
                {
                    massNaCl += volVars.massFraction(phaseIdx, FluidSystem::NaClIdx)*volVars.density(phaseIdx)
                    * scv.volume() * volVars.saturation(phaseIdx) * volVars.porosity() * volVars.extrusionFactor();
                }
            }
        }

        //do a gnuplot
        y2_.push_back(massNaCl);

        gnuplot2_.resetPlot();
        gnuplot2_.setXRange(0, time_);
        gnuplot2_.setYRange(0, 0.0001);
        gnuplot2_.setXlabel("time [s]");
        gnuplot2_.setYlabel("mass NaCl[kg]");

        gnuplot2_.addDataSetToPlot(x_, y2_, "mass NaCl");

        gnuplot2_.plot("");

        std::cout << std::setprecision(15) << "mass of NaCl is: " << massNaCl << std::endl;

    }


    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;

        const Scalar hmin = this->fvGridGeometry().bBoxMin()[1];

        // default to Neumann
        bcTypes.setAllNeumann();

    if(globalPos[1] < this->fvGridGeometry().bBoxMin()[1] + eps_)
    {
        bcTypes.setAllDirichlet();
    }

         return bcTypes;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars.setState(bothPhases);

        const Scalar hmin = this->fvGridGeometry().bBoxMin()[1];
        const Scalar xmin = this->fvGridGeometry().bBoxMin()[0];
        const Scalar xmax = this->fvGridGeometry().bBoxMax()[0];
        Scalar density = 1000.00;
        const auto g = this->gravityAtPos(globalPos)[dimWorld-1];

        priVars[pressureIdx]   = bottomPressure_ - density*9.81*globalPos[dimWorld-1]; // Bottom boundary pressure bar
        priVars[switchIdx]     = bottomLiqSaturation_; // Saturation bottom boundary
        priVars[xwNaClIdx]     = massToMoleFrac_(bottomSalinity_);// mole fraction salt
        priVars[precipNaClIdx] = 0.0;// precipitated salt
        priVars[energyEqIdx] = bottomTemperature_;// precipitated salt

        return priVars;
    }


     NeumannFluxes neumann(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolumeFace& scvf) const
    {
        PrimaryVariables values(0.0);
        const auto& globalPos = scvf.ipGlobal();
        const auto& volVars = elemVolVars[scvf.insideScvIdx()];
        const Scalar hmax = this->fvGridGeometry().bBoxMax()[1];
        Scalar temperatureRef = getParam<Scalar>("FreeFlow.RefTemperature");

        const auto episodeLength = getParam<Scalar>("TimeLoop.EpisodeLength");
        if(time_ > episodeLength){
        if(globalPos[1] > hmax - eps_)
        {

            // get free- flow properties:
            static const Scalar moleFracRefH2O = getParam<Scalar>("FreeFlow.RefMoleFracH2O");
            static const Scalar boundaryLayerThickness = getParam<Scalar>("FreeFlow.BoundaryLayerThickness");
            static const Scalar massTransferCoefficient = getParam<Scalar>("FreeFlow.MassTransferCoefficient");

            // get porous medium values:
            Scalar moleFracH2OInside = volVars.moleFraction(gasPhaseIdx, H2OIdx);
            Scalar referencePermeability_ = getParam<Scalar>("SpatialParams.referencePermeability", 2.23e-14);

            // calculate fluxes
            // liquid phase
            Scalar evaporationRateMole = 0;
            if(moleFracH2OInside - moleFracRefH2O > 0)
            {
                evaporationRateMole = massTransferCoefficient
                                        * volVars.diffusionCoefficient(gasPhaseIdx, H2OIdx)
                                        * (moleFracH2OInside - moleFracRefH2O)
                                        / boundaryLayerThickness
                                        * volVars.molarDensity(gasPhaseIdx);
            }
            else
            {
                evaporationRateMole = massTransferCoefficient
                                        * volVars.diffusionCoefficient(gasPhaseIdx, H2OIdx)
                                        * (moleFracH2OInside - moleFracRefH2O)
                                        / boundaryLayerThickness
                                        * 1.2;

            }

            values[conti0EqIdx] = evaporationRateMole;

//             gas phase
            // gas flows in
            if (volVars.pressure(gasPhaseIdx) - 1e5 > 0) {
                values[conti1EqIdx] =
//                                         (volVars.pressure(gasPhaseIdx) + volVars.molarDensity(gasPhaseIdx) * volVars.moleFraction(gasPhaseIdx, AirIdx) * 9.81 *globalPos[1]
//                                         - 1e5 - volVars.molarDensity(gasPhaseIdx) * (1-moleFracRefH2O) * 9.81 * (globalPos[1] + boundaryLayerThickness))
                                        (volVars.pressure(gasPhaseIdx) - 1e5)
                                      /(globalPos - fvGeometry.scv(scvf.insideScvIdx()).center()).two_norm()
                                       *volVars.mobility(gasPhaseIdx)
                                       *referencePermeability_
                                       *volVars.molarDensity(gasPhaseIdx)
                                       *volVars.moleFraction(gasPhaseIdx, AirIdx);
            }
            //gas flows out
            else {
                values[conti1EqIdx] =
//                                         (volVars.pressure(gasPhaseIdx) + volVars.molarDensity(gasPhaseIdx) * volVars.moleFraction(gasPhaseIdx, AirIdx) * 9.81 *globalPos[1]
//                                         - 1e5 - volVars.molarDensity(gasPhaseIdx) * (1-moleFracRefH2O) * 9.81 * (globalPos[1] + boundaryLayerThickness))
                                        (volVars.pressure(gasPhaseIdx) - 1e5)
                                        /(globalPos - fvGeometry.scv(scvf.insideScvIdx()).center()).two_norm()
                                        *volVars.mobility(gasPhaseIdx)
                                        *referencePermeability_
                                        *volVars.molarDensity(gasPhaseIdx) * (1-moleFracRefH2O);
             }


            // energy fluxes
            values[energyEqIdx] = FluidSystem::componentEnthalpy(volVars.fluidState(), gasPhaseIdx, H2OIdx) * values[conti0EqIdx] * FluidSystem::molarMass(H2OIdx);

            values[energyEqIdx] += FluidSystem::componentEnthalpy(volVars.fluidState(), gasPhaseIdx, AirIdx)* values[conti1EqIdx] * FluidSystem::molarMass(AirIdx);

            values[energyEqIdx] += FluidSystem::thermalConductivity(elemVolVars[scvf.insideScvIdx()].fluidState(), gasPhaseIdx) * (volVars.temperature() - temperatureRef)/boundaryLayerThickness;

        }
        }
        return values;
    }






    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The global position
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {

        PrimaryVariables priVars(0.0);
        priVars.setState(bothPhases);
        Scalar density = 1000.00; //FluidSystem::density(, liquidPhaseIdx);

        priVars[pressureIdx] = bottomPressure_ - density*9.81*globalPos[dimWorld-1];
        priVars[switchIdx]   = initLiqSaturation_;                 // Sw primary variable
        priVars[xwNaClIdx]   = massToMoleFrac_(initSalinity_);     // mole fraction
        priVars[precipNaClIdx] = 0.0; // [kg/m^3]
        priVars[energyEqIdx] = temperature_; // [K]

        return priVars;
    }

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * This is the method for the case where the source term is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param values The source and sink values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^3 \cdot s )] \f$
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

        Scalar moleFracNaCl_wPhase = volVars.moleFraction(liquidPhaseIdx, NaClIdx);
        Scalar moleFracNaCl_nPhase = volVars.moleFraction(gasPhaseIdx, NaClIdx);
        Scalar massFracNaCl_Max_wPhase = this->spatialParams().solubilityLimit();
        Scalar moleFracNaCl_Max_wPhase = massToMoleFrac_(massFracNaCl_Max_wPhase);
        Scalar moleFracNaCl_Max_nPhase = moleFracNaCl_Max_wPhase / volVars.pressure(gasPhaseIdx);
        Scalar saltPorosity = this->spatialParams().minimalPorosity(element, scv);

        // liquid phase
        using std::abs;
        Scalar precipSalt = volVars.porosity() * volVars.molarDensity(liquidPhaseIdx)
                                               * volVars.saturation(liquidPhaseIdx)
                                               * abs(moleFracNaCl_wPhase - moleFracNaCl_Max_wPhase);
        if (moleFracNaCl_wPhase < moleFracNaCl_Max_wPhase)
            precipSalt *= -1;

        // gas phase
        precipSalt += volVars.porosity() * volVars.molarDensity(gasPhaseIdx)
                                         * volVars.saturation(gasPhaseIdx)
                                         * abs(moleFracNaCl_nPhase - moleFracNaCl_Max_nPhase);

        // make sure we don't dissolve more salt than previously precipitated
        if (precipSalt*timeStepSize_ + volVars.solidVolumeFraction(sPhaseIdx)* volVars.solidComponentMolarDensity(sPhaseIdx)< 0)
            precipSalt = -volVars.solidVolumeFraction(sPhaseIdx)* volVars.solidComponentMolarDensity(sPhaseIdx)/timeStepSize_;

        if (volVars.solidVolumeFraction(sPhaseIdx) >= this->spatialParams().referencePorosity(element, scv) - saltPorosity  && precipSalt > 0)
            precipSalt = 0;

        source[conti0EqIdx + NaClIdx] += -precipSalt;
        source[precipNaClEqIdx] += precipSalt;
        return source;

    }

    /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */

    const std::vector<Scalar>& getPermeability()
    {
        return permeability_;
    }

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
                    permeability_[dofIdxGlobal] = volVars.permeability();
                }
            }
        }

    Scalar extrusionFactorAtPos(const GlobalPosition &globalPos) const
    {
        // As a default, i.e. if the user's problem does not overload
        // any extrusion factor method, return 1.0
        return 0.054977871437821;
    }



private:

    /*!
     * \brief Returns the molality of NaCl (mol NaCl / kg water) for a given mole fraction
     *
     * \param XwNaCl the XwNaCl [kg NaCl / kg solution]
     */
    static Scalar massToMoleFrac_(Scalar XwNaCl)
    {
       const Scalar Mw = 18.015e-3; //FluidSystem::molarMass(H2OIdx); /* molecular weight of water [kg/mol] */ //TODO use correct link to FluidSyswem later
       const Scalar Ms = 58.44e-3;  //FluidSystem::molarMass(NaClIdx); /* molecular weight of NaCl  [kg/mol] */

       const Scalar X_NaCl = XwNaCl;
       /* XwNaCl: conversion from mass fraction to mol fraction */
       auto xwNaCl = -Mw * X_NaCl / ((Ms - Mw) * X_NaCl - Ms);
       return xwNaCl;
    }


    std::string name_;

    Scalar initSalinity_;
    Scalar initPressure_;
    Scalar initLiqSaturation_;

    Scalar bottomPressure_;
    Scalar bottomLiqSaturation_;
    Scalar bottomSalinity_;
    Scalar bottomTemperature_;
    Scalar waterDepth_;
    Scalar evaporationStartTime_;

    Scalar pressureLow_, pressureHigh_;
    Scalar temperatureLow_, temperatureHigh_;
    int nTemperature_;
    int nPressure_;

    Scalar time_ = 0.0;
    Scalar timeStepSize_ = 0.0;
    static constexpr Scalar eps_ = 1e-6;

    std::vector<double> permeability_;

    Dumux::GnuplotInterface<double> gnuplot_;
    Dumux::GnuplotInterface<double> gnuplot2_;
    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> y2_;

    Scalar temperature_;


};

} //end namespace Dumux

#endif
