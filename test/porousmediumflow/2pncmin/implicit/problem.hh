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

#include <dune/grid/yaspgrid.hh>

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

#include "spatialparams.hh"

namespace Dumux {
/*!
 * \ingroup TwoPNCMinTests
 * \brief Problem where water is injected in a for flushing precipitated salt clogging a gas reservoir.
 */
template <class TypeTag>
class DissolutionProblem;

namespace Properties {
NEW_TYPE_TAG(DissolutionTypeTag, INHERITS_FROM(TwoPNCMin));
NEW_TYPE_TAG(DissolutionBoxTypeTag, INHERITS_FROM(BoxModel, DissolutionTypeTag));
NEW_TYPE_TAG(DissolutionCCTpfaTypeTag, INHERITS_FROM(CCTpfaModel, DissolutionTypeTag));

// Set the grid type
SET_TYPE_PROP(DissolutionTypeTag, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(DissolutionTypeTag, Problem, DissolutionProblem<TypeTag>);

// Set fluid configuration
SET_PROP(DissolutionTypeTag, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::BrineAir<Scalar, Components::H2O<Scalar>>;
};

SET_PROP(DissolutionTypeTag, SolidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ComponentOne = Components::NaCl<Scalar>;
    using ComponentTwo = Components::Granite<Scalar>;
    static constexpr int numInertComponents = 1;
    using type = SolidSystems::CompositionalSolidPhase<Scalar, ComponentOne, ComponentTwo, numInertComponents>;
};

// Set the spatial parameters
SET_PROP(DissolutionTypeTag, SpatialParams)
{
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = DissolutionSpatialParams<FVGridGeometry, Scalar>;
};

//Set properties here to override the default property settings
SET_INT_PROP(DissolutionTypeTag, ReplaceCompEqIdx, 1); //!< Replace gas balance by total mass balance
SET_PROP(DissolutionTypeTag, Formulation)
{ static constexpr auto value = TwoPFormulation::p1s0; };

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

        // Indices of the phases
        liquidPhaseIdx = FluidSystem::liquidPhaseIdx,
        gasPhaseIdx = FluidSystem::gasPhaseIdx,

        // index of the solid phase
        sPhaseIdx = SolidSystem::comp0Idx,


        // Index of the primary component of G and L phase
        conti0EqIdx = Indices::conti0EqIdx,
        precipNaClEqIdx = Indices::conti0EqIdx + FluidSystem::numComponents,

        // Phase State
        bothPhases = Indices::bothPhases,

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

public:
    DissolutionProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        outerSalinity_          = getParam<Scalar>("Problem.OuterSalinity");
        temperature_            = getParam<Scalar>("Problem.Temperature");
        reservoirPressure_      = getParam<Scalar>("Problem.ReservoirPressure");
        initLiqSaturation_      = getParam<Scalar>("Problem.LiquidSaturation");
        initPrecipitatedSaltBlock_  = getParam<Scalar>("Problem.InitPrecipitatedSaltBlock");

        outerLiqSaturation_     = getParam<Scalar>("Problem.OuterLiqSaturation");
        innerLiqSaturation_     = getParam<Scalar>("Problem.InnerLiqSaturation");
        innerSalinity_          = getParam<Scalar>("Problem.InnerSalinity");
        innerPressure_          = getParam<Scalar>("Problem.InnerPressure");
        outerPressure_          = getParam<Scalar>("Problem.OuterPressure");
        reservoirSaturation_    = getParam<Scalar>("Problem.reservoirSaturation");

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

        const Scalar rmax = this->fvGridGeometry().bBoxMax()[0];
        const Scalar rmin = this->fvGridGeometry().bBoxMin()[0];

        // default to Neumann
        bcTypes.setAllNeumann();

        // Constant pressure  at reservoir boundary (Dirichlet condition)
        if(globalPos[0] > rmax - eps_)
            bcTypes.setAllDirichlet();

        // Constant pressure at well (Dirichlet condition)
        if(globalPos[0] < rmin + eps_)
            bcTypes.setAllDirichlet();

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

        const Scalar rmax = this->fvGridGeometry().bBoxMax()[0];
        const Scalar rmin = this->fvGridGeometry().bBoxMin()[0];

        if(globalPos[0] > rmax - eps_)
        {
            priVars[pressureIdx]   = outerPressure_ ; // Outer boundary pressure bar
            priVars[switchIdx]     = outerLiqSaturation_; // Saturation outer boundary
            priVars[xwNaClIdx]     = massToMoleFrac_(outerSalinity_);// mole fraction salt
            priVars[precipNaClIdx] = 0.0;// precipitated salt
        }

        if(globalPos[0] < rmin + eps_)
        {
            priVars[pressureIdx]   = innerPressure_ ; // Inner boundary pressure bar
            priVars[switchIdx]     = innerLiqSaturation_; // Saturation inner boundary
            priVars[xwNaClIdx]     = massToMoleFrac_(innerSalinity_);// mole fraction salt
            priVars[precipNaClIdx] = 0.0;// precipitated salt
        }

        return priVars;
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

        priVars[pressureIdx] = reservoirPressure_;
        priVars[switchIdx]   = initLiqSaturation_;                 // Sw primary variable
        priVars[xwNaClIdx]   = massToMoleFrac_(outerSalinity_);     // mole fraction
        if(globalPos[0] > 5.0 - eps_ && globalPos[0] < 19.0 + eps_)
            priVars[precipNaClIdx] = initPrecipitatedSaltBlock_; // [kg/m^3]
        else
            priVars[precipNaClIdx] = 0.0; // [kg/m^3]

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

    int nTemperature_;
    int nPressure_;
    std::string name_;

    Scalar pressureLow_, pressureHigh_;
    Scalar temperatureLow_, temperatureHigh_;
    Scalar outerSalinity_;
    Scalar reservoirPressure_;
    Scalar innerPressure_;
    Scalar outerPressure_;
    Scalar temperature_;
    Scalar initLiqSaturation_;
    Scalar outerLiqSaturation_;
    Scalar innerLiqSaturation_;
    Scalar initPrecipitatedSaltBlock_;
    Scalar innerSalinity_;
    Scalar time_ = 0.0;
    Scalar timeStepSize_ = 0.0;
    static constexpr Scalar eps_ = 1e-6;
    Scalar reservoirSaturation_;
    std::vector<double> permeability_;
};

} //end namespace Dumux

#endif
