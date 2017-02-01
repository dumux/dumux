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
 * \brief Problem where water is injected in a for flushing precipitated salt clogging a gas reservoir.
 */
#ifndef DUMUX_DISSOLUTION_PROBLEM_HH
#define DUMUX_DISSOLUTION_PROBLEM_HH

#include <dumux/porousmediumflow/2pncmin/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/material/fluidsystems/brineair.hh>

#include "dissolutionspatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class DissolutionProblem;

namespace Properties
{
NEW_TYPE_TAG(DissolutionProblem, INHERITS_FROM(TwoPNCMin, DissolutionSpatialparams));
NEW_TYPE_TAG(DissolutionBoxProblem, INHERITS_FROM(BoxModel, DissolutionProblem));
NEW_TYPE_TAG(DissolutionCCProblem, INHERITS_FROM(CCModel, DissolutionProblem));

// Set the grid type
SET_TYPE_PROP(DissolutionProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(DissolutionProblem, Problem, DissolutionProblem<TypeTag>);

// Set fluid configuration
SET_PROP(DissolutionProblem, FluidSystem)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef FluidSystems::BrineAir<Scalar, H2O<Scalar>, true/*useComplexrelations=*/> type;
};

// Set the spatial parameters
SET_TYPE_PROP(DissolutionProblem, SpatialParams, DissolutionSpatialparams<TypeTag>);

// Enable gravity
SET_BOOL_PROP(DissolutionProblem, ProblemEnableGravity, true);

//Set properties here to override the default property settings in the model.
SET_INT_PROP(DissolutionProblem, ReplaceCompEqIdx, 1);
SET_INT_PROP(DissolutionProblem, Formulation, TwoPNCFormulation::pnsw);
}

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
class DissolutionProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum {

        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx, //Saturation
        xlNaClIdx = FluidSystem::NaClIdx,
        precipNaClIdx = FluidSystem::numComponents,

        //Indices of the components
        wCompIdx = FluidSystem::H2OIdx,
        nCompIdx = FluidSystem::AirIdx,
        NaClIdx = FluidSystem::NaClIdx,

        //Indices of the phases
        wPhaseIdx = FluidSystem::wPhaseIdx,
        nPhaseIdx = FluidSystem::nPhaseIdx,
        sPhaseIdx = FluidSystem::sPhaseIdx,

        //Index of the primary component of G and L phase
        conti0EqIdx = Indices::conti0EqIdx,
        contiTotalMassIdx = conti0EqIdx + FluidSystem::AirIdx,
        precipNaClEqIdx = Indices::conti0EqIdx + FluidSystem::numComponents,
        contiWEqIdx = conti0EqIdx + FluidSystem::H2OIdx,

        // Phase State
        wPhaseOnly = Indices::wPhaseOnly,
        nPhaseOnly = Indices::nPhaseOnly,
        bothPhases = Indices::bothPhases,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
    };


    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    DissolutionProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {

        outerSalinity_          = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, OuterSalinity);
        temperature_            = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, Temperature);
        reservoirPressure_      = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, ReservoirPressure);
        initLiqSaturation_      = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, LiquidSaturation);
        initPrecipitatedSalt1_  = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, InitPrecipitatedSalt1);
        initPrecipitatedSalt2_  = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, InitPrecipitatedSalt2);

        outerLiqSaturation_     = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, OuterLiqSaturation);
        innerLiqSaturation_     = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, InnerLiqSaturation);
        innerSalinity_          = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, InnerSalinity);
        innerPressure_          = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, InnerPressure);
        outerPressure_          = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, OuterPressure);
        reservoirSaturation_    = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, reservoirSaturation);

        nTemperature_           = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, FluidSystem, NTemperature);
        nPressure_              = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, FluidSystem, NPressure);
        pressureLow_            = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, PressureLow);
        pressureHigh_           = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, PressureHigh);
        temperatureLow_         = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, TemperatureLow);
        temperatureHigh_        = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, TemperatureHigh);
        name_                   = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);
        freqMassOutput_         = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Output, FreqMassOutput);
        storageLastTimestep_    = 0.0;
        lastMassOutputTime_     = 0.0;

        outfile.open("evaporation.out");
        outfile << "time; evaporationRate" << std::endl;

        FluidSystem::init(/*Tmin=*/temperatureLow_,
                          /*Tmax=*/temperatureHigh_,
                          /*nT=*/nTemperature_,
                          /*pmin=*/pressureLow_,
                          /*pmax=*/pressureHigh_,
                          /*np=*/nPressure_);
    }

    ~DissolutionProblem()
    {
        outfile.close();
    }

    bool shouldWriteOutput() const
    {
        return this->timeManager().timeStepIndex() % 1 == 0 ||
               this->timeManager().episodeWillBeFinished() ||
               this->timeManager().willBeFinished();
    }


    /*!
     * \name Problem parameters
     */


    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string name() const
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
    void boundaryTypesAtPos(BoundaryTypes &bcTypes, const GlobalPosition &globalPos) const
    {
        const Scalar rmax = this->bBoxMax()[0]; // outerRadius_;
        const Scalar rmin = this->bBoxMin()[0];

        // default to Neumann
        bcTypes.setAllNeumann();

        // Constant pressure  at reservoir boundary (Dirichlet condition)
        if(globalPos[0] > rmax - eps_)
            bcTypes.setAllDirichlet();

        // Constant pressure at well (Dirichlet condition)
        if(globalPos[0] < rmin + eps_)
            bcTypes.setAllDirichlet();
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        const Scalar rmax = this->bBoxMax()[0];
        const Scalar rmin = this->bBoxMin()[0];

        if(globalPos[0] > rmax - eps_)
        {
            values[pressureIdx]   = outerPressure_ ; // Outer boundary pressure bar
            values[switchIdx]     = outerLiqSaturation_; // Saturation outer boundary
            values[xlNaClIdx]     = massTomoleFrac_(outerSalinity_);// mole fraction salt
            values[precipNaClIdx] = 0.0;// precipitated salt
        }

        if(globalPos[0] < rmin + eps_)
        {

            values[pressureIdx]   = innerPressure_ ; // Inner boundary pressure bar
            values[switchIdx]     = innerLiqSaturation_; // Saturation inner boundary
            values[xlNaClIdx]     = massTomoleFrac_(innerSalinity_);// mole fraction salt
            values[precipNaClIdx] = 0.0;// precipitated salt
        }
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     */
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const Intersection &is,
                 int scvIdx,
                 int boundaryFaceIdx) const
    {
        values = 0.0;
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    void initialAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        values[pressureIdx] = reservoirPressure_;
        values[switchIdx]   = initLiqSaturation_;                 // Sl primary variable
        values[xlNaClIdx]   = massTomoleFrac_(outerSalinity_);     // mole fraction
        if(globalPos[0] > 5.0 - eps_ && globalPos[0] < 19.0 + eps_)
            values[precipNaClIdx] = initPrecipitatedSalt2_; // [kg/m^3]
        else
            values[precipNaClIdx] = initPrecipitatedSalt1_; // [kg/m^3]
    }

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * For this method, the \a values parameter stores the rate mass
     * of a component is generated or annihilate per volume
     * unit. Positive values mean that mass is created, negative ones
     * mean that it vanishes.
     */
    void solDependentSource(PrimaryVariables &source,
                            const Element &element,
                            const FVElementGeometry &fvGeometry,
                            int scvIdx,
                            const ElementVolumeVariables &elemVolVars) const
    {
        source = 0;
        const auto& volVars = elemVolVars[scvIdx];
        Scalar moleFracNaCl_lPhase = volVars.moleFraction(wPhaseIdx, NaClIdx);
        Scalar moleFracNaCl_gPhase = volVars.moleFraction(nPhaseIdx, NaClIdx);
        Scalar massFracNaCl_Max_lPhase = this->spatialParams().SolubilityLimit();
        Scalar moleFracNaCl_Max_lPhase = massTomoleFrac_(massFracNaCl_Max_lPhase);
        Scalar moleFracNaCl_Max_gPhase = moleFracNaCl_Max_lPhase / volVars.pressure(nPhaseIdx);
        Scalar saltPorosity = this->spatialParams().porosityMin(element, fvGeometry, scvIdx);

        // liquid phase
        using std::abs;
        Scalar precipSalt = volVars.porosity() * volVars.molarDensity(wPhaseIdx)
                                               * volVars.saturation(wPhaseIdx)
                                               * abs(moleFracNaCl_lPhase - moleFracNaCl_Max_lPhase);

        if (moleFracNaCl_lPhase < moleFracNaCl_Max_lPhase)
            precipSalt *= -1;

        // gas phase
        if (moleFracNaCl_gPhase > moleFracNaCl_Max_gPhase)
            precipSalt += volVars.porosity() * volVars.molarDensity(nPhaseIdx)
                                             * volVars.saturation(nPhaseIdx)
                                             * abs(moleFracNaCl_gPhase - moleFracNaCl_Max_gPhase);

        // make sure we don't disolve more salt than previously precipitated
        if (precipSalt*this->timeManager().timeStepSize() + volVars.precipitateVolumeFraction(sPhaseIdx)* volVars.molarDensity(sPhaseIdx)< 0)
            precipSalt = - volVars.precipitateVolumeFraction(sPhaseIdx)* volVars.molarDensity(sPhaseIdx)/this->timeManager().timeStepSize();

        if (volVars.precipitateVolumeFraction(sPhaseIdx) >= volVars.initialPorosity() - saltPorosity  && precipSalt > 0)
            precipSalt = 0;

        source[conti0EqIdx + NaClIdx] += -precipSalt;
        source[precipNaClEqIdx] += precipSalt;

        Valgrind::CheckDefined(source);
    }
    /*!
     * \brief Return the initial phase state inside a control volume.
     */
    int initialPhasePresence(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             int scvIdx) const
    {
        return bothPhases;
    }

private:

    /*!
     * \brief Returns the molality of NaCl (mol NaCl / kg water) for a given mole fraction
     *
     * \param XlNaCl the XlNaCl [kg NaCl / kg solution]
     */
    static Scalar massTomoleFrac_(Scalar XlNaCl)
    {
       const Scalar Mw = 18.015e-3; /* molecular weight of water [kg/mol] */
       const Scalar Ms = 58.44e-3; /* molecular weight of NaCl  [kg/mol] */

       const Scalar X_NaCl = XlNaCl;
       /* XlNaCl: conversion from mass fraction to mol fraction */
       auto xlNaCl = -Mw * X_NaCl / ((Ms - Mw) * X_NaCl - Ms);
       return xlNaCl;
    }

    int nTemperature_;
    int nPressure_;
    int freqMassOutput_;
    PrimaryVariables storageLastTimestep_;
    Scalar lastMassOutputTime_;
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
    Scalar initPrecipitatedSalt1_;
    Scalar initPrecipitatedSalt2_;
    Scalar innerSalinity_;
    static constexpr Scalar eps_ = 1e-6;
    Scalar reservoirSaturation_;
    std::ofstream outfile;

};
} //end namespace

#endif
