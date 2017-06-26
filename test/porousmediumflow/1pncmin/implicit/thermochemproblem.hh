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
 * \brief Definition of a problem for water management in PEM fuel cells.
 */
#ifndef DUMUX_THERMOCHEM_PROBLEM_HH
#define DUMUX_THERMOCHEM_PROBLEM_HH

#include <dumux/porousmediumflow/1pncmin/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/material/fluidsystems/simplesteamaircao2h2.hh>
#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>

#include "thermochemspatialparams.hh"

#define NONISOTHERMAL 1


namespace Dumux
{

template <class TypeTag>
class ThermoChemProblem;

namespace Properties
{
#if NONISOTHERMAL
NEW_TYPE_TAG(ThermoChemProblem, INHERITS_FROM(OnePNCMinNI, ThermoChemSpatialParams));
NEW_TYPE_TAG(ThermoChemBoxProblem, INHERITS_FROM(BoxModel, ThermoChemProblem));
NEW_TYPE_TAG(ThermoChemCCProblem, INHERITS_FROM(CCModel, ThermoChemProblem));
#else
NEW_TYPE_TAG(ThermoChemProblem, INHERITS_FROM(OnePNCMin, ThermoChemSpatialParams));
NEW_TYPE_TAG(ThermoChemBoxProblem, INHERITS_FROM(BoxModel, ThermoChemProblem));
NEW_TYPE_TAG(ThermoChemCCProblem, INHERITS_FROM(CCModel, ThermoChemProblem));
#endif
// Set the grid type
SET_TYPE_PROP(ThermoChemProblem, Grid, Dune::YaspGrid<2>);
// Set the problem property
SET_TYPE_PROP(ThermoChemProblem, Problem, ThermoChemProblem<TypeTag>);
// Set fluid configuration
SET_PROP(ThermoChemProblem, FluidSystem)
{ /*private:*/
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::SteamAirCaO2H2<Scalar>;
};

// Set the transport equation that is replaced by the total mass balance
// SET_INT_PROP(ThermoChemProblem, ReplaceCompEqIdx, 1 /*3*/ /*1*/);

SET_TYPE_PROP(ThermoChemProblem, LinearSolver, UMFPackBackend<TypeTag>);

}

// Set the spatial parameters
SET_TYPE_PROP(ThermoChemProblem, SpatialParams, ThermoChemSpatialparams<TypeTag>);

/*!
 * \ingroup OnePNCModel
 * \ingroup ImplicitTestProblems
 * \brief Problem or water management in PEM fuel cells.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box2pnc</tt>
 */
template <class TypeTag>
class ThermoChemProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    using ParentType = ImplicitPorousMediaProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    enum
    {
        replaceCompEqIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx),

        numComponents = FluidSystem::numComponents,
        numSComponents = FluidSystem::numSComponents,

        // Indices of the primary variables
        pressureIdx = Indices::pressureIdx, //gas-phase pressure
        firstMoleFracIdx = Indices::firstMoleFracIdx, // mole fraction water

        CaOIdx = FluidSystem::numComponents,
        CaO2H2Idx = FluidSystem::numComponents+1,

        //Equation Indices
        conti0EqIdx = Indices::conti0EqIdx,
        firstTransportEqIdx = Indices::firstTransportEqIdx,
        solidEqIdx = Indices::conti0EqIdx + FluidSystem::numComponents,

        // Phase Indices
        phaseIdx = FluidSystem::gPhaseIdx,
        cPhaseIdx = FluidSystem::cPhaseIdx,
        hPhaseIdx = FluidSystem::hPhaseIdx,

#if NONISOTHERMAL
        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx
#endif
    };

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<dim>::Entity;
    using Intersection = typename GridView::Intersection;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using GridCreator = typename GET_PROP_TYPE(TypeTag, GridCreator);
    using DimVector = Dune::FieldVector<Scalar, dim> ;


    // Select the electrochemistry method
//     typedef typename Dumux::ElectroChemistry<TypeTag, Dumux::ElectroChemistryModel::Ochs> ElectroChemistry;
//     typedef Dumux::Constants<Scalar> Constant;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    ThermoChemProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        nTemperature_           = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, FluidSystem, NTemperature);
        nPressure_              = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, FluidSystem, NPressure);
        pressureLow_            = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, PressureLow);
        pressureHigh_           = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, PressureHigh);
        temperatureLow_         = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, TemperatureLow);
        temperatureHigh_        = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, TemperatureHigh);
        name_                   = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);


        FluidSystem::init(/*Tmin=*/temperatureLow_,
                          /*Tmax=*/temperatureHigh_,
                          /*nT=*/nTemperature_,
                          /*pmin=*/pressureLow_,
                          /*pmax=*/pressureHigh_,
                          /*np=*/nPressure_);
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
     *        used for which equation on a given boundary segment
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos( const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        values.setAllNeumann();

        if(globalPos[0] < eps_ )
        {
//             values.setDirichlet(pressureIdx);
            values.setDirichlet(firstMoleFracIdx);
            values.setDirichlet(temperatureIdx);
        }

        if (globalPos[0] > this->bBoxMax()[0] - eps_){
            values.setDirichlet(pressureIdx);
//             values.setDirichlet(firstMoleFracIdx);
//             values.setDirichlet(temperatureIdx);
            values.setOutflow(temperatureIdx);
            values.setOutflow(firstMoleFracIdx);
//             values.setDirichlet(firstMoleFracIdx);
        }
        return values;
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
      PrimaryVariables priVars(0.0);

      //input parameters
      Scalar pIn;
      Scalar pOut;
      Scalar tIn;
      Scalar tOut;
      Scalar vapor;

      // read input parameters
      if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, IsCharge)){
      pIn = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Charge, PressureIn);
      pOut = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Charge, PressureOut);
      tIn = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Charge, TemperatureIn);
      tOut = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Charge, TemperatureOut);
      vapor = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Charge, VaporIn);

      }
      else{
      pIn = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Discharge, PressureIn);
      pOut = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Discharge, PressureOut);
      tIn = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Discharge, TemperatureIn);
      tOut = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Discharge, TemperatureOut);
      vapor = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Discharge, VaporIn);
      }

        if(globalPos[0] < eps_)
        {
//           priVars[pressureIdx]   = pIn;
            priVars[firstMoleFracIdx]     = vapor; // Saturation outer boundary
            priVars[temperatureIdx] = tIn;
        }
        if(globalPos[0] > this->bBoxMax()[0] - eps_)
        {
            priVars[pressureIdx] = pOut;
//             priVars[firstMoleFracIdx] = 0.01; // Saturation inner boundary
//             priVars[temperatureIdx] = tOut;
//             priVars[firstMoleFracIdx]     = 0;
        }
        return priVars;
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

        PrimaryVariables neumannAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables priVars(0.0);

        if(globalPos[0] < eps_)
        {
        //if(isCharge == true){
        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, IsCharge)){
            priVars[pressureIdx] = -GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Charge, InFlow); //[mol/s] gas inflow; negative sign: inflow
        }
        else
        priVars[pressureIdx] = -GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Discharge, InFlow); //[mol/s] gas inflow
        }

        return priVars;
    }

    /*!
     * \name Volume terms
     */

    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * \param values Stores the initial values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
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

        //if(isCharge == true){
        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, IsCharge)){
        std::cout << "true " << "\n";
        pInit = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Charge, PressureInitial);
        tInit = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Charge, TemperatureInitial);
        h2oInit = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Charge, VaporInitial);
        CaOInit = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Charge, CaOInitial);
        CaO2H2Init = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Charge, CaO2H2Initial);
        }

        else {
        std::cout << "false " << "\n";
        pInit = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Discharge, PressureInitial);
        tInit = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Discharge, TemperatureInitial);
        h2oInit = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Discharge, VaporInitial);
        CaOInit = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Discharge, CaOInitial);
        CaO2H2Init = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Discharge, CaO2H2Initial);
        }

        priVars[pressureIdx] = pInit;
        priVars[firstMoleFracIdx]   = h2oInit;
#if NONISOTHERMAL
        priVars[temperatureIdx] = tInit;
#endif
        priVars[CaOIdx] = CaOInit;
        priVars[CaO2H2Idx]   = CaO2H2Init;

        return priVars;
    }

    /*!
     * \brief Return the initial phase state inside a sub control volume.
     *
     * \param element The element of the sub control volume
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The sub control volume index
     */
    PrimaryVariables solDependentSource(const Element &element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume &scv) const
    {

        PrimaryVariables source(0.0);
        const auto& volVars = elemVolVars[scv];

        Scalar Initial_Precipitate_Volume;
        Scalar maxPrecipitateVolumeCaO;
        Scalar maxPrecipitateVolumeCaO2H2;

        //if(isCharge == true){
        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, IsCharge)){
        Initial_Precipitate_Volume = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Charge, CaO2H2Initial);
        maxPrecipitateVolumeCaO =  0.3960;
        maxPrecipitateVolumeCaO2H2 = Initial_Precipitate_Volume;
        }
        else {
        Initial_Precipitate_Volume = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Discharge, CaOInitial);
        maxPrecipitateVolumeCaO = Initial_Precipitate_Volume;
        maxPrecipitateVolumeCaO2H2 = 0.3960;
        }

        Scalar T= volVars.temperature();
        Scalar Teq = 0;

        Scalar moleFractionVapor = 1e-3;

        if(volVars.moleFraction(firstMoleFracIdx) > 1e-3)
            moleFractionVapor = volVars.moleFraction(firstMoleFracIdx) ;
        if(volVars.moleFraction(firstMoleFracIdx) >= 1.0){
            moleFractionVapor = 1;
//           std::cout << " test vapor = " << "\n";
        }
        Scalar vaporPressure = volVars.pressure(phaseIdx) *moleFractionVapor ;
        vaporPressure *= 1.0e-5;
        Scalar pFactor = log(vaporPressure);

        Teq = -12845;
        Teq /= (pFactor - 16.508);        //the equilibrium temperature

//         if(isCharge == false) {
//             if (Teq < 573.15) {
//                 std::cout << "Teq = " << Teq<< "\n";
//                 // Teq=573.15;
//             }
//         }

        Scalar moleFracH2O_fPhase = volVars.moleFraction(firstMoleFracIdx);

        Scalar moleFracCaO_sPhase = volVars.precipitateVolumeFraction(cPhaseIdx)*volVars.molarDensity(cPhaseIdx)
                                     /(volVars.precipitateVolumeFraction(hPhaseIdx)*volVars.molarDensity(hPhaseIdx)
                                        + volVars.precipitateVolumeFraction(cPhaseIdx)*volVars.molarDensity(cPhaseIdx));

        Scalar moleFracCaO2H2_sPhase = volVars.precipitateVolumeFraction(hPhaseIdx)*volVars.molarDensity(hPhaseIdx)
                                       /(volVars.precipitateVolumeFraction(hPhaseIdx)*volVars.molarDensity(hPhaseIdx)
                                          + volVars.precipitateVolumeFraction(cPhaseIdx)*volVars.molarDensity(cPhaseIdx));

        Scalar deltaH = 112e3; // J/mol

        Scalar solidDensityAverage = moleFracCaO_sPhase*volVars.molarDensity(cPhaseIdx)
                                     + moleFracCaO2H2_sPhase* volVars.molarDensity(hPhaseIdx);

//         //discharge or hydration
        if (T < Teq){

            Scalar krh =0.08;  //0.006

            Scalar rHydration = - moleFracH2O_fPhase* (volVars.molarDensity(hPhaseIdx)- solidDensityAverage)
                                                     * krh * (T-Teq)/ Teq;

            Scalar q = - rHydration ;

            // make sure not more CaO reacts than present

            if (- q*this->timeManager().timeStepSize() + moleFracCaO_sPhase* volVars.molarDensity(cPhaseIdx) < 0 + eps_){
                q = moleFracCaO_sPhase/this->timeManager().timeStepSize();
//                 std::cout << "q_discharge = " << q << "\n";
            }

            source[conti0EqIdx+CaO2H2Idx] = q;

            source[conti0EqIdx+CaOIdx] = -q;

            source[conti0EqIdx+firstMoleFracIdx] = -q;

#if NONISOTHERMAL
            source[energyEqIdx] = q * deltaH;
#endif
        }

        // charge or dehydration
        else if(T > Teq){

            Scalar krd = 0.03; //0.05;

            Scalar rDehydration = (volVars.molarDensity(cPhaseIdx)- solidDensityAverage)
                                                     * krd * (Teq-T)/ Teq;

            Scalar q =  -rDehydration;

            if (- q*this->timeManager().timeStepSize() + moleFracCaO2H2_sPhase*volVars.molarDensity(hPhaseIdx) < 0){
                q = moleFracCaO2H2_sPhase/this->timeManager().timeStepSize();
            // std::cout << "q_charge " << q << "\n";
            }

            source[conti0EqIdx+CaO2H2Idx] = -q;

            source[conti0EqIdx+CaOIdx] = q;

            source[conti0EqIdx+firstMoleFracIdx] = q;
#if NONISOTHERMAL
            source[energyEqIdx] = -q * deltaH;

#endif
        }

        else
        source = 0.0;

        return source;
    }


    /*!
     * \brief Add problem specific vtk output for the electrochemistry
     */
//     void addOutputVtkFields()
//     {
//         // add the output field specific to the electrochemistry
//         typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;
//
//         // get the number of degrees of freedom
// //         auto numDofs = this->model().numDofs();
//
//         for (const auto& element : elements(this->gridView()))
//         {
//             FVElementGeometry fvGeometry;
//             fvGeometry.update(this->gridView(), element);
//
//             ElementVolumeVariables elemVolVars;
//             elemVolVars.update(*this,
//                                element,
//                                fvGeometry,
//                                false /* oldSol? */);
//
// //             for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
// //             {
// //                 const auto& globalPos = isBox ? element.geometry().corner(scvIdx)
// //                                               : element.geometry().center();
// //
// //                 auto dofIdxGlobal = this->model().dofMapper().subIndex(element, scvIdx, dofCodim);
// //             }
//         }
//     }

private:

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

//     PrimaryVariables storageLastTimestep_;

    std::string name_;

    Scalar pressureLow_, pressureHigh_;
    Scalar temperatureLow_, temperatureHigh_;
    Scalar reservoirPressure_;
    Scalar innerPressure_;
    Scalar outerPressure_;
    Scalar temperature_;

    static constexpr Scalar eps_ = 1e-6;
    Scalar reservoirSaturation_;
    std::ofstream outfile;

//     bool isCharge_ ;

};

} //end namespace

#endif
