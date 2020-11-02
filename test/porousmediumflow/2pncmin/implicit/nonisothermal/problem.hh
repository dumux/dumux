// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \brief Problem where brine is evaporating at the top boundary. The system is closed at the remaining boundaries.
 */
#ifndef DUMUX_SALINIZATION_PROBLEM_HH
#define DUMUX_SALINIZATION_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/2pncmin/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/fluidsystems/brineair.hh>
#include <dumux/io/gnuplotinterface.hh>

#include <dumux/material/components/nacl.hh>
#include <dumux/material/components/granite.hh>
#include <dumux/material/solidsystems/compositionalsolidphase.hh>

#include "spatialparams.hh"

namespace Dumux {
/*!
 * \ingroup TwoPNCMinTests
 * \brief Problem where brine is evaporating at the top boundary. The system is closed at the remaining boundaries.
 */
template <class TypeTag>
class SalinizationProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct Salinization { using InheritsFrom = std::tuple<TwoPNCMinNI>; };
struct SalinizationBox { using InheritsFrom = std::tuple<Salinization, BoxModel>; };
struct SalinizationCCTpfa { using InheritsFrom = std::tuple<Salinization, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Salinization>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<Scalar, 2>>;
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Salinization> { using type = SalinizationProblem<TypeTag>; };

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Salinization>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::BrineAir<Scalar, Components::H2O<Scalar>>;
};

template<class TypeTag>
struct SolidSystem<TypeTag, TTag::Salinization>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ComponentOne = Components::NaCl<Scalar>;
    using ComponentTwo = Components::Granite<Scalar>;
    static constexpr int numInertComponents = 1;
    using type = SolidSystems::CompositionalSolidPhase<Scalar, ComponentOne, ComponentTwo, numInertComponents>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Salinization>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = SalinizationSpatialParams<GridGeometry, Scalar>;
};

// Set properties here to override the default property settings
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::Salinization> { static constexpr int value = 1; }; //!< Replace gas balance by total mass balance
template<class TypeTag>
struct Formulation<TypeTag, TTag::Salinization>
{ static constexpr auto value = TwoPFormulation::p0s1; };

} // end namespace Properties

/*!
 * \ingroup TwoPNCMinTests
 * \brief Problem where water is evaporating at the top of a porous media filled
 * container saturated with brine and air, which causes precipitation at the top.
 *
 * Problem with a porous media in a container, which is open to the atmosphere
 * at the top boundary. The container has dimensions of 0.2m by 0.2m. Neumann
 * no-flow boundaries are applied at the left, right and bottom boundary.
 * The grid is refined towards the upper boundary to capter the relevant processes.
 *
 * Initially the porous medium is 50 % saturated with brine. Evaporation takes
 * place at the top boundary and hence the temperature and liquid saturation
 * decreases first at the top, then in the whole system, whereas the sodium
 * chloride (NaCl) concentration increases. This results in precipitaion of NaCl
 * at the top as the solubility limit is exceeded. Due to the low liquid saturation
 * the top after some time, top temperature rises again.
 *
 * The model uses mole fractions of dissolved components and volume fractions of
 * precipitated salt as primary variables. Make sure that the according units
 * are used in the problem set-up.
 *
 * This problem uses the \ref TwoPNCMinModel.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_2pncminni_salinization</tt>
 */
template <class TypeTag>
class SalinizationProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;

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

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
    };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GlobalPosition = typename SubControlVolume::GlobalPosition;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

public:
    SalinizationProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        //Fluidsystem
        nTemperature_           = getParam<int>("FluidSystem.NTemperature");
        nPressure_              = getParam<int>("FluidSystem.NPressure");
        pressureLow_            = getParam<Scalar>("FluidSystem.PressureLow");
        pressureHigh_           = getParam<Scalar>("FluidSystem.PressureHigh");
        temperatureLow_         = getParam<Scalar>("FluidSystem.TemperatureLow");
        temperatureHigh_        = getParam<Scalar>("FluidSystem.TemperatureHigh");
        name_                   = getParam<std::string>("Problem.Name");

        //problem
        name_ = getParam<std::string>("Problem.Name");
        temperature_            = getParam<Scalar>("Problem.Temperature");

        //inital conditions
        initPressure_      = getParam<Scalar>("Problem.InitialPressure");
        initGasSaturation_      = getParam<Scalar>("Problem.InitialGasSaturation");
        initSalinity_          = getParam<Scalar>("Problem.InitialSalinity");

        unsigned int codim = GetPropType<TypeTag, Properties::GridGeometry>::discMethod == DiscretizationMethod::box ? dim : 0;

        permeability_.resize(gridGeometry->gridView().size(codim));
        FluidSystem::init(/*Tmin=*/temperatureLow_,
                          /*Tmax=*/temperatureHigh_,
                          /*nT=*/nTemperature_,
                          /*pmin=*/pressureLow_,
                          /*pmax=*/pressureHigh_,
                          /*np=*/nPressure_);
    }

    /*!
     * \brief The current time.
     */
    void setTime( Scalar time )
    {
        time_ = time;
    }

    /*!
     * \brief The time step size.
     *
     * This is used to calculate the source term.
     */
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

        // default to Neumann
        bcTypes.setAllNeumann();

        return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars.setState(bothPhases);

        return priVars;
    }


    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * This is the method for the case where the Neumann condition is
     * potentially solution dependent
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param elemFluxVarsCache Flux variables caches for all faces in stencil
     * \param scvf The sub control volume face
     *
     * Negative values mean influx.
     * E.g. for the mass balance that would be the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        PrimaryVariables values(0.0);

        const auto& globalPos = scvf.ipGlobal();
        const auto& volVars = elemVolVars[scvf.insideScvIdx()];
        const Scalar hmax = this->gridGeometry().bBoxMax()[1];

        static const Scalar temperatureRef = getParam<Scalar>("FreeFlow.RefTemperature");

        if (globalPos[1] > hmax - eps_)
        {
            // get free-flow properties:
            static const Scalar moleFracRefH2O = getParam<Scalar>("FreeFlow.RefMoleFracH2O");
            static const Scalar boundaryLayerThickness = getParam<Scalar>("FreeFlow.BoundaryLayerThickness");
            static const Scalar massTransferCoefficient = getParam<Scalar>("FreeFlow.MassTransferCoefficient");

            // get porous medium values:
            const Scalar moleFracH2OInside = volVars.moleFraction(gasPhaseIdx, H2OIdx);
            static const Scalar referencePermeability = getParam<Scalar>("SpatialParams.referencePermeability", 2.23e-14);

            // calculate fluxes
            // liquid phase
            Scalar evaporationRateMole = 0;
            if (moleFracH2OInside - moleFracRefH2O > 0)
            {
                evaporationRateMole = massTransferCoefficient
                                        * volVars.diffusionCoefficient(gasPhaseIdx, AirIdx, H2OIdx)
                                        * (moleFracH2OInside - moleFracRefH2O)
                                        / boundaryLayerThickness
                                        * volVars.molarDensity(gasPhaseIdx);
            }
            else
            {
                evaporationRateMole = massTransferCoefficient
                                        * volVars.diffusionCoefficient(gasPhaseIdx, AirIdx, H2OIdx)
                                        * (moleFracH2OInside - moleFracRefH2O)
                                        / boundaryLayerThickness
                                        * 1.2;

            }

            values[conti0EqIdx] = evaporationRateMole;

            // gas phase
            // gas flows in
            if (volVars.pressure(gasPhaseIdx) - 1e5 > 0) {
                values[conti1EqIdx] = (volVars.pressure(gasPhaseIdx) - 1e5)
                                      /(globalPos - fvGeometry.scv(scvf.insideScvIdx()).center()).two_norm()
                                      *volVars.mobility(gasPhaseIdx)
                                      *referencePermeability
                                      *volVars.molarDensity(gasPhaseIdx)
                                      *volVars.moleFraction(gasPhaseIdx, AirIdx);
            }
            //gas flows out
            else {
                values[conti1EqIdx] = (volVars.pressure(gasPhaseIdx) - 1e5)
                                      /(globalPos - fvGeometry.scv(scvf.insideScvIdx()).center()).two_norm()
                                      *volVars.mobility(gasPhaseIdx)
                                      *referencePermeability
                                      *volVars.molarDensity(gasPhaseIdx) * (1-moleFracRefH2O);
            }

            // energy fluxes
            values[energyEqIdx] = FluidSystem::componentEnthalpy(volVars.fluidState(), gasPhaseIdx, H2OIdx) * values[conti0EqIdx] * FluidSystem::molarMass(H2OIdx);

            values[energyEqIdx] += FluidSystem::componentEnthalpy(volVars.fluidState(), gasPhaseIdx, AirIdx)* values[conti1EqIdx] * FluidSystem::molarMass(AirIdx);

            values[energyEqIdx] += FluidSystem::thermalConductivity(elemVolVars[scvf.insideScvIdx()].fluidState(), gasPhaseIdx) * (volVars.temperature() - temperatureRef)/boundaryLayerThickness;
        }
        return values;
    }

    /*!
     * \brief Evaluates the initial value for a control volume.
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

        priVars[pressureIdx] = initPressure_ - density*9.81*globalPos[dimWorld-1];
        priVars[switchIdx]   = initGasSaturation_;                 // Sg primary variable
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
     * \brief Evaluates the source term for all phases within a given
     *        sub-controlvolume.
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
     * Positive values mean that the conserved quantity is created, negative ones mean that it vanishes.
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
        Scalar massFracNaCl_Max_wPhase = this->spatialParams().solubilityLimit();
        Scalar moleFracNaCl_Max_wPhase = massToMoleFrac_(massFracNaCl_Max_wPhase);
        Scalar saltPorosity = this->spatialParams().minimalPorosity(element, scv);

        // precipitation of amount of salt whic hexeeds the solubility limit
        using std::abs;
        Scalar precipSalt = volVars.porosity() * volVars.molarDensity(liquidPhaseIdx)
                                               * volVars.saturation(liquidPhaseIdx)
                                               * abs(moleFracNaCl_wPhase - moleFracNaCl_Max_wPhase)
                                               / timeStepSize_;
        if (moleFracNaCl_wPhase < moleFracNaCl_Max_wPhase)
            precipSalt *= -1;

        // make sure we don't dissolve more salt than previously precipitated
        if (precipSalt*timeStepSize_ + volVars.solidVolumeFraction(sPhaseIdx)* volVars.solidComponentMolarDensity(sPhaseIdx)< 0)
            precipSalt = -volVars.solidVolumeFraction(sPhaseIdx)* volVars.solidComponentMolarDensity(sPhaseIdx)/timeStepSize_;

        // make sure there is still pore space available for precipitation
        if (volVars.solidVolumeFraction(sPhaseIdx) >= this->spatialParams().referencePorosity(element, scv) - saltPorosity  && precipSalt > 0)
            precipSalt = 0;

        source[conti0EqIdx + NaClIdx] += -precipSalt;
        source[precipNaClEqIdx] += precipSalt;
        return source;
    }

    /*!
     * \brief Adds additional VTK output data to the VTKWriter.
     *
     * Function is called by the output module on every write.
     */

    const std::vector<Scalar>& getPermeability()
    {
        return permeability_;
    }

    void updateVtkOutput(const SolutionVector& curSol)
    {
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            const auto elemSol = elementSolution(element, curSol, this->gridGeometry());

            auto fvGeometry = localView(this->gridGeometry());
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

    Scalar extrusionFactorAtPos(const GlobalPosition& globalPos) const
    {
        return 0.054977871437821;
    }

private:

    /*!
     * \brief Returns the molality of NaCl (mol NaCl / kg water) for a given mole fraction.
     *
     * \param XwNaCl The XwNaCl [kg NaCl / kg solution]
     */
    static Scalar massToMoleFrac_(Scalar XwNaCl)
    {
       const Scalar Mw = 18.015e-3; //FluidSystem::molarMass(H2OIdx); /* molecular weight of water [kg/mol] */
       const Scalar Ms = 58.44e-3;  //FluidSystem::molarMass(NaClIdx); /* molecular weight of NaCl  [kg/mol] */

       const Scalar X_NaCl = XwNaCl;
       /* XwNaCl: conversion from mass fraction to mol fraction */
       auto xwNaCl = -Mw * X_NaCl / ((Ms - Mw) * X_NaCl - Ms);
       return xwNaCl;
    }


    std::string name_;

    Scalar initPressure_;
    Scalar initGasSaturation_;
    Scalar initSalinity_;

    Scalar temperature_;

    Scalar pressureLow_, pressureHigh_;
    Scalar temperatureLow_, temperatureHigh_;
    int nTemperature_;
    int nPressure_;

    Scalar time_ = 0.0;
    Scalar timeStepSize_ = 0.0;
    static constexpr Scalar eps_ = 1e-6;

    std::vector<Scalar> permeability_;

    Dumux::GnuplotInterface<Scalar> gnuplot_;
    Dumux::GnuplotInterface<Scalar> gnuplot2_;
    std::vector<Scalar> x_;
    std::vector<Scalar> y_;
    std::vector<Scalar> y2_;
};

} // end namespace Dumux

#endif
