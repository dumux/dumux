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
#ifndef DUMUX_MICP_COLUMN_SIMPLE_CHEM_PROBLEM_HH
#define DUMUX_MICP_COLUMN_SIMPLE_CHEM_PROBLEM_HH

// ## Initial and boundary conditions (`problem.hh`)
//
// This file contains the __problem class__ which defines the initial and boundary
// conditions for the biominearalization example.
//
// [[content]]
//
// ### Include files
// [[codeblock]]
// This header contains the PorousMediumFlowProblem class
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/common/boundarytypes.hh> // for `BoundaryTypes`
#include <dumux/common/properties.hh> // GetPropType
#include <dumux/material/components/ammonia.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/io/container.hh>
#include <algorithm> // for std::transform
// [[/codeblock]]

// ### The problem class
// The problem class defines the boundary and initial conditions.
// As this is a biomineralization problem in porous media, we inherit from the porous-media-flow problem

// [[codeblock]]
namespace Dumux {

template <class TypeTag>
class MICPColumnProblemSimpleChemistry : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NH3 = Components::Ammonia<Scalar>; // for molar mass of Ammonia
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    // We define some indices for convenience to be used later when defining the initial and boundary conditions
    enum {
        numComponents = FluidSystem::numComponents,

        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx, //Saturation
        xwNaIdx = FluidSystem::NaIdx,
        xwClIdx = FluidSystem::ClIdx,
        xwCaIdx = FluidSystem::CaIdx,
        xwUreaIdx = FluidSystem::UreaIdx,
        xwO2Idx = FluidSystem::O2Idx,
        xwBiosubIdx = FluidSystem::GlucoseIdx,
        xwSuspendedBiomassIdx = FluidSystem::SuspendedBiomassIdx,
        phiBiofilmIdx = numComponents,
        phiCalciteIdx = numComponents + 1,

        //Indices of the components
        wCompIdx = FluidSystem::wCompIdx,
        nCompIdx = FluidSystem::nCompIdx,
        NaIdx = FluidSystem::NaIdx,
        ClIdx = FluidSystem::ClIdx,
        CaIdx = FluidSystem::CaIdx,
        UreaIdx = FluidSystem::UreaIdx,
        O2Idx = FluidSystem::O2Idx,
        BiosubIdx = FluidSystem::GlucoseIdx,
        SuspendedBiomassIdx = FluidSystem::SuspendedBiomassIdx,

        wPhaseIdx = FluidSystem::wPhaseIdx,
        conti0EqIdx = Indices::conti0EqIdx,

        // Phase State
        wPhaseOnly = Indices::firstPhaseOnly,
        bothPhases = Indices::bothPhases
    };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using BoundaryTypes = Dumux::BoundaryTypes<PrimaryVariables::size()>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    static constexpr int dim = GridView::dimension;
// [[/codeblock]]
//
// We define an enum class to distinguish between different injection processes.
// We will use these later in the definition of the Neumann boundary conditions.
    enum class InjectionProcess : int
    {
        noInjection = -99,
        noFlow = -9,
        rinse = -1,
        resuscitation = 3,
        inoculation = 2,
        mineralization = 1
    };
//
// ### Reading of parameters specified in the params.input file
// [[details]] parameters
// [[codeblock]]
public:
    // This is the constructor of our problem class.
    MICPColumnProblemSimpleChemistry(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        // We read the parameters from the params.input file.
        name_  = getParam<std::string>("Problem.Name");
        temperature_ = getParam<Scalar>("Problem.Temperature");

        // biomass parameters
        ca1_ = getParam<Scalar>("BioCoefficients.Ca1");
        ca2_ = getParam<Scalar>("BioCoefficients.Ca2");
        cd1_ = getParam<Scalar>("BioCoefficients.Cd1");
        dc0_ = getParam<Scalar>("BioCoefficients.Dc0");
        kmue_  = getParam<Scalar>("BioCoefficients.Kmue");
        f_ = getParam<Scalar>("BioCoefficients.F");
        ke_ = getParam<Scalar>("BioCoefficients.Ke");
        ks_ = getParam<Scalar>("BioCoefficients.Ks");
        yield_ = getParam<Scalar>("BioCoefficients.Yield");

        // ureolysis kinetic parameters
        kub_ = getParam<Scalar>("UreolysisCoefficients.Kub");
        kurease_ = getParam<Scalar>("UreolysisCoefficients.Kurease");
        ku_ = getParam<Scalar>("UreolysisCoefficients.Ku");

        // initial values
        densityW_ = getParam<Scalar>("Initial.DensityW");
        initPressure_ = getParam<Scalar>("Initial.Pressure");

        initxwTC_ = getParam<Scalar>("Initial.XwTC");
        initxwNa_ = getParam<Scalar>("Initial.XwNa");
        initxwCl_ = getParam<Scalar>("Initial.XwCl");
        initxwCa_ = getParam<Scalar>("Initial.XwCa");
        initxwUrea_ = getParam<Scalar>("Initial.XwUrea");
        initxwTNH_ = getParam<Scalar>("Initial.XwTNH");
        initxwO2_ = getParam<Scalar>("Initial.XwO2");
        initxwBiosub_ = getParam<Scalar>("Initial.XwSubstrate");
        initxwBiosusp_ = getParam<Scalar>("Initial.XwSuspendedBiomass");
        initCalcite_ = getParam<Scalar>("Initial.Calcite");
        initBiofilm_ = getParam<Scalar>("Initial.Biofilm");

        xwNaCorr_ = getParam<Scalar>("Initial.XwNaCorr");
        xwClCorr_ = getParam<Scalar>("Initial.XwClCorr");

        // injection values
        injQ_ = getParam<Scalar>("Injection.FlowRate");

        injTC_ = getParam<Scalar>("Injection.MassFracTC");
        injNa_ = getParam<Scalar>("Injection.ConcentrationNa");
        injCa_ = getParam<Scalar>("Injection.ConcentrationCa");
        injUrea_ = getParam<Scalar>("Injection.ConcentrationUrea");
        injTNH_ = getParam<Scalar>("Injection.ConcentrationTNH");
        injO2_ = getParam<Scalar>("Injection.ConcentrationO2");
        injSub_ = getParam<Scalar>("Injection.ConcentrationSubstrate");
        injBiosusp_= getParam<Scalar>("Injection.ConcentrationSuspendedBiomass");
        injNaCorr_ = getParam<Scalar>("Injection.ConcentrationNaCorr");
        // [[/codeblock]]
        // [[/details]]
        //
        // ### Reading of the temporal sequence of the different types of injection
        // Here, the number of injections and the corresponding parameters are defined based on what is specified in the input files.
        // [[codeblock]]
        // We get the number of injections and the injection data file name from params.input
        numInjections_ = getParam<int>("Injection.NumInjections");

        // We resize the permeability vector contaning the permeabilities for the additional output
        permeability_.resize(gridGeometry->numDofs());

        // We read from the injection data file which injection type we have in each episode.
        // We will use this in the Neumann boundary condition to set time dependend, changing boundary conditions.
        // We do this similarly to the episode ends in the main file.
        const auto injType = readFileToContainer<std::vector<int>>("injection_type.dat");
        // translate integer to InjectionProcess type
        std::transform(injType.begin(), injType.end(), std::back_inserter(injectionType_),
               [](int n){ return static_cast<InjectionProcess>(n); });


       // [[/codeblock]]

        // We check the injection data against the number of injections specified in the parameter file
        // and print an error message if the test fails
        // [[codeblock]]
        if (injectionType_.size() != numInjections_)
            DUNE_THROW(Dune::IOError, "numInjections from the parameterfile and the number of injection types "
                << "specified in the injection data file do not match!\n"
                << "numInjections from parameter file: " << numInjections_ << "\n"
                << "numInjTypes from injection data file: "<< injectionType_.size());

        // Initialize the fluidsystem
        FluidSystem::init(/*startTemp=*/temperature_ -5.0, /*endTemp=*/temperature_ +5.0, /*tempSteps=*/5,
                          /*startPressure=*/1e4, /*endPressure=*/1e6, /*pressureSteps=*/500);
    }
    // [[/codeblock]]

    // In the follwing, functions to set the time, time step size and the index of the episode
    // are declared which are used the time loop in main.cc
    // [[codeblock]]
    void setTime(const Scalar t)
    { time_ = t; }

    // We need the time step size to regularize the reactive source terms.
    void setTimeStepSize(const Scalar dt)
    { timeStepSize_ = dt; }

    // We need the episode index to choose the right Neumann boundary condition for each episode based on the injection data.
    void setEpisodeIdx(const Scalar epIdx)
    { episodeIdx_ = epIdx; }
    // [[/codeblock]]

    // Here, functions to return the injectionType, the name of the problem and the temperature
    // are defined
    // [[codeblock]]
    int injectionType(int episodeIdx) const
    { return injectionType_[episodeIdx]; }

    // Get the problem name. It is used as a prefix for files generated by the simulation.
    const std::string& name() const
    { return name_; }

    // Return the temperature
    Scalar temperature() const
    { return temperature_; }
    // [[/codeblock]]

    // #### Boundary conditions
    // With the following function we define the type of boundary conditions depending on the location. Two types of boundary conditions
    // can be specified: Dirichlet or Neumann boundary condition. On a Dirichlet boundary, the values of the
    // primary variables need to be fixed. On a Neumann boundary condition, values for derivatives need to be fixed.

    // [[codeblock]]
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        // We set all to Neumann, except for the top boundary, which is set to Dirichlet.
        const Scalar zmax = this->gridGeometry().bBoxMax()[dim - 1];
        bcTypes.setAllNeumann();
        if (globalPos[dim - 1] > zmax - eps_)
            bcTypes.setAllDirichlet();

        return bcTypes;
    }
    // [[/codeblock]]

    // We define the Dirichlet boundary conditions
    // [[codeblock]]
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars.setState(wPhaseOnly);
        // We recycle the initial conditions, but additionally enforce that substrate and oxygen necessary for biomass growth are zero.
        priVars = initial_(globalPos);
        priVars[xwBiosubIdx] = 0.0;
        priVars[xwO2Idx] = 0.0;
        return priVars;
    }
    // [[/codeblock]]

    // The injections can be described with a Neumann boundary condition.
    // Since we have multiple injections durng the whole simulation period, the Neumann boundary conditions in this case are time or rather episode depended.
    // [[codeblock]]
    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    {
        // We only have injection at the bottom, negative values for injection
        if (globalPos[dim - 1] > eps_)
            return NumEqVector(0.0);

        // We calculate the injected water velocity based on the volume flux injected into a 1 inch diameter column.
        const Scalar diameter = 0.0254;
        const Scalar waterFlux = injQ_/(3.14*diameter*diameter/4.); //[m/s]
        const auto injProcess = injectionType_[episodeIdx_];

        // We also have no-injection, no-flow periods, which are coded for by -99 or 9.
        // Thus, for no flow, we set the Neumann BC to zero for all components.
        if (injProcess == InjectionProcess::noInjection || injProcess == InjectionProcess::noFlow)
            return NumEqVector(0.0);

        // We set the injection BC to the basic rinse injection and later only change the BC for those components that are different.
        NumEqVector values(0.0);

        values[conti0EqIdx + wCompIdx] = -waterFlux * 996/FluidSystem::molarMass(wCompIdx);
        values[conti0EqIdx + nCompIdx] = -waterFlux * injTC_*996 /FluidSystem::molarMass(nCompIdx);
        values[conti0EqIdx + xwCaIdx] = 0;
        values[conti0EqIdx + xwSuspendedBiomassIdx] = 0;
        values[conti0EqIdx + xwBiosubIdx] = -waterFlux * injSub_ /FluidSystem::molarMass(xwBiosubIdx);
        values[conti0EqIdx + xwO2Idx] = -waterFlux * injO2_ /FluidSystem::molarMass(O2Idx);
        values[conti0EqIdx + xwUreaIdx] = 0;
        values[conti0EqIdx + phiCalciteIdx] = 0;
        values[conti0EqIdx + phiBiofilmIdx] = 0;
        values[conti0EqIdx + xwNaIdx] = -waterFlux * (injNa_ + injNaCorr_) /FluidSystem::molarMass(NaIdx);
        values[conti0EqIdx + xwClIdx] = -waterFlux *injTNH_ /NH3::molarMass()     //NH4Cl --->  mol Cl = mol NH4
                                        -waterFlux *injNa_ /FluidSystem::molarMass(NaIdx);      //NaCl ---> mol Cl = mol Na
        // rinse, used as standard injection fluid
        if (injProcess == InjectionProcess::rinse)
        {
            return values; // do not change anything.
        }

        // injProcess == 1 codes for an injection of mineralization medium containing urea and calcium chloride.
        // Thus, we add BC terms for those components.
        // Additionally, we need to adjust the amount of water injected due to the high concentration of other components injected.
        // Finally, we need to adapt the injected NaCorr concentration to account fo the lower pH.
        else if (injProcess == InjectionProcess::mineralization)
        {
            values[conti0EqIdx + wCompIdx] = - waterFlux * 0.8716 * densityW_ /FluidSystem::molarMass(wCompIdx);    //0.8716 factor accounts for less water per volume due to relatively high solute concentrations!
            values[conti0EqIdx + nCompIdx] = - waterFlux * injTC_ * densityW_ /FluidSystem::molarMass(nCompIdx);
            values[conti0EqIdx + xwCaIdx] = - waterFlux * injCa_/FluidSystem::molarMass(CaIdx);
            values[conti0EqIdx + xwUreaIdx] = - waterFlux * injUrea_ /FluidSystem::molarMass(UreaIdx);
            values[conti0EqIdx + xwNaIdx] = - waterFlux * injNa_ /FluidSystem::molarMass(NaIdx)
                                            - waterFlux * injNaCorr_ /FluidSystem::molarMass(NaIdx)* 0.032;
            values[conti0EqIdx + xwClIdx] = - waterFlux * injTNH_ /NH3::molarMass()               //NH4Cl --->  mol Cl = mol NH4
                                            - waterFlux * 2 * injCa_/FluidSystem::molarMass(CaIdx)             //+CaCl2 --->  mol Cl = mol Ca*2
                                            -waterFlux *injNa_ /FluidSystem::molarMass(NaIdx);      //NaCl ---> mol Cl = mol Na
            return values;
        }

        // injProcess == 3 codes for a resuscitation injection to regrow biomass.
        // It is similar to a rinse injection, but with added urea, which is what we add to the basic rinse
        else if (injProcess == InjectionProcess::resuscitation )
        {
            values[conti0EqIdx + xwUreaIdx] = - waterFlux * injUrea_ /FluidSystem::molarMass(UreaIdx);
            return values;
        }

        // injProcess == 2 codes for a inoculation or biomass injection.
        // It is similar to a rinse injection, but with added suspended biomass, which is what we add to the basic rinse
        else if (injProcess == InjectionProcess::inoculation)
        {
            values[conti0EqIdx + xwSuspendedBiomassIdx] = -waterFlux * injBiosusp_ /FluidSystem::molarMass(xwSuspendedBiomassIdx);
            return values;
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid injection process " << static_cast<int>(injProcess));
    }
    // [[/codeblock]]

    // #### Initial conditions
    // We specify the initial conditions for the primary variable depending
    // on the location. Here, we set zero model fractions everywhere in the domain except for a strip
    // at the bottom of the domain where we set an initial mole fraction of $`1e-9`$.
    // [[codeblock]]
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return initial_(globalPos); }

    // [[/codeblock]]

    // #### Reactive source and sink terms
    // We calculate the reactive source and sink terms.
    // For the details, see the "simplified chemistry case" in the dissertation of Hommel available at https://elib.uni-stuttgart.de/handle/11682/8787
    // [[codeblock]]
    NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume &scv) const
    {
        const auto elemSol = elementSolution(element, elemVolVars, fvGeometry);
        const auto gradPw = evalGradients(element,
                                          element.geometry(),
                                          this->gridGeometry(),
                                          elemSol,
                                          scv.center())[pressureIdx];
        const Scalar scvPotGradNorm = gradPw.two_norm();

        //define and compute some parameters for siplicity:
        const Scalar porosity = elemVolVars[scv].porosity();
        Scalar initialPorosity = 1.0;
        constexpr int numActiveComps = SolidSystem::numComponents-SolidSystem::numInertComponents;
        for (int i = numActiveComps; i<SolidSystem::numComponents; ++i)
            initialPorosity -= elemVolVars[scv].solidVolumeFraction(i);

        const Scalar Sw = elemVolVars[scv].saturation(wPhaseIdx);
        const Scalar xlSalinity = elemVolVars[scv].moleFraction(wPhaseIdx,NaIdx)
                                  + elemVolVars[scv].moleFraction(wPhaseIdx,CaIdx)
                                  + elemVolVars[scv].moleFraction(wPhaseIdx,ClIdx);
        const Scalar densityBiofilm = elemVolVars[scv].solidComponentDensity(SolidSystem::BiofilmIdx);
        const Scalar densityCalcite = elemVolVars[scv].solidComponentDensity(SolidSystem::CalciteIdx);
        const Scalar cBio = std::max(0.0, elemVolVars[scv].moleFraction(wPhaseIdx, SuspendedBiomassIdx)
                                          * elemVolVars[scv].molarDensity(wPhaseIdx)
                                          * FluidSystem::molarMass(SuspendedBiomassIdx)); //[kg_suspended_Biomass/m³_waterphase]
        const Scalar volFracCalcite = std::max(0.0, elemVolVars[scv].solidVolumeFraction(SolidSystem::CalciteIdx));
        const Scalar volFracBiofilm = std::max(0.0, elemVolVars[scv].solidVolumeFraction(SolidSystem::BiofilmIdx));
        const Scalar massBiofilm = densityBiofilm * volFracBiofilm;
        const Scalar cSubstrate = std::max(0.0, elemVolVars[scv].moleFraction(wPhaseIdx, BiosubIdx)
                                                * elemVolVars[scv].molarDensity(wPhaseIdx)
                                                * FluidSystem::molarMass(BiosubIdx));  //[kg_substrate/m³_waterphase]
        const Scalar cO2 = std::max(0.0, elemVolVars[scv].moleFraction(wPhaseIdx, O2Idx)
                                         * elemVolVars[scv].molarDensity(wPhaseIdx)
                                         * FluidSystem::molarMass(O2Idx));                 //[kg_oxygen/m³_waterphase]
        const Scalar mUrea = std::max(0.0, moleFracToMolality(elemVolVars[scv].moleFraction(wPhaseIdx,UreaIdx),
                                                              xlSalinity,
                                                              elemVolVars[scv].moleFraction(wPhaseIdx,nCompIdx)));  //[mol_urea/kg_H2O]
        // compute rate of ureolysis:
        const Scalar vmax = kurease_;
        const Scalar Zub = kub_ *  massBiofilm;   // [kgurease/m³]
        const Scalar rurea = vmax * Zub * mUrea / (ku_ + mUrea); //[mol/m³s]

        // compute precipitation rate of calcite, no dissolution! Simplification: rprec = rurea
        // additionally regularize the precipitation rate that we do not precipitate more calcium than available.
        const Scalar maxPrecipitationRate = elemVolVars[scv].moleFraction(wPhaseIdx,CaIdx) * Sw * porosity
                                            * elemVolVars[scv].molarDensity(wPhaseIdx) / timeStepSize_;
        const Scalar rprec = std::min(rurea, maxPrecipitationRate);

        //compute biomass growth coefficient and rate
        const Scalar mue = kmue_ * cSubstrate / (ks_ + cSubstrate) * cO2 / (ke_ + cO2);// [1/s]
        const Scalar rgf = mue * massBiofilm;                //[kg/m³s]
        const Scalar rgb = mue * porosity * Sw * cBio;   //[kg/m³s]

        // compute biomass decay coefficient and rate:
        const Scalar dcf = dc0_ + (rprec * SolidSystem::molarMass(SolidSystem::CalciteIdx)
                                   /(densityCalcite * (initialPorosity - volFracCalcite)));
        const Scalar dcb = dc0_;       //[1/s]
        const Scalar rdcf = dcf * massBiofilm; //[kg/m³s]
        const Scalar rdcb = dcb * porosity * Sw * cBio;      //[kg/m³s]

        // compute attachment coefficient and rate:
        const Scalar ka = ca1_ * volFracBiofilm + ca2_;          //[1/s]
        const Scalar ra = ka * porosity * Sw * cBio;             //[kg/m³s]

        // compute detachment coefficient and rate:
        const Scalar cd2 = volFracBiofilm / (initialPorosity - volFracCalcite);      //[-]
        const Scalar kd = cd1_ * std::pow((porosity * Sw * scvPotGradNorm),0.58) + cd2 * mue;  //[1/s]
        const Scalar rd = kd * massBiofilm;                      //[kg/m³s]

        // rprec[mol/m³s]
        // rurea[mol/m³s]
        // rgb + rdcb + ra + rd [kg/m³s]
        // source[kg/m³s]
        NumEqVector source(0.0);
        source[wCompIdx] += 0;
        source[nCompIdx] += rurea - rprec;
        source[NaIdx] += 0;
        source[ClIdx] += 0;
        source[CaIdx] += - rprec;
        source[UreaIdx] += - rurea;
        source[O2Idx] += -(rgf + rgb) *f_/yield_ / FluidSystem::molarMass(O2Idx);
        source[BiosubIdx] += -(rgf + rgb) / yield_ / FluidSystem::molarMass(BiosubIdx);
        source[SuspendedBiomassIdx] += (rgb - rdcb - ra + rd) / FluidSystem::molarMass(SuspendedBiomassIdx);
        source[phiBiofilmIdx] += (rgf - rdcf + ra - rd) / SolidSystem::molarMass(SolidSystem::BiofilmIdx);
        source[phiCalciteIdx] += + rprec;

        return source;
    }
    // [[/codeblock]]

    // The permeability is added to the vtk output
    // [[codeblock]]
    // Function to return the permeability for additional vtk output
    const std::vector<Scalar>& getPermeability()
    { return permeability_; }

    // Function to update the permeability for additional vtk output
    template<class SolutionVector>
    void updateVtkOutput(const SolutionVector& curSol)
    {
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            const auto elemSol = elementSolution(element, curSol, this->gridGeometry());
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                VolumeVariables volVars;
                volVars.update(elemSol, *this, element, scv);
                const auto dofIdxGlobal = scv.dofIndex();
                permeability_[dofIdxGlobal] = volVars.permeability();
            }
        }
    }
    // [[/codeblock]]

// ### Declaring all necessary variables and private fuctions
// The internal methods are defined here
// [[codeblock]]
private:
    // Internal method for the initial condition reused for the dirichlet conditions.
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars.setState(wPhaseOnly);
        priVars[pressureIdx] = initPressure_ ; //70e5; // - (maxHeight - globalPos[1])*densityW_*9.81; //p_atm + rho*g*h
        priVars[switchIdx] = initxwTC_;
        priVars[xwNaIdx] = initxwNa_ + xwNaCorr_;
        priVars[xwClIdx] = initxwCl_ + initxwTNH_ + 2*initxwCa_ + xwClCorr_;
        priVars[xwCaIdx] = initxwCa_;
        priVars[xwUreaIdx] = initxwUrea_;
        priVars[xwO2Idx] = initxwO2_;
        priVars[xwBiosubIdx] = initxwBiosub_;
        priVars[xwSuspendedBiomassIdx] = initxwBiosusp_;
        priVars[phiBiofilmIdx] = initBiofilm_; // [m^3/m^3]
        priVars[phiCalciteIdx] = initCalcite_; // [m^3/m^3]
        return priVars;
    }

    // Internal method to calculate the molality of a component based on its mole fraction.
    static Scalar moleFracToMolality(Scalar moleFracX, Scalar moleFracSalinity, Scalar moleFracCTot)
    {
        const Scalar molalityX = moleFracX / (1 - moleFracSalinity - moleFracCTot)
                                / FluidSystem::molarMass(FluidSystem::H2OIdx);
        return molalityX;
    }
    // [[/codeblock]]

    // The remainder of the class contains an epsilon value used for floating point comparisons
    // and parameters needed to describe the chemical processess.
    // Additionally the problem name, the peremability vector as well as some time-parameters are declared
    // [[details]] private members
    // eps is used as a small value for the definition of the boundary conditions
    static constexpr Scalar eps_ = 1e-6;

    // initial condition parameters
    Scalar initPressure_;
    Scalar densityW_;//1087; // rhow=1087;
    Scalar initxwTC_;//2.3864e-7;       // [mol/mol]
    Scalar initxwNa_;//0;
    Scalar initxwCl_;//0;
    Scalar initxwCa_;//0;
    Scalar initxwUrea_;//0;
    Scalar initxwTNH_;//3.341641e-3;
    Scalar initxwO2_;//4.4686e-6;
    Scalar initxwBiosub_;//2.97638e-4;
    Scalar initxwBiosusp_;//0;
    Scalar xwNaCorr_;//2.9466e-6;
    Scalar xwClCorr_;//0;
    Scalar initBiofilm_;
    Scalar initCalcite_;
    Scalar temperature_;

    // biomass parameters for source/sink calculations
    Scalar ca1_;
    Scalar ca2_;
    Scalar cd1_;
    Scalar dc0_;
    Scalar kmue_ ;
    Scalar f_;
    Scalar ke_;
    Scalar ks_;
    Scalar yield_;
    // urease parameters for source/sink calculations
    Scalar kub_;
    Scalar kurease_;
    Scalar ku_;

    // injection parameters
    Scalar injQ_;
    Scalar injTC_;   // [kg/kg]
    Scalar injNa_;   // [kg/m³]
    Scalar injCa_;   // [kg/m³]      //computed from CaCl2
    Scalar injUrea_; // [kg/m³]
    Scalar injTNH_;  // [kg/m³]      //computed from NH4Cl
    Scalar injO2_;   // [kg/m³]
    Scalar injSub_;  // [kg/m³]
    Scalar injBiosusp_;  // [kg/m³]
    Scalar injNaCorr_;   // [kg/m³]
    int numInjections_;
    std::vector<InjectionProcess> injectionType_;

    // the problem name
    std::string name_;
    // the permeability for output
    std::vector<Scalar> permeability_;

    // timing parameters
    Scalar time_ = 0.0;
    Scalar timeStepSize_ = 0.0;
    int episodeIdx_ = 0;
};
} // end namespace Dumux
// [[/details]]
// [[/content]]
#endif
