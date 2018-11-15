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
 * \ingroup MPNCTests
 *
 * \brief Problem showcasing the capabilities of the kinetic model.
 *
 *        The whole domain is porous medium, but the upper half has properties mimicing the ones of a free-flow domain.
 *        This way a poor man's coupling approach is accomplished: Without the complications of coupling,
 *        the main characteristics a porous and a free-flow domain are depicted.
 *
 *        The porous domain is bypassed with dry air. This way the equilibration process on top of the porous domain can be studied.
 *
 * \author Philipp Nuske
 */
#ifndef DUMUX_EVAPORATION_ATMOSPHERE_PROBLEM_HH
#define DUMUX_EVAPORATION_ATMOSPHERE_PROBLEM_HH

// this sets that the relation using pc_max is used.
// i.e. - only parameters for awn, ans are given,
//      - the fit for ans involves the maximum value for pc, where Sw, awn are zero.
// setting it here, because it impacts volume variables and spatialparameters
#define USE_PCMAX 1

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/box/properties.hh>

#include <dumux/porousmediumflow/mpnc/model.hh>
#include <dumux/porousmediumflow/mpnc/pressureformulation.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/material/fluidsystems/h2on2kinetic.hh>
#include <dumux/io/gnuplotinterface.hh>
#include "plotoverline2d.hh"
#include <dumux/material/components/constant.hh>

#include "spatialparams.hh"

// material laws for interfacial area
#include <dumux/material/fluidmatrixinteractions/2pia/efftoabslawia.hh>
#include <dumux/material/fluidmatrixinteractions/2pia/awnsurfacepolynomial2ndorder.hh>
#include <dumux/material/fluidmatrixinteractions/2pia/awnsurfacepolynomialedgezero2ndorder.hh>
#include <dumux/material/fluidmatrixinteractions/2pia/awnsurfaceexpfct.hh>
#include <dumux/material/fluidmatrixinteractions/2pia/awnsurfacepcmaxfct.hh>
#include <dumux/material/fluidmatrixinteractions/2pia/awnsurfaceexpswpcto3.hh>

namespace Dumux {
/*!
 * \ingroup MPNCTests
 *
 * \brief Problem showcasing the capabilities of the kinetic model.
 */
template <class TypeTag>
class EvaporationAtmosphereProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct EvaporationAtmosphere { using InheritsFrom = std::tuple<MPNCNonequil>; };
struct EvaporationAtmosphereBox { using InheritsFrom = std::tuple<EvaporationAtmosphere, BoxModel>; };
} // end namespace TTag

// Set the grid type
SET_TYPE_PROP(EvaporationAtmosphere, Grid, Dune::YaspGrid<2, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);

// Set the problem property
SET_TYPE_PROP(EvaporationAtmosphere, Problem, EvaporationAtmosphereProblem<TypeTag>);

// Set fluid configuration
SET_TYPE_PROP(EvaporationAtmosphere,
              FluidSystem,
              FluidSystems::H2ON2Kinetic<typename GET_PROP_TYPE(TypeTag, Scalar), FluidSystems::H2ON2DefaultPolicy</*fastButSimplifiedRelations=*/true>>);

//! Set the default pressure formulation: either pw first or pn first
SET_PROP(EvaporationAtmosphere, PressureFormulation)
{
public:
    static const MpNcPressureFormulation value = MpNcPressureFormulation::leastWettingFirst;
};

// Set the type used for scalar values
SET_TYPE_PROP(EvaporationAtmosphere, Scalar, double);

// Set the fluid system
SET_PROP(EvaporationAtmosphere, SolidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using InertComponent = Components::Constant<1, Scalar>;
    using type = SolidSystems::InertSolidPhase<Scalar, InertComponent>;
};

// Set the spatial parameters
SET_TYPE_PROP(EvaporationAtmosphere, SpatialParams, EvaporationAtmosphereSpatialParams<TypeTag>);

// Set the interfacial area relation: wetting -- non-wetting
SET_PROP(EvaporationAtmosphere, AwnSurface)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, SpatialParams)::MaterialLaw;
    using MaterialLawParams = typename MaterialLaw::Params;
    using EffectiveIALaw = AwnSurfacePcMaxFct<Scalar>;
public:
    using type = EffToAbsLawIA<EffectiveIALaw, MaterialLawParams>;
};


// Set the interfacial area relation: wetting -- solid
SET_PROP(EvaporationAtmosphere, AwsSurface)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, SpatialParams)::MaterialLaw;
    using MaterialLawParams = typename MaterialLaw::Params;
    using EffectiveIALaw = AwnSurfacePolynomial2ndOrder<Scalar>;
public:
    using type = EffToAbsLawIA<EffectiveIALaw, MaterialLawParams>;
};

// Set the interfacial area relation: non-wetting -- solid
SET_PROP(EvaporationAtmosphere, AnsSurface)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, SpatialParams)::MaterialLaw;
    using MaterialLawParams = typename MaterialLaw::Params;
    using EffectiveIALaw = AwnSurfaceExpSwPcTo3<Scalar>;
public:
    using type = EffToAbsLawIA<EffectiveIALaw, MaterialLawParams>;
};
} // end namespace Properties

/*!
 * \ingroup MpNcBoxproblems
 *
 * \brief Problem that simulates the coupled heat and mass transfer processes resulting form the evaporation of liquid water from
 *        a porous medium sub-domain into a gas filled "quasi-freeflow" sub-domain.
 */
template <class TypeTag>
class EvaporationAtmosphereProblem: public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);
    using ParameterCache = typename FluidSystem::ParameterCache;
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);

    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using Indices = typename ModelTraits::Indices;

    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases       = ModelTraits::numPhases() };
    enum { numComponents   = ModelTraits::numComponents() };
    enum { s0Idx = Indices::s0Idx };
    enum { p0Idx = Indices::p0Idx };
    enum { conti00EqIdx    = Indices::conti0EqIdx };
    enum { energyEq0Idx    = Indices::energyEqIdx };
    enum { liquidPhaseIdx       = FluidSystem::liquidPhaseIdx };
    enum { gasPhaseIdx       = FluidSystem::gasPhaseIdx };
    enum { wCompIdx        = FluidSystem::H2OIdx };
    enum { nCompIdx        = FluidSystem::N2Idx };
    enum { numEnergyEqFluid = ModelTraits::numEnergyEqFluid() };
    enum { numEnergyEqSolid = ModelTraits::numEnergyEqSolid() };

    static constexpr bool enableChemicalNonEquilibrium = GET_PROP_VALUE(TypeTag, EnableChemicalNonEquilibrium);
    using ConstraintSolver = MiscibleMultiPhaseComposition<Scalar, FluidSystem>;

    // formulations
    static constexpr auto pressureFormulation = ModelTraits::pressureFormulation();
    static constexpr auto mostWettingFirst = MpNcPressureFormulation::mostWettingFirst;
    static constexpr auto leastWettingFirst = MpNcPressureFormulation::leastWettingFirst;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    EvaporationAtmosphereProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
        : ParentType(fvGridGeometry)
    {
        percentOfEquil_         = getParam<Scalar>("BoundaryConditions.percentOfEquil");
        nTemperature_           = getParam<Scalar>("FluidSystem.nTemperature");
        nPressure_              = getParam<Scalar>("FluidSystem.nPressure");
        outputName_             = getParam<std::string>("Problem.Name");
        TInitial_               = getParam<Scalar>("InitialConditions.TInitial");
        SwPMInitial_            = getParam<Scalar>("InitialConditions.SwPMInitial");
        SwFFInitial_            = getParam<Scalar>("InitialConditions.SwFFInitial");
        pnInitial_              = getParam<Scalar>("InitialConditions.pnInitial");
        pnInjection_            = getParam<Scalar>("InitialConditions.pnInjection");
        TInject_                = getParam<Scalar>("BoundaryConditions.TInject");
        massFluxInjectedPhase_  = getParam<Scalar>("BoundaryConditions.massFluxInjectedPhase");

        // initialize the tables of the fluid system
        FluidSystem::init(TInitial_ - 15.0, 453.15, nTemperature_, // T_min, T_max, n_T
                          0.75*pnInitial_, 2.25*pnInitial_, nPressure_); // p_min, p_max, n_p
    }

    void setGridVariables(std::shared_ptr<GridVariables> gridVariables)
    { gridVariables_ = gridVariables; }

     const GridVariables& gridVariables() const
    { return *gridVariables_; }

   /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string name() const
    { return outputName_ ; }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param bTypes The boundarytypes types for the conservation equations
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        // Default: Neumann
        bcTypes.setAllNeumann();

        // Put a dirichlet somewhere: we need this for convergence
        if(onRightBoundary_(globalPos) && globalPos[1] > this->spatialParams().heightPM() + eps_)
        {
            bcTypes.setAllDirichlet();
        }
        return bcTypes;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param priVars Stores the Dirichlet values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} ] \f$
     * \param globalPos The global position
     *
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
       return initial_(globalPos);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann
     *        boundary segment in dependency on the current solution.
     *
     * \param priVars Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param intersection The intersection between element and boundary
     * \param scvIdx The local index of the sub-control volume
     * \param boundaryFaceIdx The index of the boundary face
     *
     * This method is used for cases, when the Neumann condition depends on the
     * solution and requires some quantities that are specific to the fully-implicit method.
     * The \a values store the mass flux of each phase normal to the boundary.
     * Negative values indicate an inflow.
     */
    NumEqVector neumann(const Element &element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);
        const auto& globalPos = fvGeometry.scv(scvf.insideScvIdx()).dofPosition();
        const Scalar massFluxInjectedPhase = massFluxInjectedPhase_ ;

        ParameterCache dummyCache;
        FluidState fluidState;

        for (int phaseIdx=0; phaseIdx<numPhases; phaseIdx++)
        {
            fluidState.setPressure(phaseIdx, pnInitial_);
        }

        fluidState.setTemperature(gasPhaseIdx, TInject_ );
        fluidState.setTemperature(liquidPhaseIdx, TInitial_ ); // this value is a good one, TInject does not work

        // This solves the system of equations defining x=x(p,T)
        ConstraintSolver::solve(fluidState,
                                dummyCache) ;

        // Now let's make the air phase less than fully saturated with water
        fluidState.setMoleFraction(gasPhaseIdx, wCompIdx, fluidState.moleFraction(gasPhaseIdx, wCompIdx)*percentOfEquil_ ) ;
        fluidState.setMoleFraction(gasPhaseIdx, nCompIdx, 1.-fluidState.moleFraction(gasPhaseIdx, wCompIdx) ) ;

        // compute density of injection phase
        const Scalar density = FluidSystem::density(fluidState,
                                                    dummyCache,
                                                    gasPhaseIdx);
        fluidState.setDensity(gasPhaseIdx, density);
        const Scalar molarDensity = FluidSystem::molarDensity(fluidState,
                                                              dummyCache,
                                                              gasPhaseIdx);
        fluidState.setMolarDensity(gasPhaseIdx, molarDensity);

        for(int phaseIdx=0; phaseIdx<numPhases; phaseIdx++)
        {
                const Scalar h = FluidSystem::enthalpy(fluidState,
                                                       dummyCache,
                                                       phaseIdx);
                fluidState.setEnthalpy(phaseIdx, h);
        }

        const Scalar molarFlux = massFluxInjectedPhase / fluidState.averageMolarMass(gasPhaseIdx);

        // actually setting the fluxes
        if (onLeftBoundary_(globalPos) && this->spatialParams().inFF_(globalPos))
        {
            values[conti00EqIdx + gasPhaseIdx * numComponents + wCompIdx]
             = -molarFlux * fluidState.moleFraction(gasPhaseIdx, wCompIdx);
            values[conti00EqIdx + gasPhaseIdx * numComponents + nCompIdx]
             = -molarFlux * fluidState.moleFraction(gasPhaseIdx, nCompIdx);
            // energy equations are specified mass specifically
            values[energyEq0Idx + gasPhaseIdx] = - massFluxInjectedPhase
                                                    * fluidState.enthalpy(gasPhaseIdx) ;
        }
        return values;
    }

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param priVars Stores the initial solution for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        return initial_(globalPos);
    }

    /*!
     * \brief Evaluate the source term for all balance equations within a given
     *        sub-control-volume.
     *
     * \param priVars Stores the solution for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} / (m^\textrm{dim} \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     *
     *      Positive values mean that mass is created, negative ones mean that it vanishes.
     */
    NumEqVector source(const Element &element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume &scv) const
    {
        NumEqVector values(0.0);
        return values;

    }

private:
    // the internal method for the initial condition
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        const Scalar T = TInitial_;
        Scalar S[numPhases];

        if (this->spatialParams().inPM_(globalPos)){
            S[liquidPhaseIdx]    = SwPMInitial_;
            S[gasPhaseIdx]    = 1. - S[liquidPhaseIdx] ;
        }
        else if (this->spatialParams().inFF_(globalPos)){
            S[liquidPhaseIdx]    = SwFFInitial_;
            S[gasPhaseIdx]    = 1. - S[liquidPhaseIdx] ;
        }
        else
            DUNE_THROW(Dune::InvalidStateException,
                       "You should not be here: x=" << globalPos[0] << " y= "<< globalPos[dimWorld-1]);

        for (int i = 0; i < numPhases - 1; ++i)
            priVars[s0Idx + i] = S[i];

        // capillary pressure Params
        FluidState equilibriumFluidState;

        //set saturation to inital values, this needs to be done in order for the fluidState to tell me pc
        for (int phaseIdx = 0; phaseIdx < numPhases ; ++phaseIdx)
        {
            equilibriumFluidState.setSaturation(phaseIdx, S[phaseIdx]);
            equilibriumFluidState.setTemperature(phaseIdx, TInitial_ );
        }

        const auto &materialParams =
            this->spatialParams().materialLawParamsAtPos(globalPos);
        Scalar capPress[numPhases];

        //obtain pc according to saturation
        using MaterialLaw = typename ParentType::SpatialParams::MaterialLaw;
        MaterialLaw::capillaryPressures(capPress, materialParams, equilibriumFluidState);

        Scalar p[numPhases];
        if (this->spatialParams().inPM_(globalPos)){
            // Use homogenous pressure in the domain and let the newton find the pressure distribution
            using std::abs;
            p[liquidPhaseIdx] = pnInitial_  - abs(capPress[liquidPhaseIdx]);
            p[gasPhaseIdx] = p[liquidPhaseIdx] + abs(capPress[liquidPhaseIdx]);
        }
        else if (this->spatialParams().inFF_(globalPos)){
            p[gasPhaseIdx] = pnInitial_ ;
            p[liquidPhaseIdx] = pnInitial_  ;
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "You should not be here: x=" << globalPos[0] << " y= "<< globalPos[dimWorld-1]);

        if(pressureFormulation == mostWettingFirst){
            // This means that the pressures are sorted from the most wetting to the least wetting-1 in the primary variables vector.
            // For two phases this means that there is one pressure as primary variable: pw
            priVars[p0Idx] = p[liquidPhaseIdx];
        }
        else if(pressureFormulation == leastWettingFirst){
            // This means that the pressures are sorted from the least wetting to the most wetting-1 in the primary variables vector.
            // For two phases this means that there is one pressure as primary variable: pn
            priVars[p0Idx] = p[gasPhaseIdx];
        }
        else DUNE_THROW(Dune::InvalidStateException, "EvaporationAtmosphereProblem does not support the chosen pressure formulation.");

        for (int energyEqIdx=0; energyEqIdx< numEnergyEqFluid+numEnergyEqSolid; ++energyEqIdx)
                priVars[energyEq0Idx + energyEqIdx] = T;

        for (int phaseIdx=0; phaseIdx<numPhases; phaseIdx++)
             equilibriumFluidState.setPressure(phaseIdx, p[phaseIdx]);

         // This solves the system of equations defining x=x(p,T)
        ParameterCache dummyCache;
        ConstraintSolver::solve(equilibriumFluidState,
                                dummyCache) ;

        FluidState dryFluidState(equilibriumFluidState);
        // Now let's make the air phase less than fully saturated with vapor
        dryFluidState.setMoleFraction(gasPhaseIdx, wCompIdx, dryFluidState.moleFraction(gasPhaseIdx, wCompIdx) * percentOfEquil_ ) ;
        dryFluidState.setMoleFraction(gasPhaseIdx, nCompIdx, 1.0-dryFluidState.moleFraction(gasPhaseIdx, wCompIdx) ) ;

        /* Difference between kinetic and MPNC:
         * number of component related primVar and how they are calculated (mole fraction, fugacities, resp.)
         */
        if(enableChemicalNonEquilibrium){
            for(int phaseIdx=0; phaseIdx < numPhases; ++ phaseIdx)
            {
                for(int compIdx=0; compIdx <numComponents; ++compIdx){
                    int offset = compIdx + phaseIdx * numComponents  ;

                    if (this->spatialParams().inPM_(globalPos)){
                        priVars[conti00EqIdx + offset] = equilibriumFluidState.moleFraction(phaseIdx,compIdx) ;
                    }
                    else if (this->spatialParams().inFF_(globalPos)){
                        priVars[conti00EqIdx + offset] = dryFluidState.moleFraction(phaseIdx,compIdx) ;
                    }
                    else
                        DUNE_THROW(Dune::InvalidStateException, "You should not be here: x=" << globalPos[0] << " y= "<< globalPos[dimWorld-1]);
                }
            }
        }
        else
        {
            // in the case I am using the "standard" mpnc model, the variables to be set are the "fugacities"
            const Scalar fugH2O = FluidSystem::H2O::vaporPressure(T) ;
            const Scalar fugN2 = p[gasPhaseIdx] - fugH2O ;

            priVars[conti00EqIdx + FluidSystem::N2Idx] = fugN2 ;
            priVars[conti00EqIdx + FluidSystem::H2OIdx] = fugH2O ;

            Scalar xl[numComponents];
            Scalar beta[numComponents];

            const Scalar Henry              = BinaryCoeff::H2O_N2::henry(TInitial_);
            const Scalar satVapPressure     = FluidSystem::H2O::vaporPressure(TInitial_);
            xl[FluidSystem::H2OIdx]         = x_[liquidPhaseIdx][wCompIdx];
            xl[FluidSystem::N2Idx]          = x_[liquidPhaseIdx][nCompIdx];
            beta[FluidSystem::H2OIdx]       = satVapPressure ;
            beta[FluidSystem::N2Idx]        = Henry ;

            for (int i = 0; i < numComponents; ++i)
                priVars[conti00EqIdx + i] = xl[i]*beta[i]; // this should be really fug0Idx but the compiler only knows one or the other
        }
        return priVars;
    }

    /*!
     * \brief Give back whether the tested position (input) is a specific region (left) in the domain
     */
    bool onLeftBoundary_(const GlobalPosition & globalPos) const
    {  return globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + eps_;   }

    /*!
     * \brief Give back whether the tested position (input) is a specific region (right) in the domain
     */
    bool onRightBoundary_(const GlobalPosition & globalPos) const
    { return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_;    }

    /*!
     * \brief Give back whether the tested position (input) is a specific region (down, (gravityDir)) in the domain
     */
    bool onLowerBoundary_(const GlobalPosition & globalPos) const
    { return globalPos[dimWorld-1] < this->fvGridGeometry().bBoxMin()[dimWorld-1] + eps_;    }

private:
    static constexpr Scalar eps_ = 1e-6;
    Scalar percentOfEquil_ ;
    int nTemperature_;
    int nPressure_;
    std::string outputName_;
    Scalar heatIntoSolid_;
    Scalar TInitial_ ;
    Scalar SwPMInitial_ ;
    Scalar SwFFInitial_ ;
    Scalar SnInitial_;
    Scalar pnInitial_;
    Scalar pnInjection_;
    Dune::ParameterTree inputParameters_;
    Scalar x_[numPhases][numComponents] ;

    Scalar TInject_;

    Scalar massFluxInjectedPhase_ ;

    std::shared_ptr<GridVariables> gridVariables_;

public:

    Dune::ParameterTree getInputParameters() const
    { return inputParameters_; }
};

} //end namespace

#endif
