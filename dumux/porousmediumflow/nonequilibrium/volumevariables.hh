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
 * \ingroup PorousmediumNonEquilibriumModel
 * \brief This class contains the volume variables required for the
 *        modules which require the specific interfacial area between
 *        fluid phases.
 *
 * This files contains all specializations which use 'real'
 * interfacial areas.
 */
#ifndef DUMUX_NONEQUILIBRIUM_VOLUME_VARIABLES__HH
#define DUMUX_NONEQUILIBRIUM_VOLUME_VARIABLES__HH

#include <dumux/common/dimensionlessnumbers.hh>

#include <dumux/porousmediumflow/volumevariables.hh>

#include <dumux/material/fluidstates/nonequilibrium.hh>
#include <dumux/material/constraintsolvers/misciblemultiphasecomposition.hh>

namespace Dumux {

/*!
 * \file
 * \ingroup PorousmediumNonEquilibriumModel
 * \brief This class contains the volume variables required for the
 *        modules which require the specific interfacial area between
 *        fluid phases.
 */
// forward declaration
template <class TypeTag, bool enableChemicalNonEquilibrium ,bool enableThermalNonEquilibrium>
class NonEquilibriumVolumeVariablesImplementation;

template <class TypeTag>
using NonEquilibriumVolumeVariables =
NonEquilibriumVolumeVariablesImplementation<TypeTag, GET_PROP_VALUE(TypeTag, EnableChemicalNonEquilibrium), GET_PROP_VALUE(TypeTag, EnableThermalNonEquilibrium)>;

template <class TypeTag>
class NonEquilibriumVolumeVariablesImplementation<TypeTag, false/*enableChemicalNonEquilibrium*/, false/*enableThermalNonEquilibrium*/>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);
    using ParameterCache = typename FluidSystem::ParameterCache;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
public:
    template<class ElementSolution>
    void updateInterfacialArea(const ElementSolution& elemSol,
                               const FluidState & fluidState,
                               const ParameterCache &paramCache,
                               const Problem &problem,
                               const Element & element,
                               const SubControlVolume& scv)
    {}

    template<class ElementSolution>
    void updateTemperatures(const ElementSolution& elemSol,
                             const Problem &problem,
                             const Element& element,
                             const SubControlVolume& scv,
                             FluidState& fluidState)
    {}


    void updateMoleFraction(FluidState & actualFluidState,
                            ParameterCache & paramCache,
                            const PrimaryVariables & priVars)
    {}

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// // specialization for the case of NO kinetic mass but kinetic energy transfer of a fluid mixture and solid
// /////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <class TypeTag>
class NonEquilibriumVolumeVariablesImplementation<TypeTag, /*enableChemicalNonEquilibrium*/false, /*enableThermalNonEquilibrium*/true>
{
    using BaseType = PorousMediumFlowVolumeVariables<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);
    using ParameterCache = typename FluidSystem::ParameterCache;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { nusseltFormulation = GET_PROP_VALUE(TypeTag, NusseltFormulation)} ;
    enum { numEnergyEqFluid = GET_PROP_VALUE(TypeTag, NumEnergyEqFluid) };
    enum { numEnergyEqSolid = GET_PROP_VALUE(TypeTag, NumEnergyEqSolid) };
    enum { temperature0Idx = Indices::temperature0Idx };

    static_assert((numEnergyEqFluid < 2),
                   "This model is a specialization for a energy transfer of a fluid  mixture and a solid");

    using DimLessNum = DimensionlessNumbers<Scalar>;
public:
    /*!
     * \brief Updates the volume specific interfacial area [m^2 / m^3] between the phases.
     */
    template<class ElementSolution>
    void updateInterfacialArea(const ElementSolution& elemSol,
                               const FluidState & fluidState,
                               const ParameterCache &paramCache,
                               const Problem &problem,
                               const Element & element,
                               const SubControlVolume& scv)
    {
        factorMassTransfer_   = problem.spatialParams().factorMassTransfer(element, scv, elemSol);
        factorEnergyTransfer_   = problem.spatialParams().factorEnergyTransfer(element, scv, elemSol);
        characteristicLength_   = problem.spatialParams().characteristicLength(element, scv, elemSol);

        // setting the dimensionless numbers.
        // obtaining the respective quantities.
         const unsigned int vIdxGlobal = scv.dofIndex();
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            const Scalar darcyMagVelocity     = problem.gridVariables().volumeDarcyMagVelocity(phaseIdx, vIdxGlobal);
            const Scalar dynamicViscosity     = fluidState.viscosity(phaseIdx);
            const Scalar density              = fluidState.density(phaseIdx);
            const Scalar kinematicViscosity   = dynamicViscosity / density;
            const Scalar heatCapacity         = FluidSystem::heatCapacity(fluidState,
                                                                          paramCache,
                                                                          phaseIdx);
            const Scalar thermalConductivity  = FluidSystem::thermalConductivity(fluidState,
                                                                           paramCache,
                                                                           phaseIdx);
            const Scalar porosity = problem.spatialParams().porosity(element,
                                                                   scv,
                                                                   elemSol);

            reynoldsNumber_[phaseIdx]   = DimLessNum::reynoldsNumber(darcyMagVelocity,
                                                                     characteristicLength_,
                                                                     kinematicViscosity);

            prandtlNumber_[phaseIdx]    = DimLessNum::prandtlNumber(dynamicViscosity,
                                                                    heatCapacity,
                                                                    thermalConductivity);


            nusseltNumber_[phaseIdx]    = DimLessNum::nusseltNumberForced(reynoldsNumber_[phaseIdx],
                                                                          prandtlNumber_[phaseIdx],
                                                                          porosity,
                                                                          nusseltFormulation);
        }
    }

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElementSolution>
    void updateTemperatures(const ElementSolution& elemSol,
                            const Problem &problem,
                            const Element& element,
                            const SubControlVolume& scv,
                            FluidState& fluidState)
    {
        if (numEnergyEqFluid >1)
        for(int phaseIdx=0; phaseIdx < numEnergyEqFluid; ++phaseIdx)
        {
            // retrieve temperature from solution vector
            const Scalar T = BaseType::extractDofPriVars(elemSol, scv)[temperature0Idx + phaseIdx];
            fluidState.setTemperature(phaseIdx, T);
        }
        else
        {
            const Scalar T = BaseType::extractDofPriVars(elemSol, scv)[temperature0Idx];
            fluidState.setTemperature(T);
        }
        for(int solidPhaseIdx = numEnergyEqFluid; solidPhaseIdx < numEnergyEqFluid+numEnergyEqSolid; ++solidPhaseIdx)
        {
            temperatureSolid_ = BaseType::extractDofPriVars(elemSol, scv)[temperature0Idx + solidPhaseIdx];

        }
    }

    void updateMoleFraction(FluidState & actualFluidState,
                            ParameterCache & paramCache,
                            const PrimaryVariables & priVars)
    { }


    /*!
     * \brief Returns the temperature in fluid / solid phase(s)
     *        the sub-control volume.
     * \param phaseIdx The local index of the phases
     */
    Scalar temperatureSolid() const
    { return temperatureSolid_; }

    //! access function Reynolds Number
    const Scalar reynoldsNumber(const unsigned int phaseIdx) const
    { return reynoldsNumber_[phaseIdx]; }

    //! access function Prandtl Number
    const Scalar prandtlNumber(const unsigned int phaseIdx) const
    { return prandtlNumber_[phaseIdx]; }

    //! access function Nusselt Number
    const Scalar nusseltNumber(const unsigned int phaseIdx) const
    { return nusseltNumber_[phaseIdx]; }

    //! access function characteristic length
    const Scalar characteristicLength() const
    { return characteristicLength_; }

    //! access function pre factor energy transfer
    const Scalar factorEnergyTransfer() const
    { return factorEnergyTransfer_; }

    //! access function pre factor mass transfer
    const Scalar factorMassTransfer() const
    { return factorMassTransfer_; }

    /*!
     * \brief If running in valgrind this makes sure that all
     *        quantities in the volume variables are defined.
     */
    void checkDefined() const
    {
#if !defined NDEBUG && HAVE_VALGRIND
        Valgrind::CheckDefined(reynoldsNumber_);
        Valgrind::CheckDefined(prandtlNumber_);
        Valgrind::CheckDefined(nusseltNumber_);
        Valgrind::CheckDefined(characteristicLength_);
        Valgrind::CheckDefined(factorEnergyTransfer_);
        Valgrind::CheckDefined(factorMassTransfer_);
        Valgrind::CheckDefined(temperatureSolid_);
#endif
    }

private:
    //! dimensionless numbers
    Scalar reynoldsNumber_[numPhases];
    Scalar prandtlNumber_[numPhases];
    Scalar nusseltNumber_[numPhases];
    Scalar characteristicLength_;
    Scalar factorEnergyTransfer_;
    Scalar factorMassTransfer_;
    Scalar solidSurface_ ;
    Scalar temperatureSolid_;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
// // specialization for the case of (only) kinetic mass transfer. Be careful, we do not test this!
// ////////////////////////////////////////////////////////////////////////////////////////////////////
template <class TypeTag>
class NonEquilibriumVolumeVariablesImplementation<TypeTag, true/*enableChemicalNonEquilibrium*/, false/*enableThermalNonEquilibrium*/>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);
    using ParameterCache = typename FluidSystem::ParameterCache;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);

    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { wPhaseIdx = FluidSystem::wPhaseIdx };
    enum { nPhaseIdx = FluidSystem::nPhaseIdx };
    enum { nCompIdx = FluidSystem::nCompIdx } ;
    enum { wCompIdx = FluidSystem::wCompIdx } ;
    enum { dim = GridView::dimension};
    enum { sherwoodFormulation = GET_PROP_VALUE(TypeTag, SherwoodFormulation)} ;
    enum { moleFrac00Idx = Indices::conti0EqIdx };

    using AwnSurface = typename GET_PROP_TYPE(TypeTag, AwnSurface);
    using AwnSurfaceParams = typename  AwnSurface::Params;

    using DimLessNum = DimensionlessNumbers<Scalar>;
    using ConstraintSolver = MiscibleMultiPhaseComposition<Scalar, FluidSystem>;
public:
    /*!
     * \brief Updates the volume specific interfacial area [m^2 / m^3] between the phases.
     */
    template<class ElementSolution>
    void updateInterfacialArea(const ElementSolution& elemSol,
                               const FluidState & fluidState,
                               const ParameterCache &paramCache,
                               const Problem &problem,
                               const Element & element,
                               const SubControlVolume& scv)
    {
        //obtain parameters for awnsurface
        const AwnSurfaceParams & awnSurfaceParams = problem.spatialParams().aWettingNonWettingSurfaceParams(element, scv, elemSol) ;

        // obtain (standard) material parameters (needed for the residual saturations)
        const auto &materialParams  = problem.spatialParams().materialLawParams(element, scv, elemSol) ;

        Valgrind::CheckDefined(awnSurfaceParams);
        const Scalar Sw = fluidState.saturation(wPhaseIdx) ;
        const Scalar pc = fluidState.pressure(nPhaseIdx) - fluidState.pressure(wPhaseIdx);

        // so far there is only a model for kinetic mass transfer between fluid phases
        interfacialArea_ = AwnSurface::interfacialArea(awnSurfaceParams, materialParams, Sw, pc );

        Valgrind::CheckDefined(interfacialArea_);

        factorMassTransfer_   = problem.spatialParams().factorMassTransfer(element, scv, elemSol);

        characteristicLength_   = problem.spatialParams().characteristicLength(element, scv, elemSol);
        // setting the dimensionless numbers.
        // obtaining the respective quantities.
        int globalVertexIdx = problem.vertexMapper().subIndex(element, scv, dim);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            const Scalar darcyMagVelocity     = problem.gridVariables().volumeDarcyMagVelocity(phaseIdx, globalVertexIdx);
            const Scalar dynamicViscosity     = fluidState.viscosity(phaseIdx);
            const Scalar density              = fluidState.density(phaseIdx);
            const Scalar kinematicViscosity   = dynamicViscosity / density;

            // diffusion coefficient of non-wetting component in wetting phase
            const Scalar diffCoeff = FluidSystem::binaryDiffusionCoefficient(fluidState, paramCache, phaseIdx, wCompIdx, nCompIdx);

            reynoldsNumber_[phaseIdx]   = DimLessNum::reynoldsNumber(darcyMagVelocity,
                                                                     characteristicLength_,
                                                                     kinematicViscosity);

            schmidtNumber_[phaseIdx]    = DimLessNum::schmidtNumber(dynamicViscosity,
                                                                    density,
                                                                    diffCoeff);
            sherwoodNumber_[phaseIdx]   = DimLessNum::sherwoodNumber(reynoldsNumber_[phaseIdx],
                                                                     schmidtNumber_[phaseIdx],
                                                                     sherwoodFormulation);
        }
    }

     /*!
     * \brief Update composition of all phases in the mutable
     *        parameters from the primary variables.
     *
     *        \param actualFluidState Container for all the secondary variables concerning the fluids
     *        \param paramCache Container for cache parameters
     *        \param priVars The primary Variables
     *        \param *hint the volume variables, usable for initial guess of composition
     */
    void updateMoleFraction(FluidState & actualFluidState,
                            ParameterCache & paramCache,
                            const PrimaryVariables & priVars)
    {
        // setting the mole fractions of the fluid state
        for(int phaseIdx=0; phaseIdx<numPhases; ++phaseIdx)
        {
                // set the component mole fractions
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    actualFluidState.setMoleFraction(phaseIdx,
                           compIdx,
                           priVars[moleFrac00Idx +
                                   phaseIdx*numComponents +
                                   compIdx]);
                }
        }

//          // For using the ... other way of calculating equilibrium
//          THIS IS ONLY FOR silencing Valgrind but is not used in this model
            for(int phaseIdx=0; phaseIdx<numPhases; ++phaseIdx)
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    const Scalar phi = FluidSystem::fugacityCoefficient(actualFluidState,
                                                                        paramCache,
                                                                        phaseIdx,
                                                                        compIdx);
                    actualFluidState.setFugacityCoefficient(phaseIdx,
                                                      compIdx,
                                                      phi);
            }

            FluidState equilFluidState; // the fluidState *on the interface* i.e. chemical equilibrium
            equilFluidState.assign(actualFluidState) ;
            ConstraintSolver::solve(equilFluidState,
                                    paramCache,
                                    /*setViscosity=*/false,
                                    /*setEnthalpy=*/false) ;

            // Setting the equilibrium composition (in a kinetic model not necessarily the same as the actual mole fraction)
            for(int phaseIdx=0; phaseIdx<numPhases; ++phaseIdx){
                for (int compIdx=0; compIdx< numComponents; ++ compIdx){
                    xEquil_[phaseIdx][compIdx] = equilFluidState.moleFraction(phaseIdx, compIdx);
                }
            }

            // compute densities of all phases
            for(int phaseIdx=0; phaseIdx<numPhases; ++phaseIdx){
                const Scalar rho = FluidSystem::density(actualFluidState, paramCache, phaseIdx);
                actualFluidState.setDensity(phaseIdx, rho);
            }

        }

    /*!
     * \brief The mole fraction we would have in the case of chemical equilibrium /
     *        on the interface.
     *
     *     \param phaseIdx The index of the fluid phase
     *     \param compIdx The local index of the component
     */
    const Scalar xEquil(const unsigned int phaseIdx, const unsigned int compIdx) const
    {
        return xEquil_[phaseIdx][compIdx] ;
    }

    /*!
     * \brief The specific interfacial area between two fluid phases [m^2 / m^3]
     */
    const Scalar interfacialArea(const unsigned int phaseIIdx, const unsigned int phaseJIdx) const
    {
        // so far there is only a model for kinetic mass transfer between fluid phases
        assert((phaseIIdx == nPhaseIdx && phaseJIdx == wPhaseIdx)
              || (phaseIIdx == wPhaseIdx && phaseJIdx == nPhaseIdx));
        return interfacialArea_;
    }

    //! access function Reynolds Number
    const Scalar reynoldsNumber(const unsigned int phaseIdx) const
    { return reynoldsNumber_[phaseIdx]; }

    //! access function Schmidt Number
    const Scalar schmidtNumber(const unsigned int phaseIdx) const
    { return schmidtNumber_[phaseIdx]; }

    //! access function Sherwood Number
    const Scalar sherwoodNumber(const unsigned int phaseIdx) const
    { return sherwoodNumber_[phaseIdx]; }

    //! access function characteristic length
    const Scalar characteristicLength() const
    { return characteristicLength_; }

    //! access function pre factor mass transfer
    const Scalar factorMassTransfer() const
    { return factorMassTransfer_; }

    /*!
     * \brief If running in valgrind this makes sure that all
     *        quantities in the volume variables are defined.
     */
    void checkDefined() const
    {
        Valgrind::CheckDefined(interfacialArea_);
        Valgrind::CheckDefined(characteristicLength_);
        Valgrind::CheckDefined(factorMassTransfer_);
        Valgrind::CheckDefined(reynoldsNumber_);
        Valgrind::CheckDefined(schmidtNumber_);
        Valgrind::CheckDefined(sherwoodNumber_);
    }

private:
    Scalar characteristicLength_;
    Scalar factorMassTransfer_;
    Scalar solidSurface_ ;
    Scalar interfacialArea_ ;
    Scalar sherwoodNumber_[numPhases] ;
    Scalar schmidtNumber_[numPhases] ;
    Scalar reynoldsNumber_[numPhases] ;
    Scalar xEquil_[numPhases][numComponents];
};

//this is a specialization where everything is in non-equilibrium. we have to do our own stuff for the interfacial area but can use the rest from the others
template <class TypeTag>
class NonEquilibriumVolumeVariablesImplementation<TypeTag, true/*enableChemicalNonEquilibrium*/, true/*enableThermalNonEquilibrium*/>
{
    using BaseType = PorousMediumFlowVolumeVariables<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);
    using ParameterCache = typename FluidSystem::ParameterCache;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);

    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { wPhaseIdx = FluidSystem::wPhaseIdx };
    enum { nPhaseIdx = FluidSystem::nPhaseIdx };
    enum { sPhaseIdx = FluidSystem::sPhaseIdx };
    enum { nCompIdx = FluidSystem::nCompIdx } ;
    enum { wCompIdx = FluidSystem::wCompIdx } ;
    enum { nusseltFormulation = GET_PROP_VALUE(TypeTag, NusseltFormulation)} ;
    enum { sherwoodFormulation = GET_PROP_VALUE(TypeTag, SherwoodFormulation)} ;
    enum { numEnergyEqFluid = GET_PROP_VALUE(TypeTag, NumEnergyEqFluid) };
    enum { numEnergyEqSolid = GET_PROP_VALUE(TypeTag, NumEnergyEqSolid) };
    enum { temperature0Idx = Indices::temperature0Idx };
    enum { moleFrac00Idx = Indices::conti0EqIdx };

    static_assert((numEnergyEqFluid > 1),
                   "This model only deals with energy transfer between two fluids and one solid phase");
    using DimLessNum = DimensionlessNumbers<Scalar>;

    using AwnSurface = typename GET_PROP_TYPE(TypeTag, AwnSurface);
    using AwnSurfaceParams = typename  AwnSurface::Params;
    using AwsSurface = typename GET_PROP_TYPE(TypeTag, AwsSurface);
    using AwsSurfaceParams = typename  AwsSurface::Params;
    using AnsSurface = typename GET_PROP_TYPE(TypeTag, AnsSurface);
    using AnsSurfaceParams = typename  AnsSurface::Params;

    using ConstraintSolver = MiscibleMultiPhaseComposition<Scalar, FluidSystem>;


public:
    /*!
     * \brief Updates the volume specific interfacial area [m^2 / m^3] between the phases.
     */
    template<class ElementSolution>
    void updateInterfacialArea(const ElementSolution& elemSol,
                                const FluidState & fluidState,
                                const ParameterCache &paramCache,
                                const Problem &problem,
                                const Element & element,
                                const SubControlVolume& scv)
    {
        // obtain (standard) material parameters (needed for the residual saturations)
       const auto& materialParams =
            problem.spatialParams().materialLawParams(element, scv, elemSol);

        //obtain parameters for interfacial area constitutive relations
        const AwnSurfaceParams & aWettingNonWettingSurfaceParams
               = problem.spatialParams().aWettingNonWettingSurfaceParams(element, scv, elemSol);
        const AnsSurfaceParams & aNonWettingSolidSurfaceParams
               = problem.spatialParams().aNonWettingSolidSurfaceParams(element, scv, elemSol);

        Valgrind::CheckDefined(aWettingNonWettingSurfaceParams);
        Valgrind::CheckDefined(aNonWettingSolidSurfaceParams);

        const Scalar pc = fluidState.pressure(nPhaseIdx) - fluidState.pressure(wPhaseIdx);
        const Scalar Sw = fluidState.saturation(wPhaseIdx);
        Valgrind::CheckDefined(Sw);
        Valgrind::CheckDefined(pc);
        Scalar awn;

#define AwnRegul 0
        // This regularizes the interfacial area between the fluid phases.
        // This makes sure, that
        // a) some saturation cannot be lost: Never leave two phase region.
        // b) We cannot leave the fit region: no crazy (e.g. negative) values possible
//        const Scalar Swr =  aWettingNonWettingSurfaceParams.Swr() ;
//        const Scalar Snr =  aWettingNonWettingSurfaceParams.Snr() ;

        // this just leads to a stalling newton error as soon as this kicks in.
        // May be a spline or sth like this would help, but I do not which derivatives
        // to specify.
#if AwnRegul
        if(Sw < 5e-3 ) // or Sw > (1.-1e-5 )
        {
            awn = 0. ; // 10.; //
        }
        else
#endif
        awn = AwnSurface::interfacialArea(aWettingNonWettingSurfaceParams, materialParams, Sw, pc ); // 10.; //

        interfacialArea_[wPhaseIdx][nPhaseIdx] = awn ; //10. ;//
        interfacialArea_[nPhaseIdx][wPhaseIdx] = interfacialArea_[wPhaseIdx][nPhaseIdx];
        interfacialArea_[wPhaseIdx][wPhaseIdx] = 0. ;

        Scalar ans = AnsSurface::interfacialArea(aNonWettingSolidSurfaceParams, materialParams,Sw, pc ); // 10.; //
//        if (ans <0 )
//            ans = 0 ;

// Switch for using a a_{wn} relations that has some "maximum capillary pressure" as parameter.
// That value is obtained by regularization of the pc(Sw) function.
#if USE_PCMAX
       const Scalar pcMax = problem.spatialParams().pcMax(element,
                                                          scv,
                                                          elemSol);
        // I know the solid surface from the pore network. But it is more consistent to use the fit value.
        solidSurface_   = AnsSurface::interfacialArea(aNonWettingSolidSurfaceParams, materialParams, /*Sw=*/0., pcMax );
        Valgrind::CheckDefined(solidSurface_);
#endif
        interfacialArea_[nPhaseIdx][sPhaseIdx] = ans ; //10. ; //
        interfacialArea_[sPhaseIdx][nPhaseIdx] = interfacialArea_[nPhaseIdx][sPhaseIdx];
        interfacialArea_[nPhaseIdx][nPhaseIdx] = 0. ;

#if USE_PCMAX
        const Scalar aws = solidSurface_ - ans ;
        interfacialArea_[wPhaseIdx][sPhaseIdx] = aws ; //10. ; //
        interfacialArea_[sPhaseIdx][wPhaseIdx] = interfacialArea_[wPhaseIdx][sPhaseIdx];
        interfacialArea_[sPhaseIdx][sPhaseIdx] = 0. ;
#else
        const AwsSurfaceParams & aWettingSolidSurfaceParams
               = problem.spatialParams().aWettingSolidSurfaceParams();
        Valgrind::CheckDefined(aWettingSolidSurfaceParams);
        const Scalar aws = AwsSurface::interfacialArea(aWettingSolidSurfaceParams,materialParams, Sw, pc ); // 10.; //
        interfacialArea_[wPhaseIdx][sPhaseIdx] = aws ;
        interfacialArea_[sPhaseIdx][wPhaseIdx] = interfacialArea_[wPhaseIdx][sPhaseIdx];
        interfacialArea_[sPhaseIdx][sPhaseIdx] = 0. ;
#endif

        Valgrind::CheckDefined(interfacialArea_);

        factorMassTransfer_   = problem.spatialParams().factorMassTransfer(element, scv, elemSol);

        factorEnergyTransfer_   = problem.spatialParams().factorEnergyTransfer(element, scv, elemSol);

        characteristicLength_   = problem.spatialParams().characteristicLength(element, scv, elemSol);


        const unsigned int vIdxGlobal = scv.dofIndex();
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            const Scalar darcyMagVelocity     = problem.gridVariables().volumeDarcyMagVelocity(phaseIdx, vIdxGlobal);
            const Scalar dynamicViscosity     = fluidState.viscosity(phaseIdx);
            const Scalar density              = fluidState.density(phaseIdx);
            const Scalar kinematicViscosity   = dynamicViscosity / density;
            const Scalar heatCapacity         = FluidSystem::heatCapacity(fluidState, paramCache, phaseIdx);
            const Scalar thermalConductivity  = FluidSystem::thermalConductivity(fluidState, paramCache, phaseIdx);


            // diffusion coefficient of non-wetting component in wetting phase
            const Scalar diffCoeff =  FluidSystem::binaryDiffusionCoefficient(fluidState, paramCache, phaseIdx, wCompIdx, nCompIdx);
            const Scalar porosity = problem.spatialParams().porosity(element,
                                                                    scv,
                                                                   elemSol);

            reynoldsNumber_[phaseIdx]   = DimLessNum::reynoldsNumber(darcyMagVelocity,
                                                                     characteristicLength_,
                                                                     kinematicViscosity);

            prandtlNumber_[phaseIdx]    = DimLessNum::prandtlNumber(dynamicViscosity,
                                                                    heatCapacity,
                                                                    thermalConductivity);

            nusseltNumber_[phaseIdx]    = DimLessNum::nusseltNumberForced(reynoldsNumber_[phaseIdx],
                                                                          prandtlNumber_[phaseIdx],
                                                                          porosity,
                                                                          nusseltFormulation);

            schmidtNumber_[phaseIdx]    = DimLessNum::schmidtNumber(dynamicViscosity,
                                                                    density,
                                                                    diffCoeff);

            // If Diffusion is not enabled, Sherwood is divided by zero
            sherwoodNumber_[phaseIdx]   = DimLessNum::sherwoodNumber(reynoldsNumber_[phaseIdx],
                                                                      schmidtNumber_[phaseIdx],
                                                                      sherwoodFormulation);
        }
    }

    void updateMoleFraction(FluidState & actualFluidState,
                            ParameterCache & paramCache,
                            const PrimaryVariables & priVars)
    {
        // setting the mole fractions of the fluid state
        for(int phaseIdx=0; phaseIdx<numPhases; ++phaseIdx){
                // set the component mole fractions
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    actualFluidState.setMoleFraction(phaseIdx,
                           compIdx,
                           priVars[moleFrac00Idx +
                                   phaseIdx*numComponents +
                                   compIdx]);
                }
            }

//            // For using the ... other way of calculating equilibrium
//             THIS IS ONLY FOR silencing Valgrind but is not used in this model
            for(int phaseIdx=0; phaseIdx<numPhases; ++phaseIdx)
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    const Scalar phi = FluidSystem::fugacityCoefficient(actualFluidState,
                                                                        paramCache,
                                                                        phaseIdx,
                                                                        compIdx);
                    actualFluidState.setFugacityCoefficient(phaseIdx,
                                                      compIdx,
                                                      phi);
            }

            FluidState equilFluidState; // the fluidState *on the interface* i.e. chemical equilibrium
            equilFluidState.assign(actualFluidState) ;
            ConstraintSolver::solve(equilFluidState,
                                    paramCache,
                                    /*setViscosity=*/false,
                                    /*setEnthalpy=*/false) ;

            // Setting the equilibrium composition (in a kinetic model not necessarily the same as the actual mole fraction)
            for(int phaseIdx=0; phaseIdx<numPhases; ++phaseIdx){
                for (int compIdx=0; compIdx< numComponents; ++ compIdx){
                    xEquil_[phaseIdx][compIdx] = equilFluidState.moleFraction(phaseIdx, compIdx);
                }
            }

            // compute densities of all phases
            for(int phaseIdx=0; phaseIdx<numPhases; ++phaseIdx){
                const Scalar rho = FluidSystem::density(actualFluidState, paramCache, phaseIdx);
                actualFluidState.setDensity(phaseIdx, rho);
            }

     }

    template<class ElementSolution>
    void updateTemperatures(const ElementSolution& elemSol,
                            const Problem &problem,
                            const Element& element,
                            const SubControlVolume& scv,
                            FluidState& fluidState)
    {
        for(int phaseIdx=0; phaseIdx < numEnergyEqFluid; ++phaseIdx)
        {
            // retrieve temperature from solution vector
            const Scalar T = BaseType::extractDofPriVars(elemSol, scv)[temperature0Idx + phaseIdx];
            fluidState.setTemperature(phaseIdx, T);
        }
        for(int solidPhaseIdx = numEnergyEqFluid; solidPhaseIdx < numEnergyEqFluid+numEnergyEqSolid; ++solidPhaseIdx)
        {
            temperatureSolid_ = BaseType::extractDofPriVars(elemSol, scv)[temperature0Idx + solidPhaseIdx];
        }
    }

    /*!
     * \brief Returns the temperature in fluid / solid phase(s)
     *        the sub-control volume.
     * \param phaseIdx The local index of the phases
     */
    const Scalar temperatureSolid() const
    { return temperatureSolid_; }
    /*!
     * \brief The mole fraction we would have in the case of chemical equilibrium /
     *        on the interface.
     *
     *     \param phaseIdx The index of the fluid phase
     *     \param compIdx The local index of the component
     */
    const Scalar xEquil(const unsigned int phaseIdx, const unsigned int compIdx) const
    {
        return xEquil_[phaseIdx][compIdx] ;
    }

    /*!
     * \brief The specific interfacial area between two fluid phases [m^2 / m^3]
     *
     * This is _only_ required by the kinetic mass/energy modules
     *
     */
    const Scalar interfacialArea(const unsigned int phaseIIdx, const unsigned int phaseJIdx) const
    {
        // there is no interfacial area between a phase and itself
        assert(phaseIIdx not_eq phaseJIdx);
        return interfacialArea_[phaseIIdx][phaseJIdx];
    }

    //! access function Reynolds Number
    const Scalar reynoldsNumber(const unsigned int phaseIdx) const
    { return reynoldsNumber_[phaseIdx]; }

    //! access function Prandtl Number
    const Scalar prandtlNumber(const unsigned int phaseIdx) const
    { return prandtlNumber_[phaseIdx]; }

    //! access function Nusselt Number
    const Scalar nusseltNumber(const unsigned int phaseIdx) const
    { return nusseltNumber_[phaseIdx]; }

    //! access function Schmidt Number
    const Scalar schmidtNumber(const unsigned int phaseIdx) const
    { return schmidtNumber_[phaseIdx]; }

    //! access function Sherwood Number
    const Scalar sherwoodNumber(const unsigned int phaseIdx) const
    { return sherwoodNumber_[phaseIdx]; }

    //! access function characteristic length
    const Scalar characteristicLength() const
    { return characteristicLength_; }

    //! access function pre factor energy transfer
    const Scalar factorEnergyTransfer() const
    { return factorEnergyTransfer_; }

    //! access function pre factor mass transfer
    const Scalar factorMassTransfer() const
    { return factorMassTransfer_; }

    /*!
     * \brief If running in valgrind this makes sure that all
     *        quantities in the volume variables are defined.
     */
    void checkDefined() const
    {
#if !defined NDEBUG && HAVE_VALGRIND
        Valgrind::CheckDefined(reynoldsNumber_);
        Valgrind::CheckDefined(prandtlNumber_);
        Valgrind::CheckDefined(nusseltNumber_);
        Valgrind::CheckDefined(schmidtNumber_);
        Valgrind::CheckDefined(sherwoodNumber_);
        Valgrind::CheckDefined(characteristicLength_);
        Valgrind::CheckDefined(factorEnergyTransfer_);
        Valgrind::CheckDefined(factorMassTransfer_);
        Valgrind::CheckDefined(interfacialArea_);
#endif
    }

private:
    //! dimensionless numbers
    Scalar reynoldsNumber_[numPhases];
    Scalar prandtlNumber_[numPhases];
    Scalar nusseltNumber_[numPhases];
    Scalar schmidtNumber_[numPhases];
    Scalar sherwoodNumber_[numPhases];
    Scalar characteristicLength_;
    Scalar factorEnergyTransfer_;
    Scalar factorMassTransfer_;
    Scalar solidSurface_ ;
    Scalar interfacialArea_[numPhases+1][numPhases+1];
    Scalar xEquil_[numPhases][numComponents];
    Scalar temperatureSolid_;

};
//
} // namespace Dumux

#endif
