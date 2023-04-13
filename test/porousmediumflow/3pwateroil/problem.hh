// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ThreePWaterOilTests
 * \brief Non-isothermal steam-assisted gravity
 *        drainage (SAGD) problem.
 */

#ifndef DUMUX_SAGDPROBLEM_HH
#define DUMUX_SAGDPROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \file
 * \ingroup ThreePWaterOilTests
 * \brief Non-isothermal steam-assisted gravity
 *        drainage (SAGD) problem.
 */
template <class TypeTag >
class SagdProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    enum {
        pressureIdx = Indices::pressureIdx,
        switch1Idx = Indices::switch1Idx,
        switch2Idx = Indices::switch2Idx,

        contiWEqIdx = Indices::conti0EqIdx + FluidSystem::wCompIdx,
        contiNEqIdx = Indices::conti0EqIdx + FluidSystem::nCompIdx,
        energyEqIdx = Indices::energyEqIdx,

        // phase indices
        wPhaseIdx = FluidSystem::wPhaseIdx,
        nPhaseIdx = FluidSystem::nPhaseIdx,

        // phase state
        wnPhaseOnly = Indices::wnPhaseOnly,

        // world dimension
        dimWorld = GridView::dimensionworld
    };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;

public:

    SagdProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry), pOut_(4e6)
    {
        maxDepth_ = 400.0; // [m]
        FluidSystem::init();
    }

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position for which the boundary types are evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        // on bottom
        if (globalPos[1] <  eps_)
            bcTypes.setAllNeumann();

        // on top
        else if (globalPos[1] > 40.0 - eps_)
            bcTypes.setAllNeumann();

        // on bottom other than corners
        else if (globalPos[0] > 60 - eps_ )
            bcTypes.setAllDirichlet();

        // on Left
        else
            bcTypes.setAllNeumann();
        return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
       return initial_(globalPos);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the box scheme
     * \param elemVolVars The element volume variables
     * \param elemFluxVarsCache Flux variables caches for all faces in stencil
     * \param scvf The sub-control volume face
     *
     * Negative values mean influx.
     */
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& globalPos = insideScv.dofPosition();

        // negative values for injection at injection well
        if (globalPos[1] > 8.5 - eps_ && globalPos[1] < 9.5 + eps_)
        {
            values[contiNEqIdx] = -0.0;
            values[contiWEqIdx] = -0.193;//*0.5;   // (55.5 mol*12.5)/3600 mol/s m = 0.193
            values[energyEqIdx] = -9132;//*0.5;   // J/sec m 9132
        }
        else if (globalPos[1] > 2.5 - eps_ && globalPos[1] < 3.5 + eps_) // production well
        {

            const Scalar elemPressW = elemVolVars[scvf.insideScvIdx()].pressure(wPhaseIdx);            //Pressures
            const Scalar elemPressN = elemVolVars[scvf.insideScvIdx()].pressure(nPhaseIdx);

            const Scalar densityW = elemVolVars[scvf.insideScvIdx()].fluidState().density(wPhaseIdx);  //Densities
            const Scalar densityN = elemVolVars[scvf.insideScvIdx()].fluidState().density(nPhaseIdx);

            const Scalar elemMobW = elemVolVars[scvf.insideScvIdx()].mobility(wPhaseIdx);      //Mobilities
            const Scalar elemMobN = elemVolVars[scvf.insideScvIdx()].mobility(nPhaseIdx);

            const Scalar enthW = elemVolVars[scvf.insideScvIdx()].enthalpy(wPhaseIdx);      //Enthalpies
            const Scalar enthN = elemVolVars[scvf.insideScvIdx()].enthalpy(nPhaseIdx);

            const Scalar wellRadius = 0.50 * 0.3048; // 0.50 ft as specified by SPE9


            const Scalar gridHeight_ = 0.5;
            const Scalar effectiveRadius_ = 0.208 * gridHeight_;  //Peaceman's Well Model

            using std::log;
            // divided by molarMass() of water to convert from kg/m s to mol/m s
            const Scalar qW = (((2*3.1415*0.5*4e-14)/(log(effectiveRadius_/wellRadius))) *
                                densityW * elemMobW * ( elemPressW-pOut_))/0.01801528;
            // divided by molarMass() of HeavyOil to convert from kg/m s to mol/m s
            const Scalar qN = (((2*3.1415*0.5*4e-14)/(log(effectiveRadius_/wellRadius))) *
                                densityN * elemMobN  * (elemPressN-pOut_))/0.35;

            Scalar qE;
            // without cooling:
            // qE = qW*0.018*enthW + qN*enthN*0.350;

            // with cooling: see Diplomarbeit Stefan Roll, Sept. 2015
            Scalar wT = elemVolVars[scvf.insideScvIdx()].temperature(); // well temperature
            if ( wT > 495. )
            {
              qE = qW*0.018*enthW + qN*enthN*0.350 + (wT-495.)*5000.; // ~3x injected enthalpy
              std::cout<< "Cooling now! Extracted enthalpy: " << qE << std::endl;
            } else {
              qE = qW*0.018*enthW + qN*enthN*0.350;
              }


            values[contiWEqIdx] = qW;
            values[contiNEqIdx] = qN;
            values[energyEqIdx] = qE;
        }
        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     */
     PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        return initial_(globalPos);
    }

private:
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values.setState(wnPhaseOnly);
        Scalar densityW = 1000.0;
        values[pressureIdx] = 101300.0 + (maxDepth_ - globalPos[1])*densityW*9.81;

        values[switch1Idx] = 295.13;   // temperature
        values[switch2Idx] = 0.3;   // NAPL saturation
        return values;
    }

    Scalar maxDepth_;
    static constexpr Scalar eps_ = 1e-6;
    Scalar pIn_;
    Scalar pOut_;

    std::string name_;
};
} // end namespace Dumux

#endif
