// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/**
 * \file
 * \ingroup OnePNCTests
 * \brief Definition of a problem for the 1pnc problem:
 * Component transport of nitrogen dissolved in the water phase.
 */

#ifndef DUMUX_1P2C_TEST_PROBLEM_HH
#define DUMUX_1P2C_TEST_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/porousmediumflow/problem.hh>


namespace Dumux {

/*!
 * \ingroup OnePNCTests
 * \brief Definition of a problem for the 1pnc problem:
 *  Component transport of nitrogen dissolved in the water phase.
 *
 * Nitrogen is dissolved in the water phase and is transported with the
 * water flow from the left side to the right.
 *
 * The model domain is specified in the input file and
 * we use homogeneous soil properties.
 * Initially, the domain is filled with pure water.
 *
 * At the left side, a Dirichlet condition defines the nitrogen mole fraction.
 * The water phase flows from the left side to the right if the applied pressure
 * gradient is >0. The nitrogen is transported with the water flow
 * and leaves the domain at the boundary, where again Dirichlet boundary
 * conditions are applied.
 *
 * This problem uses the \ref OnePNCModel model.
 */
template <class TypeTag>
class OnePTwoCTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;

    // copy some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    enum
    {
        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,

        // component indices
        H2OIdx = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::H2OIdx),
        N2Idx = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::N2Idx),

        // indices of the equations
        contiH2OEqIdx = Indices::conti0EqIdx + H2OIdx,
        contiN2EqIdx = Indices::conti0EqIdx + N2Idx
    };

    //! Property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();
    static const bool isBox = GetPropType<TypeTag, Properties::GridGeometry>::discMethod == DiscretizationMethods::box;

    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;

public:
    OnePTwoCTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry), useNitscheTypeBc_(getParam<bool>("Problem.UseNitscheTypeBc", false))
    {
        //initialize fluid system
        FluidSystem::init();

        // stating in the console whether mole or mass fractions are used
        if(useMoles)
            std::cout<<"problem uses mole fractions"<<std::endl;
        else
            std::cout<<"problem uses mass fractions"<<std::endl;
    }

    /*!
     * \name Problem parameters
     */
    // \{

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        if(globalPos[0] < eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values = initial_(globalPos);

        // condition for the N2 molefraction at left boundary
        if (globalPos[0] < eps_ )
        {
            values[N2Idx] = 2.0e-5;
        }

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a neumann
     *        boundary segment (implementation for the box scheme).
     *
     * This is the method for the case where the Neumann condition is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param elemFluxVarsCache Flux variables caches for all faces in stencil
     * \param scvf The sub-control volume face
     *
     * For this method, the \a values parameter stores the flux
     * in normal direction of each phase. Negative values mean influx.
     * E.g. for the mass balance that would be the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    template<bool useBox = isBox, std::enable_if_t<useBox, int> = 0>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        // no-flow everywhere except at the right boundary
        NumEqVector values(0.0);
        const auto xMax = this->gridGeometry().bBoxMax()[0];
        const auto& ipGlobal = scvf.ipGlobal();
        if (ipGlobal[0] < xMax - eps_)
            return values;

        // set a fixed pressure on the right side of the domain
        static const Scalar dirichletPressure = 1e5;
        const auto& volVars = elemVolVars[scvf.insideScvIdx()];
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];

        // if specified in the input file, use a Nitsche type boundary condition for the box model.
        if (useNitscheTypeBc_)
        {
            values[contiH2OEqIdx] = (volVars.pressure() - dirichletPressure)*1e7;
            values[contiN2EqIdx] = values[contiH2OEqIdx] * (useMoles ? volVars.moleFraction(0, N2Idx)
                                                                     : volVars.massFraction(0, N2Idx));
            return values;
        }

        // evaluate the pressure gradient
        GlobalPosition gradP(0.0);
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto xIp = scv.dofPosition()[0];
            auto tmp = fluxVarsCache.gradN(scv.localDofIndex());
            tmp *= xIp > xMax - eps_ ? dirichletPressure
                                     : elemVolVars[scv].pressure();
            gradP += tmp;
        }

        // compute flux
        auto phaseFlux = vtmv(scvf.unitOuterNormal(), volVars.permeability(), gradP);
        phaseFlux *= useMoles ? volVars.molarDensity() : volVars.density();
        phaseFlux *= volVars.mobility();

        // set Neumann bc values
        values[contiH2OEqIdx] = phaseFlux;
        // emulate an outflow condition for the component transport on the right side
        values[contiN2EqIdx] = phaseFlux * ( useMoles ? volVars.moleFraction(0, N2Idx) : volVars.massFraction(0, N2Idx) );

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a neumann
     *        boundary segment (implementation for the tpfa scheme).
     */
     template<bool useBox = isBox, std::enable_if_t<!useBox, int> = 0>
     NumEqVector neumann(const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const ElementFluxVariablesCache& elemFluxVarsCache,
                         const SubControlVolumeFace& scvf) const
     {
        // no-flow everywhere except at the right boundary
        NumEqVector values(0.0);
        const auto xMax = this->gridGeometry().bBoxMax()[0];
        const auto& ipGlobal = scvf.ipGlobal();
        if (ipGlobal[0] < xMax - eps_)
            return values;

        // set a fixed pressure on the right side of the domain
        static const Scalar dirichletPressure = 1e5;
        const auto& volVars = elemVolVars[scvf.insideScvIdx()];

        auto d = ipGlobal - element.geometry().center();
        d /= d.two_norm2();

        auto upwindTerm = useMoles ? volVars.molarDensity() : volVars.density();
        upwindTerm *= volVars.mobility();

        const auto tij = vtmv(scvf.unitOuterNormal(), volVars.permeability(), d);
        const auto phaseFlux = -1.0*upwindTerm*tij*(dirichletPressure - volVars.pressure());

        // set Neumann bc values
        values[contiH2OEqIdx] = phaseFlux;
        // emulate an outflow condition for the component transport on the right side
        values[contiN2EqIdx] = phaseFlux * (useMoles ? volVars.moleFraction(0, N2Idx) : volVars.massFraction(0, N2Idx));

        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the source term for all phases within a given
     *        sub-control volume.
     *
     * For this method, the \a priVars parameter stores the rate mass
     * of a component is generated or annihilated per volume
     * unit. Positive values mean that mass is created, negative ones
     * mean that it vanishes.
     *
     * The units must be according to either using mole or mass fractions (mole/(m^3*s) or kg/(m^3*s)).
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    { return NumEqVector(0.0); }

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return initial_(globalPos); }

    // \}

private:
    // the internal method for the initial condition
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars;
        priVars[pressureIdx] = 2e5 - 1e5*globalPos[0]; // initial condition for the pressure
        priVars[N2Idx] = 0.0;  // initial condition for the N2 molefraction
        return priVars;
    }
        static constexpr Scalar eps_ = 1e-6;
        bool useNitscheTypeBc_;
    };

} // end namespace Dumux

#endif
