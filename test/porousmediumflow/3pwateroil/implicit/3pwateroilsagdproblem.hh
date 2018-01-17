// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later vesion.                                      *
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
 * \ingroup ThreePWaterOilTests
 * \brief Non-isothermal SAGD problem
 */
#ifndef DUMUX_SAGDPROBLEM_HH
#define DUMUX_SAGDPROBLEM_HH

#include <dumux/porousmediumflow/problem.hh>

#include <dumux/discretization/box/properties.hh>
#include <dumux/porousmediumflow/3pwateroil/model.hh>
#include <dumux/material/fluidsystems/h2oheavyoilfluidsystem.hh>
#include "3pwateroilsagdspatialparams.hh"

namespace Dumux
{
/*!
 * \file
 * \ingroup ThreePWaterOilTests
 * \brief Non-isothermal SAGD problem
 */
template <class TypeTag>
class SagdProblem;

namespace Properties
{
NEW_TYPE_TAG(SagdProblem, INHERITS_FROM(ThreePWaterOilNI, SagdSpatialParams));
NEW_TYPE_TAG(ThreePWaterOilSagdBoxProblem, INHERITS_FROM(BoxModel, SagdProblem));

// Set the grid type
SET_TYPE_PROP(SagdProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(SagdProblem, Problem, Dumux::SagdProblem<TypeTag>);

// Set the fluid system
SET_TYPE_PROP(SagdProblem,
              FluidSystem,
              Dumux::FluidSystems::H2OHeavyOil<typename GET_PROP_TYPE(TypeTag, Scalar)>);

SET_BOOL_PROP(SagdProblem, OnlyGasPhaseCanDisappear, true);

SET_BOOL_PROP(SagdProblem, UseMoles, true);
}


/*!
 * \ingroup ThreePWaterOilBoxModel
 * \ingroup ImplicitTestProblems
 * \brief Non-isothermal problem where ...
 *
 * This problem uses the \ref ThreePWaterOilModel.
 *
 *  */
template <class TypeTag >
class SagdProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    enum {
        pressureIdx = Indices::pressureIdx,
        switch1Idx = Indices::switch1Idx,
        switch2Idx = Indices::switch2Idx,

        contiWEqIdx = Indices::contiWEqIdx,
        contiNEqIdx = Indices::contiNEqIdx,
        energyEqIdx = Indices::energyEqIdx,

        // phase indices
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,

        // phase state
        wnPhaseOnly = Indices::wnPhaseOnly,

        // world dimension
        dimWorld = GridView::dimensionworld
    };

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NeumannFluxes = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:

    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */

    SagdProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry), pOut_(4e6)
    {

        maxDepth_ = 400.0; // [m]
        FluidSystem::init();
        totalMassProducedOil_ =0;
        totalMassProducedWater_ =0;

        name_ = getParam<std::string>("Problem.Name");
    }


    void episodeEnd()
    {
        // Start new episode if episode is over
        // for first 10 year episode length is 1 year
        this->timeManager().startNextEpisode(3600*24*1.);   //episode length sent to 1 day
            std::cout<<"Episode index is set to: "<<this->timeManager().episodeIndex()<<std::endl;
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string name() const
    { return name_; }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param bcTypes The boundary types for the conservation equations
     * \param globalPos The position for which the bc type should be evaluated
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
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The position for which the bc type should be evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
       return initial_(globalPos);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values The neumann values for the conservation equations
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the box scheme
     * \param is The intersection between element and boundary
     * \param scvIdx The local vertex index
     * \param boundaryFaceIdx The index of the boundary face
     * \param elemVolVars Element volume variables
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    NeumannFluxes neumann(const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const ElementVolumeVariables& elemVolVars,
                          const SubControlVolumeFace& scvf) const
    {
        NeumannFluxes values(0.0);

         const auto& globalPos = scvf.ipGlobal();
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
            //divided by molarMass() of water to convert from kg/m s to mol/m s
            const Scalar qW = (((2*3.1415*0.5*4e-14)/(log(effectiveRadius_/wellRadius))) *
                                densityW * elemMobW * ( elemPressW-pOut_))/0.01801528;
            //divided by molarMass() of HeavyOil to convert from kg/m s to mol/m s
            const Scalar qN = (((2*3.1415*0.5*4e-14)/(log(effectiveRadius_/wellRadius))) *
                                densityN * elemMobN  * (elemPressN-pOut_))/0.35;

            Scalar qE;
            //without cooling:
            // qE = qW*0.018*enthW + qN*enthN*0.350;

            //with cooling: see Diplomarbeit Stefan Roll, Sept. 2015
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
            massProducedOil_ = qN;
            massProducedWater_ = qW;
        }
        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The initial values for the primary variables
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
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
        values[switch2Idx] = 0.3;   //NAPL saturation
        return values;
    }

    Scalar maxDepth_;
    static constexpr Scalar eps_ = 1e-6;
    Scalar pIn_;
    Scalar pOut_;
    Scalar totalMassProducedOil_;
    Scalar totalMassProducedWater_;

    mutable Scalar massProducedOil_;
    mutable Scalar massProducedWater_;

    std::string name_;

    std::ofstream massBalance;
};
} //end namespace

#endif
