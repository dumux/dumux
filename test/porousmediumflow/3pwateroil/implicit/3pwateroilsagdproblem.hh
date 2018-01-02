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
 *
 * \brief Non-isothermal SAGD problem
 *
 */
#ifndef DUMUX_SAGDPROBLEM_HH
#define DUMUX_SAGDPROBLEM_HH

#include <dumux/porousmediumflow/implicit/problem.hh>

#include <dumux/porousmediumflow/3pwateroil/model.hh>
#include <dumux/material/fluidsystems/h2oheavyoil.hh>
#include "3pwateroilsagdspatialparams.hh"

namespace Dumux
{
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
              Dumux::FluidSystems::H2OHeavyOil<TypeTag, typename GET_PROP_TYPE(TypeTag, Scalar)>);


// Enable gravity
SET_BOOL_PROP(SagdProblem, ProblemEnableGravity, true);

// Use forward differences instead of central differences
SET_INT_PROP(SagdProblem, ImplicitNumericDifferenceMethod, +1);

// Write newton convergence
SET_BOOL_PROP(SagdProblem, NewtonWriteConvergence, false);

SET_BOOL_PROP(SagdProblem, UseSimpleModel, true);

SET_BOOL_PROP(SagdProblem, UseMassOutput, true);
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
class SagdProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::Grid Grid;

    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        pressureIdx = Indices::pressureIdx,
        switch1Idx = Indices::switch1Idx,
        switch2Idx = Indices::switch2Idx,

        energyEqIdx = Indices::energyEqIdx,

        // phase and component indices
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,
        wCompIdx = Indices::wCompIdx,
        nCompIdx = Indices::nCompIdx,

        // Phase State
        wPhaseOnly = Indices::wPhaseOnly,
        wnPhaseOnly = Indices::wnPhaseOnly,
        wgPhaseOnly = Indices::wgPhaseOnly,
        threePhases = Indices::threePhases,

        //contiWEqIdx = Indices::contiWEqIdx,
        //contiNEqIdx = Indices::contiNEqIdx,
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };


    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };

public:

    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */

    SagdProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView), pOut_(4e6)
    {

        maxDepth_ = 400.0; // [m]
        FluidSystem::init();
        totalMassProducedOil_ =0;
        totalMassProducedWater_ =0;

        this->timeManager().startNextEpisode(86400);

        name_ = GET_RUNTIME_PARAM(TypeTag, std::string, Problem.Name);
    }

    bool shouldWriteRestartFile() const
    {
        return 0;
    }


    /*!
     * \brief Called directly after the time integration.
     */
    void postTimeStep()
    {
        double time = this->timeManager().time();
        double dt = this->timeManager().timeStepSize();

        // Calculate storage terms
        PrimaryVariables storage;
        this->model().globalStorage(storage);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0)
        {
            std::cout<<"Storage: " << storage << "Time: " << time+dt << std::endl;
            massBalance.open ("massBalance.txt", std::ios::out | std::ios::app );
                        massBalance << "         Storage       " << storage
                                    << "         Time           " << time+dt
                                    << std::endl;
                        massBalance.close();

        }

            // Calculate storage terms
        PrimaryVariables storageW, storageN;
        //Dune::FieldVector<Scalar, 2> flux(0.0);
        this->model().globalPhaseStorage(storageW, wPhaseIdx);
        this->model().globalPhaseStorage(storageN, nPhaseIdx);

        //mass of Oil
        const Scalar newMassProducedOil_ = massProducedOil_;
        std::cout<<" newMassProducedOil_ : "<< newMassProducedOil_ << " Time: " << time+dt << std::endl;

        totalMassProducedOil_ += newMassProducedOil_;
        std::cout<<" totalMassProducedOil_ : "<< totalMassProducedOil_ << " Time: " << time+dt << std::endl;
        //mass of Water
        const Scalar newMassProducedWater_ = massProducedWater_;
        //std::cout<<" newMassProducedWater_ : "<< newMassProducedWater_ << " Time: " << time+dt << std::endl;

        totalMassProducedWater_ += newMassProducedWater_;
        //std::cout<<" totalMassProducedWater_ : "<< totalMassProducedWater_ << " Time: " << time+dt << std::endl;


        const int timeStepIndex = this->timeManager().timeStepIndex();

        if (timeStepIndex == 0 ||
            timeStepIndex % 100 == 0 ||   //after every 1000000 secs
            this->timeManager().episodeWillBeFinished() ||
            this->timeManager().willBeFinished())
        {
            std::cout<<" totalMassProducedOil_ : "<< totalMassProducedOil_ << " Time: " << time+dt << std::endl;
            std::cout<<" totalMassProducedWater_ : "<< totalMassProducedWater_ << " Time: " << time+dt << std::endl;
        }
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


    void sourceAtPos(PrimaryVariables &values,
                     const GlobalPosition &globalPos) const
    {
        values = 0.0;
    }

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
   void boundaryTypesAtPos(BoundaryTypes &bcTypes,
            const GlobalPosition &globalPos) const
    {
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
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
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
    void solDependentNeumann(PrimaryVariables &values,
                      const Element &element,
                      const FVElementGeometry &fvGeometry,
                      const Intersection &is,
                      const int scvIdx,
                      const int boundaryFaceIdx,
                      const ElementVolumeVariables &elemVolVars) const
    {
        values = 0;

        GlobalPosition globalPos;
        if (isBox)
           globalPos = element.geometry().corner(scvIdx);
        else
           globalPos = is.geometry().center();

        // negative values for injection at injection well
        if (globalPos[1] > 8.5 - eps_ && globalPos[1] < 9.5 + eps_)
        {
            values[Indices::contiNEqIdx] = -0.0;
            values[Indices::contiWEqIdx] = -0.193;//*0.5;   // (55.5 mol*12.5)/3600 mol/s m = 0.193
            values[Indices::energyEqIdx] = -9132;//*0.5;   // J/sec m 9132
        }
        else if (globalPos[1] > 2.5 - eps_ && globalPos[1] < 3.5 + eps_) // production well
        {

            const Scalar elemPressW = elemVolVars[scvIdx].pressure(wPhaseIdx);            //Pressures
            const Scalar elemPressN = elemVolVars[scvIdx].pressure(nPhaseIdx);

            const Scalar densityW = elemVolVars[scvIdx].fluidState().density(wPhaseIdx);  //Densities
            const Scalar densityN = elemVolVars[scvIdx].fluidState().density(nPhaseIdx);

            const Scalar elemMobW = elemVolVars[scvIdx].mobility(wPhaseIdx);      //Mobilities
            const Scalar elemMobN = elemVolVars[scvIdx].mobility(nPhaseIdx);

            const Scalar enthW = elemVolVars[scvIdx].enthalpy(wPhaseIdx);      //Enthalpies
            const Scalar enthN = elemVolVars[scvIdx].enthalpy(nPhaseIdx);

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
            Scalar wT = elemVolVars[scvIdx].temperature(); // well temperature
            if ( wT > 495. )
            {
              qE = qW*0.018*enthW + qN*enthN*0.350 + (wT-495.)*5000.; // ~3x injected enthalpy
              std::cout<< "Cooling now! Extracted enthalpy: " << qE << std::endl;
            } else {
              qE = qW*0.018*enthW + qN*enthN*0.350;
              }


            values[Indices::contiWEqIdx] = qW;
            values[Indices::contiNEqIdx] = qN;
            values[Indices::energyEqIdx] = qE;
            massProducedOil_ = qN;
            massProducedWater_ = qW;
        }
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
    void initialAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
    }

    /*!
     * \brief Return the initial phase state inside a control volume.
     *
     * \param globalPos The global position
     */
    int initialPhasePresenceAtPos(const GlobalPosition &globalPos) const
    {
        return wnPhaseOnly;
    }

private:
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        Scalar densityW = 1000.0;
        values[pressureIdx] = 101300.0 + (maxDepth_ - globalPos[1])*densityW*9.81;

        values[switch1Idx] = 295.13;   // temperature
        values[switch2Idx] = 0.3;   //NAPL saturation
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
