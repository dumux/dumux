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
 * \brief Non-isothermal steam injection problem where superheated water is injected under high pressure and vaporizes after some distance from the injection well
 *
 */
#ifndef DUMUX_SUPERHEATEDPROBLEM_HH
#define DUMUX_SUPERHEATEDPROBLEM_HH

#include <dumux/implicit/2p1cni/2p1cmodel.hh>
#include <dumux/implicit/common/implicitporousmediaproblem.hh>

#include <dune/grid/uggrid.hh>

#include "conceptspatialparams.hh"
#include <dumux/material/fluidsystems/purewaterfluidsystem.hh>


#include <dumux/linear/amgbackend.hh>

#define ISOTHERMAL 0


namespace Dumux
{
template <class TypeTag>
class InjectionProblem;

namespace Properties
{
NEW_TYPE_TAG(InjectionProblem, INHERITS_FROM(TwoPOneCNI, InjectionProblemSpatialParams));
NEW_TYPE_TAG(InjectionBoxProblem, INHERITS_FROM(BoxModel, InjectionProblem));

 SET_PROP(InjectionProblem, Grid)
{
    typedef Dune::UGGrid<3> type;
};



// Set the problem property
SET_PROP(InjectionProblem, Problem)
{
    typedef Dumux::InjectionProblem<TypeTag> type;
};

// Set the fluid system
SET_TYPE_PROP(InjectionProblem,
              FluidSystem,
              Dumux::FluidSystems::PureWaterSimpleFluidSystem<typename GET_PROP_TYPE(TypeTag, Scalar), /*useComplexRelations*/true>);


// Enable gravity
SET_BOOL_PROP(InjectionProblem, ProblemEnableGravity, true);

// Use forward differences instead of central differences
SET_INT_PROP(InjectionProblem, ImplicitNumericDifferenceMethod, +1);

// Write newton convergence
SET_BOOL_PROP(InjectionProblem, NewtonWriteConvergence, false);

//Define whether spurious cold-water flow into the steam is blocked
SET_BOOL_PROP(InjectionProblem, UseBlockingOfSpuriousFlow, true);

// Linear solver settings
SET_TYPE_PROP(InjectionProblem, LinearSolver, Dumux::AMGBackend<TypeTag> );

}


/*!
 * \ingroup TwoPOneCNIBoxModel
 * \ingroup ImplicitTestProblems
 * \brief Non-isothermal problem where superheated water is injected under high pressure and vaporizes after some distance from the injection well.
 *
 * This problem uses the \ref TwoPOneCNIModel.
 *
 *  */
template <class TypeTag >
class InjectionProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::Grid Grid;
    typedef Dumux::H2O<Scalar> IapwsH2O;

    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        pressureIdx = Indices::pressureIdx,
        switch1Idx = Indices::switch1Idx,

        energyEqIdx = Indices::energyEqIdx,

        // phase and component indices
        wPhaseIdx = Indices::wPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,


        // Phase State
        wPhaseOnly = Indices::wPhaseOnly,
        gPhaseOnly = Indices::gPhaseOnly,
        twoPhases = Indices::twoPhases,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };


    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) LocalResidual;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;
     typedef Dune::FieldVector<Scalar,dim> Vector;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef Dune::array<Scalar,7> EpisodeArray;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, 4> Flux;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };

public:

    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */

    InjectionProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView), eps_(1e-6)/*, pOut_(4e6)*/
    {

    name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                             std::string,
                                             Problem,
                                             Name);


    episodes_[0] = 1.0; //hours 6
    episodes_[1] = 2.0; //10
    episodes_[2] = 3.0; // 12
    episodes_[3] = 6.0;// 22
    episodes_[4] = 12.0; // 24
    episodes_[5] = 24.0; // 32
    episodes_[6] = 480.0; // 36


    writeAllSteps_ = GET_RUNTIME_PARAM(TypeTag, bool, Problem.WriteAllTimeSteps);

    preheatZone_ = GET_RUNTIME_PARAM(TypeTag, std::string, Problem.PreheatZone);

    massRate_ = GET_RUNTIME_PARAM(TypeTag, Scalar, Problem.MassRate);

   myMassrate_ = 0;
   myP_ = 0;
   myEnth_ = 0;

    FluidSystem::init(
        /*tempMin=*/273.15,
             /*tempMax=*/623.15,
             /*numTemp=*/300,
             /*pMin=*/0.0,
             /*pMax=*/20e6,
             /*numP=*/300);

    logfile_.open("protokoll.txt"); //for the logfile
    logfile_<<"Logfile: "  << std::endl;
    logfile_<<"Time [s]:  Enthalpy   Mass Rate   Pressure  massCalc  "  << std::endl;

    this->timeManager().startNextEpisode(episodes_[0]*3600); //set the first episode lenght

    } //end of constructor


     ~InjectionProblem() //The deconstructor
   {
     logfile_.close();
   }


 bool shouldWriteOutput() const //define output
    {
        if(this->timeManager().episodeWillBeOver())
          {
            return true;
          }

        else if(writeAllSteps_)
            {
              return true;
            }

        else
        {
            return  (this->timeManager().timeStepIndex() == 0) ||
                    this->timeManager().time() == episodes_[5] ||
                    this->timeManager().willBeFinished();
        }

    }


 void episodeEnd() //Episode Management
      {
                std::cout<<"Episode index before end of current episode is set to: "<<this->timeManager().episodeIndex()<<std::endl;
                Scalar episodeLenght = episodes_[this->timeManager().episodeIndex()]*3600 - this->timeManager().time();
                this->timeManager().startNextEpisode(episodeLenght);
                std::cout<<"Episode index of new episode is set to: "<<this->timeManager().episodeIndex()<<std::endl;
                std::cout<<"New episode length is set to: "<<episodeLenght<<std::endl;
      }

 bool shouldWriteRestartFile() const
{
return false;
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


    //Writes out storage terms

       void postTimeStep()
    {

     // Calculate storage terms
        PrimaryVariables storage;
        Flux flux(0.0);

        this->model().globalStorage(storage);

        //calculating flux across the layer of the coordinate 2 (Z axis)
        //intercepted at coordinate value = 999m (just below the top of the cap rock)
        calculateFluxAcrossLayer(flux, 0, /*depthFluxMeasurement_*/0.135);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout<<"Storage: mass=[" << storage[0]  << "]"
            << " energy=[" << storage[1] << "]"<< "flux=[" <<flux << "] \n";
        }

        std::cout << "Enthalpy [J/s] " << myEnth_ << std::endl;
        std::cout << "MassRate [kg/s] " << myMassrate_ << std::endl;
        std::cout << "Pressure [Pa] " << myP_ << std::endl;
        std::cout << "Mass Calc [kg/s] = " << flux[0]+flux[1] << std::endl;

        logfile_ << this->timeManager().time()+ this->timeManager().timeStepSize() << "  " << myEnth_ << "  " << myMassrate_ << "  " << myP_  << "  " << flux[0]+flux[1] << std::endl;

    }


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
     * \param values The boundary types for the conservation equations
     * \param globalPos The position for which the bc type should be evaluated
     */
   void boundaryTypesAtPos(BoundaryTypes &bcTypes,
            const GlobalPosition &globalPos) const
    {
         if (onBackBoundary_(globalPos) || onTopBoundary_(globalPos)) {
           bcTypes.setAllDirichlet();
        }
        else {
           bcTypes.setAllNeumann();
        }
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
     * \param fvGeomtry The finite-volume geometry in the box scheme
     * \param is The intersection between element and boundary
     * \param scvIdx The local vertex index
     * \param boundaryFaceIdx The index of the boundary face
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */

        void solDependentNeumann (PrimaryVariables &values,
                  const Element &element,
                  const FVElementGeometry &fvGeometry,
                  const Intersection &intersection,
                  const int scvIdx,
                  const int boundaryFaceIdx,
                  const ElementVolumeVariables &elemVolVars) const
    {

        const GlobalPosition &globalPos =
            fvGeometry.boundaryFace[boundaryFaceIdx].ipGlobal;

        values = 0.0;

           if(onInjectionBoundary_(globalPos))
             {
               Scalar pElemW = elemVolVars[scvIdx].pressure(wPhaseIdx);
               Scalar pWell = 2e5;

               Scalar enthalpy = IapwsH2O::liquidEnthalpy(273+110, pWell);
//                std::cout << "Enthalpy = " << enthalpy << std::endl;

               Scalar massRate = (pWell - pElemW);
//                std::cout << "MassRate = " << massRate << std::endl;

               Scalar area = 1.0 * 2 * 3.141592 * 0.135 /360; //injection area [m^2] (circumference * height)
               values[Indices::conti0EqIdx] = -massRate/area; //[kg/(s m^2)]

               values[Indices::energyEqIdx] = -massRate/area * enthalpy; //[J/s m^2] (mass * spec. enthalpy of steam)

    myMassrate_ = massRate ;
    myP_ = pElemW;
    myEnth_ = enthalpy*massRate;

    if (this->gridView().comm().size() > 1)
    {
     myMassrate_ = this->gridView().comm().sum(myMassrate_);
     myP_ = this->gridView().comm().sum(myP_);
     myEnth_ = this->gridView().comm().sum(myEnth_);
    }

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
     * \param vert The vertex
     * \param globalIdx The index of the global vertex
     * \param globalPos The global position
     */
    int initialPhasePresence(const Vertex &vert,
                             int &globalIdx,
                             const GlobalPosition &globalPos) const
    {
        return wPhaseOnly;
      //return gPhaseOnly;
      //return twoPhases;
    }





private:
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        Scalar densityW = 1000.0;
        values[pressureIdx] = 101300. + (this->bBoxMax()[dim-1] - globalPos[2])*densityW*9.81; //non-wetting phase pressure (is equal to wetting phase pressure for Sw = 1 if VanGenuchten Law is used)
        values[switch1Idx] = 283.15;

        //choose the preheat zone via input file parameter
        if(preheatZone_ == "full") //preheat everything
          {
           values[switch1Idx] = 273.15 + 90.0; //90 °C preheating temperature
          }


        else if (preheatZone_ == "half") //preheat the lower half of the domain
         {
            if(globalPos[2]< 3.5 +eps_ )
              {
               values[switch1Idx] = 273.15 + 90.0;
              }

            else {}
         }

        else if(preheatZone_ == "stripe") //preheat a stripe
        {
            if(globalPos[2]< 3.0 +eps_ && globalPos[2] > 2.0 - eps_)
              {
               values[switch1Idx] = 273.15 + 90.0;
              }

            else {}
        }

        else {}
    }

     GlobalPosition convertToPolarCoordinates_(const GlobalPosition &globalPos) const
    {
        GlobalPosition polarPos;
        polarPos[0] = sqrt(pow(globalPos[0],2) + pow(globalPos[1],2)); //radius
        polarPos[1] = acos(globalPos[0]/polarPos[0]); //angle
        polarPos[2] = globalPos[2]; //z-coordinate
        return polarPos;
    }

        bool onSideBoundary_(const GlobalPosition &globalPos) const
    {
        Scalar angle = convertToPolarCoordinates_(globalPos)[1];
        return angle < eps_ || angle > 30 - eps_;
    }

    bool onTopBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[dim-1] > this->bBoxMax()[dim-1] - eps_ ;
    }

    bool onInjectionBoundary_(const GlobalPosition &globalPos) const
    {
        Scalar radius = convertToPolarCoordinates_(globalPos)[0];
        return radius < 0.135 + eps_  && globalPos[dim-1] > 2.0 -eps_ && globalPos[dim-1] < 3.0 + eps_ ;
    }

    bool onBackBoundary_(const GlobalPosition &globalPos) const
    {
        Scalar radius = convertToPolarCoordinates_(globalPos)[0];
        return radius > 5.0 - eps_;
    }

         /*!
      * \brief Calculate the fluxes across a certain layer in the domain.
      * The layer is situated perpendicular to the coordinate axis "coord" and cuts
      * the axis at the value "coordVal".
      *
      * \param globalSol The global solution vector
      * \param flux A vector to store the flux
      * \param axis The dimension, perpendicular to which the layer is situated
      * \param coordVal The (Scalar) coordinate on the axis, at which the layer is situated
      */
     void calculateFluxAcrossLayer(Flux &flux, int coord,
             Scalar coordVal)
     {
         //Scalar massUpwindWeight = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);
         Scalar massUpwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);
         ElementVolumeVariables elemVolVars;
         FVElementGeometry fvGeometry;

         ElementIterator elemIt = this->gridView().template begin<0>();
         const ElementIterator &endit = this->gridView().template end<0> ();

         GlobalPosition globalI, globalJ;
         Flux tmpFlux(0.0);
         int sign;


         // Loop over elements
         for (; elemIt != endit; ++elemIt)
         {
            if(elemIt->partitionType() == Dune::InteriorEntity){
             fvGeometry.update(this->gridView(), *elemIt);
             elemVolVars.update(*this,
                     *elemIt,
                     fvGeometry,
                     false /* oldSol? */);
             if (elemIt->partitionType() != Dune::InteriorEntity)
                 continue;

             for (int faceId = 0; faceId< fvGeometry.numScvf; faceId++)
             {
                 int idxI = fvGeometry.subContVolFace[faceId].i;
                 int idxJ = fvGeometry.subContVolFace[faceId].j;
                 int flagI, flagJ;

                 globalI = fvGeometry.subContVol[idxI].global;
                 globalJ = fvGeometry.subContVol[idxJ].global;

                 // 2D case: give y or x value of the line over which flux is to be
                 //            calculated.
                 // up to now only flux calculation to lines or planes (3D) parallel to
                 // x, y and z axis possible

                 // Flux across plane with z = 80 numEq
                 if (globalI[coord] < coordVal + 0.5e-2  )
                     flagI = 1;
                 else
                     flagI = -1;

                 if (globalJ[coord] < coordVal)
                     flagJ = 1;
                 else
                     flagJ = -1;

                 if (flagI == flagJ)
                 {
                     sign = 0;
                 }
                 else
                 {
                     if (flagI > 0)
                         sign = -1;
                     else
                         sign = 1;
                 }
                 // get variables
                 if (flagI != flagJ)
                 {
                     FluxVariables fluxVars(*this,
                             *elemIt,
                             fvGeometry,
                             faceId,
                             elemVolVars);
                     tmpFlux = 0;


                     ////////
                     // advective fluxes of all components in all phases
                     ////////

                     // loop over all phases
                           for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                           {

                               // data attached to upstream and the downstream vertices
                               // of the current phase
                               const VolumeVariables &up = elemVolVars[fluxVars.upstreamIdx(phaseIdx)];
                               const VolumeVariables &dn = elemVolVars[fluxVars.downstreamIdx(phaseIdx)];

                               // add advective flux of current phase (mass)
                               tmpFlux[phaseIdx] +=
                                   fluxVars.volumeFlux(phaseIdx)
                                  *
                                   ((    massUpwindWeight_)*up.density(phaseIdx)
                                    +
                                    (1 - massUpwindWeight_)*dn.density(phaseIdx));

                                   //energy
                                 tmpFlux[phaseIdx+numPhases] +=
                                   fluxVars.volumeFlux(phaseIdx)
                                  *
                                   ((    massUpwindWeight_)*up.density(phaseIdx)*up.enthalpy(phaseIdx)
                                    +
                                    (1 - massUpwindWeight_)*dn.density(phaseIdx)*dn.enthalpy(phaseIdx));

                           }


                     // the face normal points into the outward direction, so we
                     // have to multiply all fluxes with -1
                     tmpFlux *= -1;
                     tmpFlux *= sign;
                     flux += tmpFlux;

                 }
             }
            }
         }
         //If parallel sum of fluxes over all processors
         if (this->gridView().comm().size() > 1)
//          flux = this->problem_().gridView().comm().sum(flux);
         flux = this->gridView().comm().sum(flux);
     }


    Scalar eps_;
    std::string name_;
   mutable Scalar myMassrate_;
   mutable  Scalar myP_;
   mutable  Scalar myEnth_;


    EpisodeArray episodes_; //Array storing the episode lengths for vtu output
    bool writeAllSteps_;

    std::ofstream logfile_;
    std::string preheatZone_;
    Scalar massRate_;


};
} //end namespace

#endif
