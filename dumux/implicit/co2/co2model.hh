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
 * \brief Adaption of the fully implicit scheme to the CO2Model model.
 */
#ifndef DUMUX_CO2_MODEL_HH
#define DUMUX_CO2_MODEL_HH

#include <dumux/implicit/2p2c/2p2cmodel.hh>

namespace Dumux
{
/*!
 * \ingroup CO2Model
 * \brief Adaption of the BOX or CC scheme to the non-isothermal two-phase two-component flow model.
 *
 *   See TwoPTwoCModel for reference to the equations used.
 *   The CO2 model is derived from the 2p2c model. In the CO2 model the phase switch criterion
 *   is different from the 2p2c model.
 *   The phase switch occurs when the equilibrium concentration
 *   of a component in a phase is exceeded, instead of the sum of the components in the virtual phase
 *   (the phase which is not present) being greater that unity as done in the 2p2c model.
 *   The CO2VolumeVariables do not use a constraint solver for calculating the mole fractions as is the
 *   case in the 2p2c model. Instead mole fractions are calculated in the FluidSystem with a given
 *   temperature, pressurem and salinity.
 *   The model is able to use either mole or mass fractions. The property useMoles can be set to either true or false in the
 * 	 problem file. Make sure that the according units are used in the problem setup. useMoles is set to false by default.
 *
 */

template<class TypeTag>
class CO2Model: public TwoPTwoCModel<TypeTag>
{
    typedef TwoPTwoCModel<TypeTag> ParentType;

     typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
     typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
     typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
     typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

     typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
     enum {
         switchIdx = Indices::switchIdx,

         wPhaseIdx = Indices::wPhaseIdx,
         nPhaseIdx = Indices::nPhaseIdx,
         wCompIdx = Indices::wCompIdx,
         nCompIdx = Indices::nCompIdx,

         wPhaseOnly = Indices::wPhaseOnly,
         nPhaseOnly = Indices::nPhaseOnly,
         bothPhases = Indices::bothPhases,

         pwsn = TwoPTwoCFormulation::pwsn,
         pnsw = TwoPTwoCFormulation::pnsw,
         formulation = GET_PROP_VALUE(TypeTag, Formulation)
     };

     typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
     typedef typename GridView::template Codim<0>::Iterator ElementIterator;
     enum {
         dim = GridView::dimension,
         dimWorld = GridView::dimensionworld
     };

     typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
     typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
     static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
     enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
     enum { dofCodim = isBox ? dim : 0 };

public:


     /*!
      * \brief Update the static data of all vertices in the grid.
      *
      * \param curGlobalSol The current global solution
      * \param oldGlobalSol The previous global solution
      */
     void updateStaticData(SolutionVector &curGlobalSol,
                           const SolutionVector &oldGlobalSol)
     {
         bool wasSwitched = false;
         int succeeded;
         try {
             for (unsigned i = 0; i < ParentType::staticDat_.size(); ++i)
                 ParentType::staticDat_[i].visited = false;

             FVElementGeometry fvGeometry;
             static VolumeVariables volVars;
             ElementIterator eIt = this->gridView_().template begin<0> ();
             const ElementIterator &eEndIt = this->gridView_().template end<0> ();
             for (; eIt != eEndIt; ++eIt)
             {
                 fvGeometry.update(this->gridView_(), *eIt);
                 for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
                 {
                     int globalIdx = this->dofMapper().map(*eIt, scvIdx, dofCodim);

                     if (ParentType::staticDat_[globalIdx].visited)
                         continue;

                     ParentType::staticDat_[globalIdx].visited = true;
                     volVars.update(curGlobalSol[globalIdx],
                                    this->problem_(),
                                    *eIt,
                                    fvGeometry,
                                    scvIdx,
                                    false);
                     const GlobalPosition &globalPos = eIt->geometry().corner(scvIdx);
                     if (primaryVarSwitch_(curGlobalSol,
                         volVars,
                         globalIdx,
                         globalPos))
                     {
                         this->jacobianAssembler().markDofRed(globalIdx);
                         wasSwitched = true;
                     }
                 }
             }
             succeeded = 1;
         }
         catch (Dumux::NumericalProblem &e)
         {
             std::cout << "\n"
             << "Rank " << this->problem_().gridView().comm().rank()
             << " caught an exception while updating the static data." << e.what()
             << "\n";
             succeeded = 0;
         }
         //make sure that all processes succeeded. If not throw a NumericalProblem to decrease the time step size.
         if (this->gridView_().comm().size() > 1)
             succeeded = this->gridView_().comm().min(succeeded);

         if (!succeeded) {
             if(this->problem_().gridView().comm().rank() == 0)
                 DUNE_THROW(NumericalProblem,
                            "A process did not succeed in updating the static data.");
                 return;
         }

         // make sure that if there was a variable switch in an
         // other partition we will also set the switch flag
         // for our partition.
         if (this->gridView_().comm().size() > 1)
             wasSwitched = this->gridView_().comm().max(wasSwitched);

         ParentType::setSwitched_(wasSwitched);
     }

 protected:


     /*!
      * \brief Performs variable switch at a vertex, returns true if a
      *        variable switch was performed.
      */
     bool primaryVarSwitch_(SolutionVector &globalSol,
                              const VolumeVariables &volVars, int globalIdx,
                              const GlobalPosition &globalPos)
       {
         typename FluidSystem::ParameterCache paramCache;
           // evaluate primary variable switch
           bool wouldSwitch = false;
           int phasePresence = ParentType::staticDat_[globalIdx].phasePresence;
           int newPhasePresence = phasePresence;

           // check if a primary var switch is necessary
           if (phasePresence == nPhaseOnly)
           {

               Scalar xnw = volVars.moleFraction(nPhaseIdx, wCompIdx);
               Scalar xnwMax = FluidSystem::equilibriumMoleFraction(volVars.fluidState(), paramCache, nPhaseIdx);

               if(xnw > xnwMax)
                   wouldSwitch = true;

               if (ParentType::staticDat_[globalIdx].wasSwitched)
                   xnwMax *= 1.02;

               //If mole fraction is higher than the equilibrium mole fraction make a phase switch
               if(xnw > xnwMax)
               {
                   // wetting phase appears
                   std::cout << "wetting phase appears at vertex " << globalIdx
                             << ", coordinates: " << globalPos << ", xnw > xnwMax: "
                             << xnw << " > "<< xnwMax << std::endl;
                   newPhasePresence = bothPhases;
                   if (formulation == pnsw)
                       globalSol[globalIdx][switchIdx] = 0.0;
                   else if (formulation == pwsn)
                       globalSol[globalIdx][switchIdx] = 1.0;
               }
           }
           else if (phasePresence == wPhaseOnly)
           {

               Scalar xwn = volVars.moleFraction(wPhaseIdx, nCompIdx);
               Scalar xwnMax = FluidSystem::equilibriumMoleFraction(volVars.fluidState(), paramCache, wPhaseIdx);

               //If mole fraction is higher than the equilibrium mole fraction make a phase switch
               if(xwn > xwnMax)
                   wouldSwitch = true;
               if (ParentType::staticDat_[globalIdx].wasSwitched)
                   xwnMax *= 1.02;


               if(xwn > xwnMax)
               {
                   // non-wetting phase appears
                   std::cout << "non-wetting phase appears at vertex " << globalIdx
                             << ", coordinates: " << globalPos << ", xwn > xwnMax: "
                             << xwn << " > "<< xwnMax << std::endl;

                   newPhasePresence = bothPhases;
                   if (formulation == pnsw)
                       globalSol[globalIdx][switchIdx] = 0.999;
                   else if (formulation == pwsn)
                       globalSol[globalIdx][switchIdx] = 0.001;
               }
           }
           else if (phasePresence == bothPhases)
           {
               Scalar Smin = 0.0;
               if (ParentType::staticDat_[globalIdx].wasSwitched)
                   Smin = -0.01;

               if (volVars.saturation(nPhaseIdx) <= Smin)
               {
                   wouldSwitch = true;
                   // nonwetting phase disappears
                   std::cout << "Nonwetting phase disappears at vertex " << globalIdx
                             << ", coordinates: " << globalPos << ", sn: "
                             << volVars.saturation(nPhaseIdx) << std::endl;
                   newPhasePresence = wPhaseOnly;

                   if(!useMoles) //mass-fraction formulation
                   {
					   globalSol[globalIdx][switchIdx]
						   = volVars.massFraction(wPhaseIdx, nCompIdx);
                   }
                   else //mole-fraction formulation
                   {
					   globalSol[globalIdx][switchIdx]
					       = volVars.moleFraction(wPhaseIdx, nCompIdx);
                   }
               }
               else if (volVars.saturation(wPhaseIdx) <= Smin)
               {
                   wouldSwitch = true;
                   // wetting phase disappears
                   std::cout << "Wetting phase disappears at vertex " << globalIdx
                             << ", coordinates: " << globalPos << ", sw: "
                             << volVars.saturation(wPhaseIdx) << std::endl;
                   newPhasePresence = nPhaseOnly;

                   if(!useMoles) //mass-fraction formulation
                   {
					   globalSol[globalIdx][switchIdx]
						   = volVars.massFraction(nPhaseIdx, wCompIdx);
                   }
                   else //mole-fraction formulation
                   {
						globalSol[globalIdx][switchIdx]
						= volVars.moleFraction(nPhaseIdx, wCompIdx);
                   }
               }
           }

           ParentType::staticDat_[globalIdx].phasePresence = newPhasePresence;
           ParentType::staticDat_[globalIdx].wasSwitched = wouldSwitch;
           return phasePresence != newPhasePresence;
       }


};

}

#endif
