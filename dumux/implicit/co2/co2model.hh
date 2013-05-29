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
 * \brief Adaption of the BOX or CC scheme to the two-phase two-component flow model without constraint solver.
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
 *
 */

template<class TypeTag>
class CO2Model: public TwoPTwoCModel<TypeTag>
{


    typedef TwoPTwoCModel<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, BaseModel) BaseType;

     typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
     typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
     typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
     typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
     typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
     typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
     typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
     typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
     enum {
         numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
         numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
     };

     typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
     enum {
         pressureIdx = Indices::pressureIdx,
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
     typedef typename GridView::ctype CoordScalar;
     typedef typename GridView::template Codim<0>::Entity Element;
     typedef typename GridView::template Codim<0>::Iterator ElementIterator;
     enum {
         dim = GridView::dimension,
         dimWorld = GridView::dimensionworld
     };
     typedef typename GridView::template Codim<dim>::Entity Vertex;
     typedef typename GridView::template Codim<dim>::Iterator VertexIterator;

     typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
     typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;
     typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
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

         for (unsigned i = 0; i < ParentType::staticDat_.size(); ++i)
             ParentType::staticDat_[i].visited = false;

         FVElementGeometry fvGeometry;
         static VolumeVariables volVars;
         ElementIterator elemIt = this->gridView_().template begin<0> ();
         const ElementIterator &elemEndIt = this->gridView_().template end<0> ();
         for (; elemIt != elemEndIt; ++elemIt)
         {
             fvGeometry.update(this->gridView_(), *elemIt);
             for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
             {
                 int globalIdx = this->dofMapper().map(*elemIt, scvIdx, dofCodim);

                 if (ParentType::staticDat_[globalIdx].visited)
                     continue;

                 ParentType::staticDat_[globalIdx].visited = true;
                 volVars.update(curGlobalSol[globalIdx],
                                this->problem_(),
                                *elemIt,
                                fvGeometry,
                                scvIdx,
                                false);
                 const GlobalPosition &globalPos = elemIt->geometry().corner(scvIdx);
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

         // make sure that if there was a variable switch in an
         // other partition we will also set the switch flag
         // for our partition.
         if (this->gridView_().comm().size() > 1)
             wasSwitched = this->gridView_().comm().max(wasSwitched);

         ParentType::setSwitched_(wasSwitched);
     }

 protected:


     /*!
      * \brief Set the old phase of all verts state to the current one.
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

               Scalar xnw = volVars.fluidState().moleFraction(nPhaseIdx, wCompIdx);
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

               Scalar xwn = volVars.fluidState().moleFraction(wPhaseIdx, nCompIdx);
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

                   globalSol[globalIdx][switchIdx]
                       = volVars.fluidState().massFraction(wPhaseIdx, nCompIdx);
               }
               else if (volVars.saturation(wPhaseIdx) <= Smin)
               {
                   wouldSwitch = true;
                   // wetting phase disappears
                   std::cout << "Wetting phase disappears at vertex " << globalIdx
                             << ", coordinates: " << globalPos << ", sw: "
                             << volVars.saturation(wPhaseIdx) << std::endl;
                   newPhasePresence = nPhaseOnly;

                   globalSol[globalIdx][switchIdx]
                       = volVars.fluidState().massFraction(nPhaseIdx, wCompIdx);
               }
           }

           ParentType::staticDat_[globalIdx].phasePresence = newPhasePresence;
           ParentType::staticDat_[globalIdx].wasSwitched = wouldSwitch;
           return phasePresence != newPhasePresence;
       }


};

}

#endif
