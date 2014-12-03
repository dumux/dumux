// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
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
 * \brief This file adds kinetic mass and energy transfer modules to
 *        the M-phase N-component model.
 */
#ifndef DUMUX_MPNC_MODEL_KINETIC_HH
#define DUMUX_MPNC_MODEL_KINETIC_HH

// equilibrium model
#include <dumux/implicit/mpnc/mpncmodel.hh>

// generic stuff
#include "mpncpropertieskinetic.hh"

namespace Dumux
{
/*!
 * \ingroup MPNCModel
 * \brief A fully implicit model for MpNc flow using
 *        vertex centered finite volumes.
 *        This is the specialization that is able to capture kinetic mass and / or energy transfer.
 *
 *        Please see the comment in the localresidualenergykinetic class about the inclusion of volume changing work.
 *
 *        Please see the comment about different possibilities of calculating phase-composition in the mpnclocalresidualmasskinetic class.
 */

template<class TypeTag>
class MPNCModelKinetic : public MPNCModel<TypeTag>
{
    typedef MPNCModel<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;
    typedef Dumux::MPNCVtkWriter<TypeTag> MPNCVtkWriter;

    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy)};
    enum { enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion)};
    enum { enableKinetic = GET_PROP_VALUE(TypeTag, EnableKinetic)};
    enum { numEnergyEquations = GET_PROP_VALUE(TypeTag, NumEnergyEquations)};
    enum { enableSmoothUpwinding = GET_PROP_VALUE(TypeTag, ImplicitEnableSmoothUpwinding)};
    enum { enablePartialReassemble = GET_PROP_VALUE(TypeTag, ImplicitEnablePartialReassemble)};
    enum { enableJacobianRecycling = GET_PROP_VALUE(TypeTag, ImplicitEnableJacobianRecycling)};
    enum { numDiffMethod = GET_PROP_VALUE(TypeTag, ImplicitNumericDifferenceMethod)};
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases)};
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents)};
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq)};
    enum { numEnergyEqs = Indices::numPrimaryEnergyVars};
    enum { dimWorld = GridView::dimensionworld};
    enum { dim = GridView::dimension};

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::BlockVector<GlobalPosition>                GlobalPositionField;
    typedef Dune::array<GlobalPositionField, numPhases>      PhaseGlobalPositionField;

    typedef std::vector<Dune::FieldVector<Scalar, 1> >  ScalarVector;
    typedef Dune::array<ScalarVector, numPhases>         PhaseVector;
    typedef Dune::FieldVector<Scalar, dim>              DimVector;
    typedef Dune::BlockVector<DimVector>                DimVectorField;
    typedef Dune::array<DimVectorField, numPhases>       PhaseDimVectorField;

public:
    /*!
     * \brief Initialize the fluid system's static parameters
     * \param problem The Problem
     */
    void init(Problem &problem)
    {
        ParentType::init(problem);
        static_assert(FluidSystem::numComponents==2, " so far kinetic mass transfer assumes a two-component system. ");
    }

    /*!
     * \brief Initialze the velocity vectors.
     *
     *        Otherwise the initial solution cannot be written to disk.
     */
    void initVelocityStuff(){
        // belongs to velocity averaging
        int numVertices = this->gridView().size(dim);

        // allocation and bringing to size
        for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            volumeDarcyVelocity_[phaseIdx].resize(numVertices);
            volumeDarcyMagVelocity_[phaseIdx].resize(numVertices);
            std::fill(volumeDarcyMagVelocity_[phaseIdx].begin(), volumeDarcyMagVelocity_[phaseIdx].end(), 0.0);
            boxSurface_.resize(numVertices);
            std::fill(boxSurface_.begin(), boxSurface_.end(), 0.0);
            volumeDarcyVelocity_[phaseIdx] = 0;
        }
    }


    /*!
     * \brief face-area weighted average of the velocities.
     *
     *      This function calculates the darcy velocities on the faces and averages them for one vertex.
     *      The weight of the average is the area of the face.
     *
     *      Called by newtonControlle in newtonUpdate if velocityAveragingInProblem is true.
     */
    void calcVelocityAverage()
    {
        Scalar numVertices = this->gridView().size(dim);

        // reset
        for (int phaseIdx =0; phaseIdx<numPhases; ++phaseIdx){
            std::fill(boxSurface_.begin(), boxSurface_.end(), 0.0);
            volumeDarcyVelocity_[phaseIdx] = 0;
            std::fill(volumeDarcyMagVelocity_[phaseIdx].begin(), volumeDarcyMagVelocity_[phaseIdx].end(), 0.0);
        }

        // loop all elements
        ElementIterator eIt     = this->gridView_().template begin<0>();
        ElementIterator eEndIt  = this->gridView_().template end<0>();
        for (; eIt != eEndIt; ++eIt){
            //obtaining the elementVolumeVariables needed for the Fluxvariables
            FVElementGeometry fvGeometry;
            ElementVolumeVariables elemVolVars;
            ElementBoundaryTypes elemBcTypes;

            fvGeometry.update(this->gridView_(), *eIt);
            elemBcTypes.update(this->problem_(), *eIt);

            this->setHints(*eIt, elemVolVars);
            elemVolVars.update(this->problem_(),
                               *eIt,
                               fvGeometry,
                               false);

            this->updateCurHints(*eIt, elemVolVars);
            for (int fIdx = 0; fIdx < fvGeometry.numScvf; ++ fIdx) {
                int i = fvGeometry.subContVolFace[fIdx].i;
                int I = this->vertexMapper().map(*eIt, i, dim);

                int j = fvGeometry.subContVolFace[fIdx].j;
                int J = this->vertexMapper().map(*eIt, j, dim);

                const Scalar scvfArea     = fvGeometry.subContVolFace[fIdx].normal.two_norm();
                boxSurface_[I]      += scvfArea;
                boxSurface_[J]      += scvfArea;

                FluxVariables fluxVars(this->problem_(),
                            *eIt,
                            fvGeometry,
                            fIdx,
                            elemVolVars);

                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    GlobalPosition faceDarcyVelocity = fluxVars.velocity(phaseIdx);
                    faceDarcyVelocity   *= scvfArea;

                    // although not yet really a volume average, see later: divide by surface
                    volumeDarcyVelocity_[phaseIdx][I] += faceDarcyVelocity;
                    volumeDarcyVelocity_[phaseIdx][J] += faceDarcyVelocity;
                } // end phases
            } // end faces
        } // end elements

        // Divide by the sum of the faces: actually producing the average (as well as it's magnitude)
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // first, divide the velocity field by the
            // respective finite volume's surface area
            for (int I = 0; I < numVertices; ++I){
             volumeDarcyVelocity_[phaseIdx][I]      /= boxSurface_[I];
             volumeDarcyMagVelocity_[phaseIdx][I]   = volumeDarcyVelocity_[phaseIdx][I].two_norm() ;
            }// end all vertices
        }// end all phases
    }// end calcVelocity

//    /*!
//     * \brief Check whether the current solution makes sense.
//     */
//    void checkPlausibility() const
//    {
//        // Looping over all elements of the domain
//        ElementIterator eEndIt = this->problem_().gridView().template end<0>();
//        for (ElementIterator eIt = this->problem_().gridView().template begin<0>() ; eIt not_eq eEndIt; ++eIt)
//        {
//            ElementVolumeVariables elemVolVars;
//            FVElementGeometry fvGeometry;
//
//            // updating the volume variables
//            fvGeometry.update(this->problem_().gridView(), *eIt);
//            elemVolVars.update(this->problem_(), *eIt, fvGeometry, false);
//
//            std::stringstream  message ;
//            // number of scv
//            const unsigned int numScv = fvGeometry.numScv; // box: numSCV, cc:1
//
//            for (unsigned int scvIdx = 0; scvIdx < numScv; ++scvIdx) {
//
//                const FluidState & fluidState = elemVolVars[scvIdx].fluidState();
//
//                // energy check
//                for(unsigned int energyEqIdx=0; energyEqIdx<numEnergyEqs; energyEqIdx++){
//                    const Scalar eps = 1e-6 ;
////                    const Scalar temperatureTest = elemVolVars[scvIdx].fluidState().temperature();
//                    const Scalar temperatureTest = elemVolVars[scvIdx].temperature(energyEqIdx);
////                    const Scalar temperatureTest = 42;
//
//                    if (not std::isfinite(temperatureTest) or temperatureTest < 0. ){
//                        message <<"\nUnphysical Value in Energy: \n";
//                        message << "\tT" <<"_"<<FluidSystem::phaseName(energyEqIdx)<<"="<< temperatureTest <<"\n";
//                    }
//                }
//
//                // mass Check
//                for(int phaseIdx=0; phaseIdx<numPhases; phaseIdx++){
//                    const Scalar eps = 1e-6 ;
//                    for (int compIdx=0; compIdx< numComponents; ++ compIdx){
//                        const Scalar xTest = fluidState.moleFraction(phaseIdx, compIdx);
//                        if (not std::isfinite(xTest) or xTest < 0.-eps or xTest > 1.+eps ){
//                            message <<"\nUnphysical Value in Mass: \n";
//
//                            message << "\tx" <<"_"<<FluidSystem::phaseName(phaseIdx)
//                                    <<"^"<<FluidSystem::componentName(compIdx)<<"="
//                                    << fluidState.moleFraction(phaseIdx, compIdx) <<"\n";
//                        }
//                    }
//                }
//
//				// interfacial area check (interfacial area between fluid as well as solid phases)
//				for(int phaseIdxI=0; phaseIdxI<numPhases+1; phaseIdxI++){
//					const Scalar eps = 1e-6 ;
//					for (int phaseIdxII=0; phaseIdxII< numPhases+1; ++ phaseIdxII){
//						if (phaseIdxI == phaseIdxII)
//							continue;
//						assert(numEnergyEqs == 3) ; // otherwise this ia call does not make sense
//						const Scalar ia = elemVolVars[scvIdx].interfacialArea(phaseIdxI, phaseIdxII);
//						if (not std::isfinite(ia) or ia < 0.-eps ) {
//							message <<"\nUnphysical Value in interfacial area: \n";
//							message << "\tia" <<FluidSystem::phaseName(phaseIdxI)
//											 <<FluidSystem::phaseName(phaseIdxII)<<"="
//									<< ia << "\n" ;
//							message << "\t S[0]=" << fluidState.saturation(0);
//							message << "\t S[1]=" << fluidState.saturation(1);
//							message << "\t p[0]=" << fluidState.pressure(0);
//							message << "\t p[1]=" << fluidState.pressure(1);
//						}
//					}
//				}
//
//                // General Check
//                for(int phaseIdx=0; phaseIdx<numPhases; phaseIdx++){
//                    const Scalar eps = 1e-6 ;
//                    const Scalar saturationTest = fluidState.saturation(phaseIdx);
//                    if (not std::isfinite(saturationTest) or  saturationTest< 0.-eps or saturationTest > 1.+eps ){
//                        message <<"\nUnphysical Value in Saturation: \n";
//                        message << "\tS" <<"_"<<FluidSystem::phaseName(phaseIdx)<<"=" << std::scientific
//                        << fluidState.saturation(phaseIdx) << std::fixed << "\n";
//                    }
//                }
//
//				// velocity Check
//                const unsigned int globalVertexIdx = this->problem_().vertexMapper().map(*eIt, scvIdx, dim);
//				for(int phaseIdx=0; phaseIdx<numPhases; phaseIdx++){
//					const Scalar eps = 1e-6 ;
//					const Scalar velocityTest = volumeDarcyMagVelocity(phaseIdx, globalVertexIdx);
//					if (not std::isfinite(velocityTest) ){
//						message <<"\nUnphysical Value in Velocity: \n";
//						message << "\tv" <<"_"<<FluidSystem::phaseName(phaseIdx)<<"=" << std::scientific
//						<< velocityTest << std::fixed << "\n";
//					}
//				}
//
//
//                // Some check wrote into the error-message, add some additional information and throw
//                if (not message.str().empty()){
//                    // Getting the spatial coordinate
//                    const GlobalPosition & globalPosCurrent = fvGeometry.subContVol[scvIdx].global;
//                    std::stringstream positionString ;
//
//                    // Add physical location
//                    positionString << "Here:";
//                    for(int i=0; i<dim; i++)
//                        positionString << " x"<< (i+1) << "="  << globalPosCurrent[i] << " "   ;
//                    message << "Unphysical value found! \n" ;
//                    message << positionString.str() ;
//                    message << "\n";
//
//                    message << " Here come the primary Variables:" << "\n" ;
//                    for(unsigned int priVarIdx =0 ; priVarIdx<numEq; ++priVarIdx){
//                        message << "priVar[" << priVarIdx << "]=" << elemVolVars[scvIdx].priVar(priVarIdx) << "\n";
//                    }
//                    DUNE_THROW(NumericalProblem, message.str());
//                }
//            } // end scv-loop
//        } // end element loop
//    }

    /*!
     * \brief Access to the averaged (magnitude of) velocity for each vertex.
     *
     * \param phaseIdx The index of the fluid phase
     * \param dofIdxGlobal The global index of the degree of freedom
     *
     */
    const Scalar volumeDarcyMagVelocity(const unsigned int phaseIdx,
                                        const unsigned int dofIdxGlobal) const
    { return volumeDarcyMagVelocity_[phaseIdx][dofIdxGlobal]; }

    /*!
     * \brief Access to the averaged velocity for each vertex.
     *
     * \param phaseIdx The index of the fluid phase
     * \param dofIdxGlobal The global index of the degree of freedom
     */
    const GlobalPosition volumeDarcyVelocity(const unsigned int phaseIdx,
                                        const unsigned int dofIdxGlobal) const
    { return volumeDarcyVelocity_[phaseIdx][dofIdxGlobal]; }

private:
    PhaseGlobalPositionField volumeDarcyVelocity_;
    PhaseVector         volumeDarcyMagVelocity_ ;
    ScalarVector        boxSurface_ ;
};

} // namespace Dumux

#include "mpncpropertydefaultskinetic.hh"

#endif
