/*****************************************************************************
 *   Copyright (C) 2009 by Markus Wolff                                      *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#ifndef DUMUX_FVVELOCITY2P_ADAPTIVE_HH
#define DUMUX_FVVELOCITY2P_ADAPTIVE_HH

/**
 * @file
 * @brief  Velocity Field from a finite volume solution of a pressure equation.
 * @author Markus Wolff
 */

#include <dumux/decoupled/2p/diffusion/fv/fvpressure2padaptive.hh>

namespace Dumux
{
//! \ingroup FV2p
//! \brief Determines the velocity from a finite volume solution of the  pressure equation of the sequential Model (IMPES).
/*! Calculates phase velocities or total velocity from a known pressure field in context of a Finite Volume implementation for the evaluation
 * of equations of the form
 * \f[\text{div}\, \boldsymbol{v}_{total} = q.\f]
 * The wetting or the non-wetting phase pressure, or the global pressure has to be given as piecewise constant cell values.
 * The phase velocities are calculated following  Darcy's law as
 * \f[\boldsymbol{v}_\alpha = \lambda_\alpha \boldsymbol{K} \left(\text{grad}\, p_\alpha + \rho_\alpha g  \text{grad}\, z\right),\f]
 * where \f$p_\alpha\f$ denotes the pressure of phase \f$_\alpha\f$ (wetting or non-wetting), \f$\boldsymbol{K}\f$ the absolute permeability, \f$\lambda_\alpha\f$ the phase mobility, \f$\rho_\alpha\f$ the phase density and \f$g\f$ the gravity constant.
 * The total velocity is either calculated as sum of the phase velocities
 * \f[\boldsymbol{v}_{total} = \boldsymbol{v}_{wetting}+\boldsymbol{v}_{non-wetting},\f]
 * or with a given global pressure
 * \f[\boldsymbol{v}_{total} = \lambda_{total} \boldsymbol{K} \left(\text{grad}\, p_{global} + \sum f_\alpha \rho_\alpha g  \text{grad}\, z\right).\f]
 *
 * \tparam TypeTag The Type Tag
 */

template<class TypeTag>
class FVVelocity2Padaptive: public FVPressure2Padaptive<TypeTag>
{
    typedef FVVelocity2Padaptive<TypeTag> ThisType;
    typedef FVPressure2Padaptive<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
     typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
     typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
     typedef typename GET_PROP_TYPE(TypeTag, PTAG(Variables)) Variables;

     typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;
     typedef typename SpatialParameters::MaterialLaw MaterialLaw;

     typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

     typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
     typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

     typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;
     typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes))::PrimaryVariables PrimaryVariables;

typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef Dune::GenericReferenceElements<Scalar, dim> ReferenceElementContainer;
    typedef Dune::GenericReferenceElements<Scalar, dim-1> ReferenceElementFaceContainer;
    typedef Dune::GenericReferenceElement<Scalar, dim> ReferenceElement;

    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureNW,
        pglobal = Indices::pressureGlobal,
        vw = Indices::velocityW,
        vn = Indices::velocityNW,
        vt = Indices::velocityTotal,
        Sw = Indices::saturationW,
        Sn = Indices::saturationNW,
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        eqIdxPress = Indices::pressEqIdx,
        eqIdxSat = Indices::satEqIdx
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
    };

    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

public:
    //! Constructs a FVVelocity2Padaptive object
    /*!
     * \param problem a problem class object
     */
    FVVelocity2Padaptive(Problem& problem)
    : FVPressure2Padaptive<TypeTag>(problem)
    {
    	// todo: kompatibilität prüfen
        if (GET_PROP_VALUE(TypeTag, PTAG(EnableCompressibility)) && velocityType_ == vt)
        {
            DUNE_THROW(Dune::NotImplemented, "Total velocity - global pressure - model cannot be used with compressible fluids!");
        }
        if (velocityType_ != vw && velocityType_ != vn && velocityType_ != vt)
        {
            DUNE_THROW(Dune::NotImplemented, "Velocity type not supported!");
        }
    }

    //! Calculate the velocity.
    /*!
     *
     *  Given the piecewise constant pressure \f$p\f$,
     *  this method calculates the velocity
     *  The method is needed in the IMPES (Implicit Pressure Explicit Saturation) algorithm which is used for a fractional flow formulation
     *  to provide the velocity field required for the solution of the saturation equation.
     */
    void calculateVelocity();

    void update()
    {
        ParentType::update();

        calculateVelocity();

        return;
    }

    //! \brief Write data files
    /*  \param name file name */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        ParentType::addOutputVtkFields(writer);

//        Dune::BlockVector<Dune::FieldVector<Scalar, dim> > &velocity = *(writer.template allocateManagedBuffer<Scalar, dim> (
//                this->problem().gridView().size(0)));
//        Dune::BlockVector<Dune::FieldVector<Scalar, dim> > &velocitySecondPhase = *(writer.template allocateManagedBuffer<Scalar, dim> (
//                this->problem().gridView().size(0)));
//
//        // compute update vector
//        ElementIterator eItEnd = this->problem().gridView().template end<0>();
//        for (ElementIterator eIt = this->problem().gridView().template begin<0>(); eIt != eItEnd; ++eIt)
//        {
//            // cell index
//            int globalIdx = this->problem().variables().index(*eIt);
//
//            Dune::FieldVector<Scalar, 2*dim> fluxW(0);
//            Dune::FieldVector<Scalar, 2*dim> fluxNW(0);
//            // run through all intersections with neighbors and boundary
//            IntersectionIterator
//            isItEnd = this->problem().gridView().iend(*eIt);
//            for (IntersectionIterator
//                    isIt = this->problem().gridView().ibegin(*eIt); isIt
//                    !=isItEnd; ++isIt)
//            {
//                int isIndex = isIt->indexInInside();
//
//                fluxW[isIndex] = isIt->geometry().volume() * (isIt->centerUnitOuterNormal() * this->problem().variables().velocityElementFace(*eIt, isIndex));
//                fluxNW[isIndex] = isIt->geometry().volume() * (isIt->centerUnitOuterNormal() * this->problem().variables().velocitySecondPhase()[this->problem().variables().index(*eIt)][isIndex]);
//            }
//
//            Dune::FieldVector<Scalar, dim> refVelocity(0);
//            refVelocity[0] = 0.5 * (fluxW[1] - fluxW[0]);
//            refVelocity[1] = 0.5 * (fluxW[3] - fluxW[2]);
//
//            const Dune::FieldVector<Scalar, dim>& localPos = ReferenceElementContainer::general(eIt->geometry().type()).position(0,
//                    0);
//
//            // get the transposed Jacobian of the element mapping
//            const FieldMatrix& jacobianInv = eIt->geometry().jacobianInverseTransposed(localPos);
//            FieldMatrix jacobianT(jacobianInv);
//            jacobianT.invert();
//
//            // calculate the element velocity by the Piola transformation
//            Dune::FieldVector<Scalar, dim> elementVelocity(0);
//            jacobianT.umtv(refVelocity, elementVelocity);
//            elementVelocity /= eIt->geometry().integrationElement(localPos);
//
//            velocity[globalIdx] = elementVelocity;
//
//            refVelocity = 0;
//            refVelocity[0] = 0.5 * (fluxNW[1] - fluxNW[0]);
//            refVelocity[1] = 0.5 * (fluxNW[3] - fluxNW[2]);
//
//            // calculate the element velocity by the Piola transformation
//            elementVelocity = 0;
//            jacobianT.umtv(refVelocity, elementVelocity);
//            elementVelocity /= eIt->geometry().integrationElement(localPos);
//
//            velocitySecondPhase[globalIdx] = elementVelocity;
//        }
//
//        //switch velocities
//        switch (velocityType_)
//        {
//        case vw:
//        {
//            writer.attachCellData(velocity, "wetting-velocity", dim);
//            writer.attachCellData(velocitySecondPhase, "non-wetting-velocity", dim);
//            break;
//        }
//        case vn:
//        {
//            writer.attachCellData(velocity, "non-wetting-velocity", dim);
//            writer.attachCellData(velocitySecondPhase, "wetting-velocity", dim);
//            break;
//        }
//        case vt:
//        {
//            writer.attachCellData(velocity, "total velocity", dim);
//            break;
//        }
//        }
//
//        return;
    }

private:
    static const int velocityType_ = GET_PROP_VALUE(TypeTag, PTAG(VelocityFormulation)); //!< gives kind of velocity used (\f$ 0 = v_w\f$, \f$ 1 = v_n\f$, \f$ 2 = v_t\f$)
};

template<class TypeTag>
void FVVelocity2Padaptive<TypeTag>::calculateVelocity()
{
    //reset velocity
    this->problem().variables().velocity() = Dune::FieldVector<Scalar, dimWorld>(0.0);
    this->problem().variables().velocitySecondPhase() = Dune::FieldVector<Scalar, dimWorld>(0.0);

    BoundaryTypes bcType;

    // compute update vector
    ElementIterator eItEnd = this->problem().gridView().template end<0>();
    for (ElementIterator eIt = this->problem().gridView().template begin<0>(); eIt != eItEnd; ++eIt)
    {
#if HAVE_MPI
        if (eIt->partitionType() == Dune::GhostEntity || eIt->partitionType() == Dune::OverlapEntity)
        {
            continue;
        }
#endif

        // cell index
        int globalIdxI = this->problem().variables().index(*eIt);

        GlobalPosition globalPos = eIt->geometry().center();

        Scalar pressI = this->problem().variables().pressure()[globalIdxI];
        Scalar pcI = this->problem().variables().capillaryPressure(globalIdxI);
        Scalar lambdaWI = this->problem().variables().mobilityWetting(globalIdxI);
        Scalar lambdaNWI = this->problem().variables().mobilityNonwetting(globalIdxI);
        Scalar fractionalWI = this->problem().variables().fracFlowFuncWetting(globalIdxI);
        Scalar fractionalNWI = this->problem().variables().fracFlowFuncNonwetting(globalIdxI);
        Scalar densityWI = this->problem().variables().densityWetting(globalIdxI);
        Scalar densityNWI = this->problem().variables().densityNonwetting(globalIdxI);

        // run through all intersections with neighbors and boundary
//        int isIndex = -1;
        IntersectionIterator isItEnd = this->problem().gridView().iend(*eIt);
        for (IntersectionIterator
                isIt = this->problem().gridView().ibegin(*eIt); isIt
                !=isItEnd; ++isIt)
        {
            // local number of facet
//        	isIndex++;

            int isIndex = isIt->indexInInside();

            Dune::FieldVector<Scalar,dimWorld> unitOuterNormal = isIt->centerUnitOuterNormal();

            // handle interior face
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer neighborPointer = isIt->outside();
                int globalIdxJ = this->problem().variables().index(*neighborPointer);

                if (neighborPointer->level()==eIt.level())
                {
					// neighbor cell center in global coordinates
					const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().center();

					// distance vector between barycenters
					Dune::FieldVector<Scalar,dimWorld> distVec = globalPosNeighbor - globalPos;

					// compute distance between cell centers
					Scalar dist = distVec * unitOuterNormal;
//					Scalar dist = distVec.two_norm();

					// compute vectorized permeabilities
					FieldMatrix meanPermeability(0);

					this->problem().spatialParameters().meanK(meanPermeability,
							this->problem().spatialParameters().intrinsicPermeability(*eIt),
							this->problem().spatialParameters().intrinsicPermeability(*neighborPointer));

					Dune::FieldVector<Scalar, dim> permeability(0);
					meanPermeability.mv(unitOuterNormal, permeability);

					Scalar pressJ = this->problem().variables().pressure()[globalIdxJ];
					Scalar pcJ = this->problem().variables().capillaryPressure(globalIdxJ);
					Scalar lambdaWJ = this->problem().variables().mobilityWetting(globalIdxJ);
					Scalar lambdaNWJ = this->problem().variables().mobilityNonwetting(globalIdxJ);
					Scalar fractionalWJ = this->problem().variables().fracFlowFuncWetting(globalIdxJ);
					Scalar fractionalNWJ = this->problem().variables().fracFlowFuncNonwetting(globalIdxJ);
					Scalar densityWJ = this->problem().variables().densityWetting(globalIdxJ);
					Scalar densityNWJ = this->problem().variables().densityNonwetting(globalIdxJ);

					//calculate potential gradients
					Scalar potentialW = 0;
					Scalar potentialNW = 0;

					potentialW = this->problem().variables().potentialWetting(globalIdxI, isIndex);
					potentialNW = this->problem().variables().potentialNonwetting(globalIdxI, isIndex);

					Scalar densityW = (potentialW> 0.) ? densityWI : densityWJ;
					Scalar densityNW = (potentialNW> 0.) ? densityNWI : densityNWJ;

					densityW = (potentialW == 0.) ? 0.5 * (densityWI + densityWJ) : densityW;
					densityNW = (potentialNW == 0.) ? 0.5 * (densityNWI + densityNWJ) : densityNW;

					switch (this->pressureType)
					{
						case pw:
						{
							potentialW = (pressI - pressJ);
							potentialNW = (pressI - pressJ+ pcI - pcJ);
							break;
						}
						case pn:
						{
							potentialW = (pressI - pressJ - pcI + pcJ);
							potentialNW = (pressI - pressJ);
							break;
						}
						case pglobal:
						{
							potentialW = (pressI - pressJ - 0.5 * (fractionalNWI+fractionalNWJ)*(pcI - pcJ));
							potentialNW = (pressI - pressJ + 0.5 * (fractionalWI+fractionalWJ)*(pcI - pcJ));
							break;
						}
					}

					potentialW += densityW * (distVec * this->gravity);//delta z/delta x in unitOuterNormal[z]
					potentialNW += densityNW * (distVec * this->gravity);

					//store potentials for further calculations (velocity, saturation, ...)
					this->problem().variables().potentialWetting(globalIdxI, isIndex) = potentialW;
					this->problem().variables().potentialNonwetting(globalIdxI, isIndex) = potentialNW;

					//do the upwinding of the mobility depending on the phase potentials
					Scalar lambdaW = (potentialW > 0.) ? lambdaWI : lambdaWJ;
					lambdaW = (potentialW == 0) ? 0.5 * (lambdaWI + lambdaWJ) : lambdaW;
					Scalar lambdaNW = (potentialNW > 0.) ? lambdaNWI : lambdaNWJ;
					lambdaNW = (potentialNW == 0) ? 0.5 * (lambdaNWI + lambdaNWJ) : lambdaNW;
					densityW = (potentialW> 0.) ? densityWI : densityWJ;
					densityW = (potentialW == 0.) ? 0.5 * (densityWI + densityWJ) : densityW;
					densityNW = (potentialNW> 0.) ? densityNWI : densityNWJ;
					densityNW = (potentialNW == 0.) ? 0.5 * (densityNWI + densityNWJ) : densityNW;

					//calculate the gravity term
					Dune::FieldVector<Scalar,dimWorld> velocityW(permeability);
					Dune::FieldVector<Scalar,dimWorld> velocityNW(permeability);
					Dune::FieldVector<Scalar,dimWorld> gravityTermW(unitOuterNormal);
					Dune::FieldVector<Scalar,dimWorld> gravityTermNW(unitOuterNormal);

					gravityTermW *= (this->gravity*permeability)*(lambdaW * densityW);
					gravityTermNW *= (this->gravity*permeability)*(lambdaNW * densityNW);

					//calculate velocity depending on the pressure used -> use pc = pn - pw
					switch (this->pressureType)
					{
						case pw:
						{
							velocityW *= lambdaW * (pressI - pressJ)/dist;
							velocityNW *= lambdaNW * (pressI - pressJ)/dist + 0.5 * (lambdaNWI + lambdaNWJ) * (pcI - pcJ) / dist;
							velocityW += gravityTermW;
							velocityNW += gravityTermNW;
							break;
						}
						case pn:
						{
							velocityW *= lambdaW * (pressI - pressJ)/dist - 0.5 * (lambdaWI + lambdaWJ) * (pcI - pcJ) / dist;
							velocityNW *= lambdaNW * (pressI - pressJ) / dist;
							velocityW += gravityTermW;
							velocityNW += gravityTermNW;
							break;
						}
						case pglobal:
						{
							this->problem().variables().velocity()[globalIdxI][isIndex] = permeability;
							this->problem().variables().velocity()[globalIdxI][isIndex]*= (lambdaW + lambdaNW)* (pressI - pressJ )/ dist;
							this->problem().variables().velocity()[globalIdxI][isIndex] += gravityTermW;
							this->problem().variables().velocity()[globalIdxI][isIndex] += gravityTermNW;
							break;
						}
					}

					//store velocities
					switch (velocityType_)
					{
						case vw:
						{
							this->problem().variables().velocity()[globalIdxI][isIndex] = velocityW;
							this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityNW;
							break;
						}
						case vn:
						{
							this->problem().variables().velocity()[globalIdxI][isIndex] = velocityNW;
							this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityW;
							break;
						}
						case vt:
						{
							switch (this->pressureType)
							{
								case pw:
								{
									this->problem().variables().velocity()[globalIdxI][isIndex] = velocityW + velocityNW;
									this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityNW;
									break;
								}
								case pn:
								{
									this->problem().variables().velocity()[globalIdxI][isIndex] = velocityW + velocityNW;
									this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityW;
									break;
								}
							}
							break;
						}
					} //end of switch (velocityType_)
                } //End of case "same level"

                if (neighborPointer->level()==eIt.level()+1)
                {
                	int globalIdxK = 0;
					ElementPointer thirdCellPointer = isIt->outside();
					bool foundK=false;
					bool foundIJ=false;
					// We are looking for two things:
					// IsIndexJ, the index of the interface from the neighbor-cell point of view
					// GlobalIdxK, the index of the third cell
					// for efficienty this is done in one IntersectionIterator-Loop

					Scalar areaWeight = 0;

					// Intersectioniterator around cell J
//					int isIndexJ = -1;
					IntersectionIterator isItEndJ = this->problem().gridView().iend(*neighborPointer);

					for (IntersectionIterator isItJ = this->problem().gridView().ibegin(*neighborPointer); isItJ != isItEndJ; ++isItJ)
					{
						// increase isIndexJJ, if it is not found yet
						if (!foundIJ)
//							isIndexJ++;
						if (isItJ->neighbor())
						{
							ElementPointer neighborPointer2 = isItJ->outside();

							// Neighbor of neighbor is Cell I?
							if (this->problem().variables().index(*neighborPointer2)==globalIdxI)
							{
								foundIJ=true;
								areaWeight = isItJ->geometry().volume();
							}
							else
							{
								if (neighborPointer2->level()==eIt.level()+1)
								{
									// To verify, we found the correct Cell K, we check
									// - is level(K)=level(J)?
									// - is distance(IJ)=distance(IK)?
									// - K is neighbor of J.
									const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().center();
									const GlobalPosition& globalPosThirdCell = neighborPointer2->geometry().center();
									Dune::FieldVector<Scalar, dimWorld> distVecIJ = globalPosNeighbor - globalPos;
									Dune::FieldVector<Scalar, dimWorld> distVecIK = globalPosThirdCell - globalPos;
									if (((distVecIK.two_norm()-distVecIJ.two_norm()))/distVecIJ.two_norm() < 0.001)
									{
										globalIdxK= this->problem().variables().index(*neighborPointer2);
										thirdCellPointer = neighborPointer2;
										foundK=true;
									}

								}
							}
						}
						if (foundIJ && foundK) break;
					}

					int isIndexJ = isIt->indexInOutside();
					areaWeight /= isIt->geometry().volume();
//					areaWeight = 1.0;

					// neighbor cell center in global coordinates
					const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().center();
					const GlobalPosition& globalPosInterface = isIt->geometry().center();

					Dune::FieldVector<Scalar,dimWorld> distVec = globalPosInterface - globalPos;
					Scalar lI= distVec*unitOuterNormal;
					distVec = globalPosNeighbor - globalPosInterface;
					Scalar lJ= distVec*unitOuterNormal;
					Scalar l=lI+lJ;

					FieldMatrix permeabilityI(0);
					FieldMatrix permeabilityJ(0);
					FieldMatrix permeabilityK(0);

					this->problem().spatialParameters().meanK(permeabilityI,
							this->problem().spatialParameters().intrinsicPermeability(*eIt));
					this->problem().spatialParameters().meanK(permeabilityJ,
							this->problem().spatialParameters().intrinsicPermeability(*neighborPointer));
					this->problem().spatialParameters().meanK(permeabilityK,
							this->problem().spatialParameters().intrinsicPermeability(*thirdCellPointer));

					// Calculate permeablity component normal to interface
					Scalar kI, kJ, kK, ng, kMean;//, gI, gJ, gK;
					Dune::FieldVector<Scalar, dim> permI(0);
					Dune::FieldVector<Scalar, dim> permJ(0);
					Dune::FieldVector<Scalar, dim> permK(0);

					permeabilityI.mv(unitOuterNormal, permI);
					permeabilityJ.mv(unitOuterNormal, permJ);
					permeabilityK.mv(unitOuterNormal, permK);

					// kI,kJ,kK = (n^T)Kn
					kI=unitOuterNormal*permI;
					kJ=unitOuterNormal*permJ;
					kK=unitOuterNormal*permK;
					kMean=kI*kJ*kK*l/(kJ*kK*lI+kI*(kJ+kK)/2*lJ);

					ng=this->gravity*unitOuterNormal;

					// find secondary variables in cell J and K
					Scalar pressJ = this->problem().variables().pressure()[globalIdxJ];
					Scalar pcJ = this->problem().variables().capillaryPressure(globalIdxJ);
					Scalar lambdaWJ = this->problem().variables().mobilityWetting(globalIdxJ);
					Scalar lambdaNWJ = this->problem().variables().mobilityNonwetting(globalIdxJ);
					Scalar densityWJ = this->problem().variables().densityWetting(globalIdxJ);
					Scalar densityNWJ = this->problem().variables().densityNonwetting(globalIdxJ);
					Scalar fractionalWJ = this->problem().variables().fracFlowFuncWetting(globalIdxJ);
					Scalar fractionalNWJ = this->problem().variables().fracFlowFuncNonwetting(globalIdxJ);

					Scalar pressK = this->problem().variables().pressure()[globalIdxK];
					Scalar pcK = this->problem().variables().capillaryPressure(globalIdxK);
					Scalar densityWK = this->problem().variables().densityWetting(globalIdxK);
					Scalar densityNWK = this->problem().variables().densityNonwetting(globalIdxK);
					Scalar fractionalWK = this->problem().variables().fracFlowFuncWetting(globalIdxK);
					Scalar fractionalNWK = this->problem().variables().fracFlowFuncNonwetting(globalIdxK);

					Scalar pressJK=(pressJ+pressK)/2;
					Scalar pcJK=(pcJ+pcK)/2;


					// calculate potential gradients
					// reuse potentials from fvpressure2padaptive

					Scalar potentialWIJ = this->problem().variables().potentialWetting(globalIdxI, isIndex);
					Scalar potentialNWIJ = this->problem().variables().potentialNonwetting(globalIdxI, isIndex);
					// We dont know the insideIndex of interface IK
					Scalar potentialWIK = potentialWIJ;
					Scalar potentialNWIK = potentialNWIJ;
					// preliminary potential. The "real" ones are found below

					// Comment: reversed weighting is plausible, too (swap lJ and lI)
					Scalar rhoMeanWIJ = (lJ*densityWI + lI*densityWJ)/l;
					Scalar rhoMeanNWIJ = (lJ*densityNWI + lI*densityNWJ)/l;
					Scalar rhoMeanWIK = (lJ*densityWI + lI*densityWK)/l;
					Scalar rhoMeanNWIK = (lJ*densityNWI + lI*densityNWK)/l;

					Scalar fMeanWIJ = (lJ*fractionalWI+lI*fractionalWJ)/l;
					Scalar fMeanNWIJ = (lJ*fractionalNWI+lI*fractionalNWJ)/l;
					Scalar fMeanWIK = (lJ*fractionalWI+lI*fractionalWK)/l;
					Scalar fMeanNWIK = (lJ*fractionalNWI+lI*fractionalNWK)/l;

					// Upwinding for finding the upwinding direction
					Scalar densityWIJ = (potentialWIJ> 0.) ? densityWI : densityWJ;
					Scalar densityNWIJ = (potentialNWIJ> 0.) ? densityNWI : densityNWJ;
					Scalar densityWIK = (potentialWIK> 0.) ? densityWI : densityWK;
					Scalar densityNWIK = (potentialNWIK> 0.) ? densityNWI : densityNWK;

					densityWIJ = (potentialWIJ == 0.) ? rhoMeanWIJ : densityWIJ;
					densityNWIJ = (potentialNWIJ == 0.) ? rhoMeanNWIJ : densityNWIJ;
					densityWIK = (potentialWIK == 0.) ? rhoMeanWIK : densityWIK;
					densityNWIK = (potentialNWIK == 0.) ? rhoMeanNWIK : densityNWIK;

					Scalar fractionalWIJ = (potentialWIJ> 0.) ? fractionalWI : fractionalWJ;
					Scalar fractionalNWIJ = (potentialNWIJ> 0.) ? fractionalNWI : fractionalNWJ;
					Scalar fractionalWIK = (potentialWIK> 0.) ? fractionalWI : fractionalWK;
					Scalar fractionalNWIK = (potentialNWIK> 0.) ? fractionalNWI : fractionalNWK;

					fractionalWIJ = (potentialWIJ == 0.) ? fMeanWIJ : fractionalWIJ;
					fractionalNWIJ = (potentialNWIJ == 0.) ? fMeanNWIJ : fractionalNWIJ;
					fractionalWIK = (potentialWIK == 0.) ? fMeanWIK : fractionalWIK;
					fractionalNWIK = (potentialNWIK == 0.) ? fMeanNWIK : fractionalNWIK;

					switch (this->pressureType)
					{
						case pw:
						{
							potentialWIJ = (pressI-pressJK)/l+
									(densityWIJ-lJ/l*(kI+kK)/kI*(densityWIK-densityWIJ)/2)*ng;
							potentialNWIJ = (pressI+pcI-(pressJK+pcJK))/l+
									(densityNWIJ-lJ/l*(kI+kK)/kI*(densityNWIK-densityNWIJ)/2)*ng;
							potentialWIK = (pressI-pressJK)/l+
									(densityWIK-lJ/l*(kI+kK)/kI*(densityWIJ-densityWIK)/2)*ng;
							potentialNWIK = (pressI+pcI-(pressJK+pcJK))/l+
									(densityNWIK-lJ/l*(kI+kK)/kI*(densityNWIJ-densityNWIK)/2)*ng;
							break;
						}
						case pn:
						{
							potentialWIJ = (pressI-pcI-(pressJK-pcJK))/l+
									(densityWIJ-lJ/l*(kI+kK)/kI*(densityWIK-densityWIJ)/2)*ng;
							potentialNWIJ = (pressI-pressJK)/l+
									(densityNWIJ-lJ/l*(kI+kK)/kI*(densityNWIK-densityNWIJ)/2)*ng;
							potentialWIK = (pressI-pcI-(pressJK-pcJK))/l+
									(densityWIK-lJ/l*(kI+kK)/kI*(densityWIJ-densityWIK)/2)*ng;
							potentialNWIK = (pressI-pressJK)/l+
									(densityNWIK-lJ/l*(kI+kK)/kI*(densityNWIJ-densityNWIK)/2)*ng;
							break;
						}
						case pglobal:
						{
							potentialWIJ = (pressI-fractionalNWIJ*pcI-(pressJK-fractionalNWIJ*pcJK))/l+
									(densityWIJ-lJ/l*(kI+kK)/kI*(densityWIK-densityWIJ)/2)*ng;
							potentialNWIJ = (pressI+fractionalWIJ*pcI-(pressJK+fractionalWIJ*pcJK))/l+
									(densityNWIJ-lJ/l*(kI+kK)/kI*(densityNWIK-densityNWIJ)/2)*ng;
							potentialWIK = (pressI-fractionalNWIK*pcI-(pressJK-fractionalNWIK*pcJK))/l+
									(densityWIK-lJ/l*(kI+kK)/kI*(densityWIJ-densityWIK)/2)*ng;
							potentialNWIK = (pressI+fractionalWIK*pcI-(pressJK+fractionalWIK*pcJK))/l+
									(densityNWIK-lJ/l*(kI+kK)/kI*(densityNWIJ-densityNWIK)/2)*ng;
							break;
						}
					}

					//store potentials for further calculations (velocity, saturation, ...)
					// these quantities only have correct sign (needed for upwinding)
					// potentials are defined slightly different for adaptive scheme
					this->problem().variables().potentialWetting(globalIdxI, isIndex) = potentialWIJ;
					this->problem().variables().potentialNonwetting(globalIdxI, isIndex) = potentialNWIJ;
					this->problem().variables().potentialWetting(globalIdxJ, isIndexJ) = -potentialWIJ;
					this->problem().variables().potentialNonwetting(globalIdxJ, isIndexJ) = -potentialNWIJ;


					//do the upwinding of the mobility depending on the phase potentials
					Scalar lambdaWIJ = (potentialWIJ > 0.) ? lambdaWI : lambdaWJ;
					lambdaWIJ = (potentialWIJ == 0) ? 0.5 * (lambdaWI + lambdaWJ) : lambdaWIJ;
					Scalar lambdaNWIJ = (potentialNWIJ > 0.) ? lambdaNWI : lambdaNWJ;
					lambdaNWIJ = (potentialNWIJ == 0) ? 0.5 * (lambdaNWI + lambdaNWJ) : lambdaNWIJ;

					densityWIJ = (potentialWIJ> 0.) ? densityWI : densityWJ;
					densityNWIJ = (potentialNWIJ> 0.) ? densityNWI : densityNWJ;
					densityWIK = (potentialWIK> 0.) ? densityWI : densityWK;
					densityNWIK = (potentialNWIK> 0.) ? densityNWI : densityNWK;

					densityWIJ = (potentialWIJ == 0.) ? rhoMeanWIJ : densityWIJ;
					densityNWIJ = (potentialNWIJ == 0.) ? rhoMeanNWIJ : densityNWIJ;
					densityWIK = (potentialWIK == 0.) ? rhoMeanWIK : densityWIK;
					densityNWIK = (potentialNWIK == 0.) ? rhoMeanNWIK : densityNWIK;

					fractionalWIJ = (potentialWIJ> 0.) ? fractionalWI : fractionalWJ;
					fractionalNWIJ = (potentialNWIJ> 0.) ? fractionalNWI : fractionalNWJ;
					fractionalWIK = (potentialWIK> 0.) ? fractionalWI : fractionalWK;
					fractionalNWIK = (potentialNWIK> 0.) ? fractionalNWI : fractionalNWK;

					fractionalWIJ = (potentialWIJ == 0.) ? fMeanWIJ : fractionalWIJ;
					fractionalNWIJ = (potentialNWIJ == 0.) ? fMeanNWIJ : fractionalNWIJ;
					fractionalWIK = (potentialWIK == 0.) ? fMeanWIK : fractionalWIK;
					fractionalNWIK = (potentialNWIK == 0.) ? fMeanNWIK : fractionalNWIK;

					//calculate velocities and the gravity term
					Dune::FieldVector<Scalar,dimWorld> velocityW(unitOuterNormal);
					Dune::FieldVector<Scalar,dimWorld> velocityNW(unitOuterNormal);
					Dune::FieldVector<Scalar,dimWorld> gravityTermW(unitOuterNormal);
					Dune::FieldVector<Scalar,dimWorld> gravityTermNW(unitOuterNormal);

					gravityTermW *= lambdaWIJ*kMean*ng;
					gravityTermW *= densityWIJ-(lJ/l)*(kI+kK)/kI*(densityWIK-densityWIJ)/2;
					gravityTermNW *= lambdaNWIJ*kMean*ng;
					gravityTermNW *= densityNWIJ-(lJ/l)*(kI+kK)/kI*(densityNWIK-densityNWIJ)/2;

					switch (this->pressureType)
					{
						case pw:
						{
							velocityW *= lambdaWIJ*kMean*(pressI-pressJK)/l;
							velocityNW *= lambdaNWIJ*kMean*(pressI+pcI-(pressJK+pcJK))/l;

							velocityW += gravityTermW;
							velocityNW += gravityTermNW;
							break;
						}
						case pn:
						{
							velocityNW *= lambdaNWIJ*kMean*(pressI-pressJK)/l;
							velocityW *= lambdaWIJ*kMean*(pressI-pcI-(pressJK-pcJK))/l;

							velocityW += gravityTermW;
							velocityNW += gravityTermNW;
							break;
						}
						case pglobal:
						{
							velocityW *= lambdaWIJ*kMean*(pressI-fractionalNWIJ*pcI-(pressJK-fractionalNWIJ*pcJK))/l;
							velocityNW *= lambdaNWIJ*kMean*(pressI+fractionalWIJ*pcI-(pressJK+fractionalWIJ*pcJK))/l;

							velocityW += gravityTermW;
							velocityNW += gravityTermNW;
							break;
						}
					}

					switch (velocityType_)
					{
						case vw:
						{
						    Dune::FieldVector<Scalar,dimWorld> vel(velocityW);
						    vel*=areaWeight;
							this->problem().variables().velocity()[globalIdxI][isIndex] += vel;
						    vel = velocityNW;
						    vel*=areaWeight;
							this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] += vel;
//                            this->problem().variables().velocity()[globalIdxI][isIndex] = velocityW;
//                            this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityNW;
							this->problem().variables().velocity()[globalIdxJ][isIndexJ] = velocityW;
							this->problem().variables().velocitySecondPhase()[globalIdxJ][isIndexJ] = velocityNW;
							break;
						}
						case vn:
						{
                            Dune::FieldVector<Scalar,dimWorld> vel(velocityNW);
                            vel*=areaWeight;
							this->problem().variables().velocity()[globalIdxI][isIndex] += (vel);
                            vel = velocityW;
                            vel*=areaWeight;
							this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] += (vel);
//                            this->problem().variables().velocity()[globalIdxI][isIndex] = velocityNW;
//                            this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityW;
							this->problem().variables().velocity()[globalIdxJ][isIndexJ] = velocityNW;
							this->problem().variables().velocitySecondPhase()[globalIdxJ][isIndexJ] = velocityW;
							break;
						}
						case vt:
						{
							switch (this->pressureType)
							{
								case pw:
								{
								    Dune::FieldVector<Scalar,dimWorld> vel(velocityW);
								    vel+=velocityNW;
									this->problem().variables().velocity()[globalIdxI][isIndex] += (vel *=areaWeight);
									vel = velocityNW;
									this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] += (vel*=areaWeight);
//                                    this->problem().variables().velocity()[globalIdxI][isIndex] = (velocityW+velocityNW);
//                                    this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityNW;
									this->problem().variables().velocity()[globalIdxJ][isIndexJ] = velocityW + velocityNW;
									this->problem().variables().velocitySecondPhase()[globalIdxJ][isIndexJ] = velocityNW;

									break;
								}
								case pn:
								{
                                    Dune::FieldVector<Scalar,dimWorld> vel(velocityW);
                                    vel+=velocityNW;
									this->problem().variables().velocity()[globalIdxI][isIndex] += (vel *=areaWeight);
									vel = velocityW;
									this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] += (vel*=areaWeight);
//								                                      this->problem().variables().velocity()[globalIdxI][isIndex] = (velocityW+velocityNW);
//								                                      this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityW;
									this->problem().variables().velocity()[globalIdxJ][isIndexJ] = velocityW + velocityNW;
									this->problem().variables().velocitySecondPhase()[globalIdxJ][isIndexJ] = velocityW;
									break;
								}
							}
							break;
						}
					}
                }
            }//end intersection with neighbor

            // handle boundary face
            if (isIt->boundary())
            {
                //get boundary type
                this->problem().boundaryTypes(bcType, *isIt);
                PrimaryVariables boundValues(0.0);

                if (bcType.isDirichlet(eqIdxPress))
                {
                    this->problem().dirichlet(boundValues, *isIt);

                    // center of face in global coordinates
                    GlobalPosition globalPosFace = isIt->geometry().center();

                    // distance vector between barycenters
                    Dune::FieldVector<Scalar,dimWorld> distVec = globalPosFace - globalPos;

                    // compute distance between cell centers
                    Scalar dist = distVec * unitOuterNormal;
//                    Scalar dist = distVec.two_norm();

                    //permeability vector at boundary
                    // compute vectorized permeabilities
                    FieldMatrix meanPermeability(0);

                    this->problem().spatialParameters().meanK(meanPermeability,
                            this->problem().spatialParameters().intrinsicPermeability(*eIt));

                    Dune::FieldVector<Scalar, dim> permeability(0);
                    meanPermeability.mv(unitOuterNormal, permeability);

                    Scalar satBound = 0;
                    if (bcType.isDirichlet(eqIdxSat))
                    {
                        satBound = boundValues[saturationIdx];
                    }
                    else
                    {
                        satBound = this->problem().variables().saturation()[globalIdxI];
                    }

                    //determine phase saturations from primary saturation variable
                    Scalar satW;
                    //Scalar satNW;
                    switch (this->saturationType)
                    {
                    case Sw:
                    {
                        satW = satBound;
                        //satNW = 1-satBound;
                        break;
                    }
                    case Sn:
                    {
                        satW = 1-satBound;
                        //satNW = satBound;
                        break;
                    }
                    default:
                    {
                        DUNE_THROW(Dune::RangeError, "saturation type not implemented");
                        break;
                    }
                    }
                    Scalar pressBound = boundValues[pressureIdx];
                    Scalar pcBound = MaterialLaw::pC(this->problem().spatialParameters().materialLawParams(*eIt), satW);

                    //determine phase pressures from primary pressure variable
                    Scalar pressW=0;
                    Scalar pressNW=0;
                    switch (this->pressureType)
                    {
                    case pw:
                    {
                        pressW = pressBound;
                        pressNW = pressBound + pcBound;
                        break;
                    }
                    case pn:
                    {
                        pressW = pressBound - pcBound;
                        pressNW = pressBound;
                        break;
                    }
                    }

                    //get temperature at current position
                    Scalar temperature = this->problem().temperature(*eIt);

                    Scalar densityWBound = 0;
                    Scalar densityNWBound = 0;
                    Scalar lambdaWBound = 0;
                    Scalar lambdaNWBound = 0;
                    if (this->compressibility)
                    {
                        FluidState fluidState;
                        fluidState.update(satW, pressW, pressNW, temperature);
                        densityWBound = FluidSystem::phaseDensity(wPhaseIdx, temperature, pressW, fluidState);
                        densityNWBound = FluidSystem::phaseDensity(nPhaseIdx, temperature, pressNW, fluidState);
                        Scalar viscosityWBound = FluidSystem::phaseViscosity(wPhaseIdx, temperature, pressW, fluidState);
                        Scalar viscosityNWBound = FluidSystem::phaseViscosity(nPhaseIdx, temperature, pressNW, fluidState);
                        lambdaWBound = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(*eIt), satW) / viscosityWBound * densityWBound;
                        lambdaNWBound = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(*eIt), satW) / viscosityNWBound * densityNWBound;
                    }
                    else
                    {
                        Scalar referencePressure =  this->problem().referencePressure(*eIt);
                        FluidState fluidState;
                        fluidState.update(satW, referencePressure, referencePressure, temperature);
                        densityWBound = FluidSystem::phaseDensity(wPhaseIdx, temperature, referencePressure, fluidState);
                        densityNWBound = FluidSystem::phaseDensity(nPhaseIdx, temperature, referencePressure, fluidState);
                        Scalar viscosityWBound = FluidSystem::phaseViscosity(wPhaseIdx, temperature, referencePressure, fluidState);
                        Scalar viscosityNWBound = FluidSystem::phaseViscosity(nPhaseIdx, temperature, referencePressure, fluidState);
                        lambdaWBound = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(*eIt), satW) / viscosityWBound;
                        lambdaNWBound = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(*eIt), satW) / viscosityNWBound;
                    }

                    Scalar potentialW = 0;
                    Scalar potentialNW = 0;

                    potentialW = this->problem().variables().potentialWetting(globalIdxI, isIndex);
                    potentialNW = this->problem().variables().potentialNonwetting(globalIdxI, isIndex);

                    Scalar densityW = (potentialW> 0.) ? densityWI : densityWBound;
                    Scalar densityNW = (potentialNW> 0.) ? densityNWI : densityNWBound;

                    densityW = (potentialW == 0.) ? 0.5 * (densityWI + densityWBound) : densityW;
                    densityNW = (potentialNW == 0.) ? 0.5 * (densityNWI + densityNWBound) : densityNW;

                    //calculate potential gradient
                    switch (this->pressureType)
                    {
                    case pw:
                    {
                        potentialW = (pressI - pressBound);
                        potentialNW = (pressI + pcI - pressBound - pcBound);
                        break;
                    }
                    case pn:
                    {
                        potentialW = (pressI - pcI - pressBound + pcBound);
                        potentialNW = (pressI - pressBound);
                        break;
                    }
                    case pglobal:
                    {
                        potentialW = (pressI - pressBound - fractionalNWI * (pcI - pcBound));
                        potentialNW = (pressI - pressBound + fractionalWI * (pcI - pcBound));
                        break;
                    }
                    }

                    potentialW += densityW * (distVec * this->gravity);
                    potentialNW += densityNW * (distVec * this->gravity);

                    //store potential gradients for further calculations
                    this->problem().variables().potentialWetting(globalIdxI, isIndex) = potentialW;
                    this->problem().variables().potentialNonwetting(globalIdxI, isIndex) = potentialNW;

                    //do the upwinding of the mobility depending on the phase potentials
                    Scalar lambdaW = (potentialW > 0.) ? lambdaWI : lambdaWBound;
                    lambdaW = (potentialW == 0) ? 0.5 * (lambdaWI + lambdaWBound) : lambdaW;
                    Scalar lambdaNW = (potentialNW > 0.) ? lambdaNWI : lambdaNWBound;
                    lambdaNW = (potentialNW == 0) ? 0.5 * (lambdaNWI + lambdaNWBound) : lambdaNW;
                    densityW = (potentialW> 0.) ? densityWI : densityWBound;
                    densityW = (potentialW == 0.) ? 0.5 * (densityWI + densityWBound) : densityW;
                    densityNW = (potentialNW> 0.) ? densityNWI : densityNWBound;
                    densityNW = (potentialNW == 0.) ? 0.5 * (densityNWI + densityNWBound) : densityNW;

                    //calculate the gravity term
                    Dune::FieldVector<Scalar,dimWorld> velocityW(permeability);
                    Dune::FieldVector<Scalar,dimWorld> velocityNW(permeability);
                    Dune::FieldVector<Scalar,dimWorld> gravityTermW(unitOuterNormal);
                    Dune::FieldVector<Scalar,dimWorld> gravityTermNW(unitOuterNormal);

                    gravityTermW *= (this->gravity*permeability)*(lambdaW * densityW);
                    gravityTermNW *= (this->gravity*permeability)*(lambdaNW * densityNW);

                    //calculate velocity depending on the pressure used -> use pc = pn - pw
                    switch (this->pressureType)
                    {
                    case pw:
                    {
                        velocityW *= lambdaW * (pressI - pressBound)/dist;
                        velocityNW *= lambdaNW * (pressI - pressBound)/dist + 0.5 * (lambdaNWI + lambdaNWBound) * (pcI - pcBound) / dist;
                        velocityW += gravityTermW;
                        velocityNW += gravityTermNW;
                        break;
                    }
                    case pn:
                    {
                        velocityW *= lambdaW * (pressI - pressBound)/dist - 0.5 * (lambdaWI + lambdaWBound) * (pcI - pcBound) / dist;
                        velocityNW *= lambdaNW * (pressI - pressBound) / dist;
                        velocityW += gravityTermW;
                        velocityNW += gravityTermNW;
                        break;
                    }
                    case pglobal:
                    {
                        this->problem().variables().velocity()[globalIdxI][isIndex] = permeability;
                        this->problem().variables().velocity()[globalIdxI][isIndex] *= (lambdaW + lambdaNW)* (pressI - pressBound )/ dist;
                        this->problem().variables().velocity()[globalIdxI][isIndex] += gravityTermW;
                        this->problem().variables().velocity()[globalIdxI][isIndex] += gravityTermNW;
                        break;
                    }
                    }

                    //store velocities
                    switch (velocityType_)
                    {
                    case vw:
                    {
                        this->problem().variables().velocity()[globalIdxI][isIndex] = velocityW;
                        this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityNW;
                        break;
                    }
                    case vn:
                    {
                        this->problem().variables().velocity()[globalIdxI][isIndex] = velocityNW;
                        this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityW;
                        break;
                    }
                    case vt:
                    {
                        switch (this->pressureType)
                        {
                        case pw:
                        {
                            this->problem().variables().velocity()[globalIdxI][isIndex] = velocityW + velocityNW;
                            this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityNW;
                            break;
                        }
                        case pn:
                        {
                            this->problem().variables().velocity()[globalIdxI][isIndex] = velocityW + velocityNW;
                            this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityW;
                            break;
                        }
                        }
                        break;
                    }
                    }
                }//end dirichlet boundary

                else if (bcType.isNeumann(eqIdxPress))
                {
                    this->problem().neumann(boundValues, *isIt);

                    Dune::FieldVector<Scalar,dimWorld> velocityW(unitOuterNormal);
                    Dune::FieldVector<Scalar,dimWorld> velocityNW(unitOuterNormal);

                    velocityW *= boundValues[wPhaseIdx];
                    velocityNW *= boundValues[nPhaseIdx];

                    if (!this->compressibility)
                    {
                        velocityW /= densityWI;
                        velocityNW /= densityNWI;
                    }

                    switch (velocityType_)
                    {
                    case vw:
                    {
                        this->problem().variables().velocity()[globalIdxI][isIndex] = velocityW;
                        this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityNW;
                        break;
                    }
                    case vn:
                    {
                        this->problem().variables().velocity()[globalIdxI][isIndex] = velocityNW;
                        this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityW;
                        break;
                    }
                    case vt:
                    {
                        switch (this->pressureType)
                        {
                        case pw:
                        {
                            this->problem().variables().velocity()[globalIdxI][isIndex] = velocityW + velocityNW;
                            this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityNW;
                            break;
                        }
                        case pn:
                        {
                            this->problem().variables().velocity()[globalIdxI][isIndex] = velocityW + velocityNW;
                            this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityW;
                            break;
                        }
                        }
                        break;
                    }
                    }
                }//end neumann boundary
                else
                {
                    DUNE_THROW(Dune::NotImplemented, "No valid boundary condition type defined for pressure equation!");
                }
            }//end boundary
        }// end all intersections
    }// end grid traversal
//                        printvector(std::cout, this->problem().variables().velocity(), "velocity", "row", 4, 1, 3);
    return;
}
}
#endif
