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
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the stokes box model.
 */

#ifndef DUMUX_STOKESNCNI_COUPLING_LOCAL_RESIDUAL_HH
#define DUMUX_STOKESNCNI_COUPLING_LOCAL_RESIDUAL_HH

#include <dumux/freeflow/stokesncni/stokesncnilocalresidual.hh>
#include <dumux/freeflow/stokesncni/stokesncnimodel.hh>

namespace Dumux
{
	/*!
   * \ingroup ImplicitLocalResidual
   * \ingroup TwoPTwoCNIStokesTwoCNIModel
	 * \brief Element-wise calculation of the Jacobian matrix for problems
	 *        using the Stokes box model.
	 *
	 * This class is also used for the stokes transport
	 * model, which means that it uses static polymorphism.
   * 
   * \todo Please doc me more!
   * This file should contain a more detailed description of the coupling conditions.
	 */
	template<class TypeTag>
	class StokesncniCouplingLocalResidual : public StokesncniLocalResidual<TypeTag>
	{
	protected:
		
		typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;
		
		typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
		
		typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
		typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
		typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
		
		
		enum {
			dim = GridView::dimension,
			dimWorld = GridView::dimensionworld,
			numEq = GET_PROP_VALUE(TypeTag, NumEq),
			numComponents = Indices::numComponents
		};
		enum {
			//indices of the equations
			massBalanceIdx = Indices::massBalanceIdx, //!< Index of the mass balance
			
			momentumXIdx = Indices::momentumXIdx, //!< Index of the x-component of the momentum balance
			momentumYIdx = Indices::momentumYIdx, //!< Index of the y-component of the momentum balance
			momentumZIdx = Indices::momentumZIdx, //!< Index of the z-component of the momentum balance
			lastMomentumIdx = Indices::lastMomentumIdx, //!< Index of the last component of the momentum balance
			transportEqIdx = Indices::transportEqIdx,//!< Index of the transport equation
			energyEqIdx = Indices::energyEqIdx
		};
		enum {
			//indices of phase and transported component
			phaseIdx = Indices::phaseIdx,
			transportCompIdx = Indices::transportCompIdx,
			temperatureIdx = Indices::temperatureIdx
		};
		enum {
			dimXIdx = Indices::dimXIdx, //!< Index for the first component of a vector
			dimYIdx = Indices::dimYIdx, //!< Index for the second component of a vector
			dimZIdx = Indices::dimZIdx //!< Index for the third component of a vector
		};
		
		typedef typename GridView::ctype CoordScalar;
		typedef Dune::ReferenceElements<Scalar, dim> ReferenceElements;
		typedef Dune::ReferenceElement<Scalar, dim> ReferenceElement;
		
		typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
		typedef Dune::FieldVector<Scalar, dim> DimVector;
		
		typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
		
		typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
		typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
		
		typedef typename GridView::IntersectionIterator IntersectionIterator;
		typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
		
		static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

	public:
        /*!
         * \brief Modified boundary treatment for the stokes model
         */
		void evalBoundary_()
		{
			assert(this->residual_.size() == this->fvGeometry_().numScv);
			const ReferenceElement &refElement = ReferenceElements::general(this->element_().geometry().type());
			
			// loop over vertices of the element
			for (int idx = 0; idx < this->fvGeometry_().numScv; idx++)
			{
				// consider only SCVs on the boundary
				if (this->fvGeometry_().subContVol[idx].inner)
					continue;
				
				// important at corners of the grid
				DimVector momentumResidual(0.0);
				DimVector averagedNormal(0.0);
				int numberOfOuterFaces = 0;
				// evaluate boundary conditions for the intersections of
				// the current element
				const BoundaryTypes &bcTypes = this->bcTypes_(idx);
				IntersectionIterator isIt = this->gridView_().ibegin(this->element_());
				const IntersectionIterator &endIt = this->gridView_().iend(this->element_());
				for (; isIt != endIt; ++isIt)
				{
					// handle only intersections on the boundary
					if (!isIt->boundary())
						continue;
					
					// assemble the boundary for all vertices of the current face
					const int fIdx = isIt->indexInInside();
					const int numFaceVertices = refElement.size(fIdx, 1, dim);
					
					// loop over the single vertices on the current face
					for (int faceVertexIdx = 0; faceVertexIdx < numFaceVertices; ++faceVertexIdx)
					{
						// only evaluate, if we consider the same face vertex as in the outer
						// loop over the element vertices
						if (refElement.subEntity(fIdx, 1, faceVertexIdx, dim)
                            != idx)
							continue;
						
						const int boundaryFaceIdx = this->fvGeometry_().boundaryFaceIndex(fIdx, faceVertexIdx);
						const FluxVariables boundaryVars(this->problem_(),
														 this->element_(),
														 this->fvGeometry_(),
														 boundaryFaceIdx,
														 this->curVolVars_(),
														 true);
						
						// the computed residual of the momentum equations is stored
						// into momentumResidual for the replacement of the mass balance
						// in case of Dirichlet conditions for the momentum balance;
						// the fluxes at the boundary are added in the second step
						if (this->momentumBalanceDirichlet_(bcTypes))
						{
							DimVector muGradVelNormal(0.);
							const DimVector &boundaryFaceNormal =
							boundaryVars.face().normal;
							
							boundaryVars.velocityGrad().umv(boundaryFaceNormal, muGradVelNormal);
							muGradVelNormal *= (boundaryVars.dynamicViscosity()
                                                + boundaryVars.dynamicEddyViscosity());
							
							for (int i=0; i < this->residual_.size(); i++)
								Valgrind::CheckDefined(this->residual_[i]);
							for (int dimIdx=0; dimIdx < dim; ++dimIdx)
								momentumResidual[dimIdx] = this->residual_[idx][momentumXIdx+dimIdx];
							
							//Sign is right!!!: boundary flux: -mu grad v n
							//but to compensate outernormal -> residual - (-mu grad v n)
							momentumResidual += muGradVelNormal;
							averagedNormal += boundaryFaceNormal;
						}
						
						// evaluate fluxes at a single boundary segment
						asImp_()->evalNeumannSegment_(isIt, idx, boundaryFaceIdx, boundaryVars);
						asImp_()->evalOutflowSegment_(isIt, idx, boundaryFaceIdx, boundaryVars);
						
						//for the corner points, the boundary flux across the vertical non-coupling boundary faces
						//has to be calculated to fulfill the mass balance
						//convert suddomain intersection into multidomain intersection and check whether it is an outer boundary
						if(!GridView::Grid::multiDomainIntersection(*isIt).neighbor()
                           && (this->boundaryHasMortarCoupling_(this->bcTypes_(idx)) || this->momentumBalanceHasNeumann_(this->bcTypes_(idx))))
						{
						    const GlobalPosition& globalPos = this->fvGeometry_().subContVol[idx].global;
                            //problem specific function, in problem orientation of interface is known
                            if(this->problem_().isInterfaceCornerPoint(globalPos))
                            {
                                 PrimaryVariables priVars(0.0);
                                //                         DimVector faceCoord = this->fvGeometry_().boundaryFace[boundaryFaceIdx].ipGlobal;
                                //                         std::cout<<faceCoord<<std::endl;

                                 const int numVertices = refElement.size(dim);
                                 bool evalBoundaryFlux = false;
                                 for(int equationIdx = 0; equationIdx < numEq; ++equationIdx)
                                 {
                                     for(int i= 0; i < numVertices; i++)
                                     {
                                         //if vertex is on boundary and not the coupling vertex: check whether an outflow condition is set
                                         if(this->model_().onBoundary(this->element_(), i) && i!=idx)
                                             if (!this->bcTypes_(i).isOutflow(equationIdx))
                                                 evalBoundaryFlux = true;
                                     }

                                     //calculate the actual boundary fluxes and add to residual (only for momentum and transport equation, mass balance already has outflow)
                                     if(evalBoundaryFlux)
                                     {
                                         asImp_()->computeFlux(priVars, boundaryFaceIdx, true/*on boundary*/);
                                         this->residual_[idx][equationIdx] += priVars[equationIdx];
                                     }
                                 }
                            }
						}
						// Beavers-Joseph condition at the coupling boundary/interface
						if(boundaryHasCoupling_(bcTypes))
						{
							evalBeaversJoseph_(isIt, idx, boundaryFaceIdx, boundaryVars);
						}
						if(boundaryHasCoupling_(bcTypes) || boundaryHasMortarCoupling_(bcTypes))
                        {
                            asImp_()->evalCouplingVertex_(isIt, idx, boundaryFaceIdx, boundaryVars);
                        }
						// count the number of outer faces to determine, if we are on
						// a corner point and if an interpolation should be done
						numberOfOuterFaces++;
					} // end loop over face vertices
				} // end loop over intersections
				
				// replace defect at the corner points of the grid
				// by the interpolation of the primary variables
				if(!bcTypes.isDirichlet(massBalanceIdx))
				{
					if (this->momentumBalanceDirichlet_(bcTypes))
						this->replaceMassbalanceResidual_(momentumResidual, averagedNormal, idx);
					else // de-stabilize (remove alpha*grad p - alpha div f
						// from computeFlux on the boundary)
						this->removeStabilizationAtBoundary_(idx);
				}
				if (numberOfOuterFaces == 2)
					this->interpolateCornerPoints_(bcTypes, idx);
			} // end loop over element vertices
			
			// evaluate the dirichlet conditions of the element
			if (this->bcTypes_().hasDirichlet())
				asImp_()->evalDirichlet_();
		}
		
		void evalBoundaryPDELab_()
		{
			// loop over vertices of the element
			for (int idx = 0; idx < this->fvGeometry_().numScv; idx++)
			{
				// consider only SCVs on the boundary
				if (this->fvGeometry_().subContVol[idx].inner)
					continue;
				
				this->removeStabilizationAtBoundary_(idx);
			} // end loop over vertices
		}
		
	protected:
		/*!
		 * \brief Evaluate one part of the Dirichlet-like coupling conditions for a single
		 *        sub-control volume face; rest is done in the local coupling operator
		 */
		void evalCouplingVertex_(const IntersectionIterator &isIt,
								 const int scvIdx,
								 const int boundaryFaceIdx,
								 const FluxVariables &boundaryVars)
		{
			const VolumeVariables &volVars = this->curVolVars_()[scvIdx];
			const BoundaryTypes &bcTypes = this->bcTypes_(scvIdx);
			
			// TODO: beaversJosephCoeff is used to determine, if we are on a coupling face
			// this is important at the corners. However, a better way should be found.
			
			// check if BJ-coeff is not zero to decide, if coupling condition
			// for the momentum balance (Dirichlet vor v.n) has to be set;
			// may be important at corners
			const Scalar beaversJosephCoeff = this->problem_().beaversJosephCoeff(this->element_(),
																				  this->fvGeometry_(),
																				  *isIt,
																				  scvIdx,
																				  boundaryFaceIdx);
			// set velocity normal to the interface
			if (bcTypes.isCouplingInflow(momentumYIdx) && beaversJosephCoeff)
				this->residual_[scvIdx][momentumYIdx] =
				volVars.velocity() *
				boundaryVars.face().normal /
				boundaryVars.face().normal.two_norm();
			Valgrind::CheckDefined(this->residual_[scvIdx][momentumYIdx]);
			
			// add pressure correction - required for pressure coupling,
			// if p.n comes from the pm
			if ((bcTypes.isCouplingOutflow(momentumYIdx) && beaversJosephCoeff) || bcTypes.isMortarCoupling(momentumYIdx))
			{
				DimVector pressureCorrection(isIt->centerUnitOuterNormal());
				pressureCorrection *= volVars.pressure(); // TODO: 3D
				this->residual_[scvIdx][momentumYIdx] += pressureCorrection[momentumYIdx]*
				boundaryVars.face().area;
			}
			
			// set mole fraction for the transported components
			for (int compIdx = 0; compIdx < numComponents; compIdx++)
			{
				int eqIdx =  dim + compIdx;
				if (eqIdx != massBalanceIdx) {
                    if (bcTypes.isCouplingOutflow(eqIdx))
                    {
                        if(useMoles)
                            this->residual_[scvIdx][eqIdx] = volVars.moleFraction(compIdx);
                        else
                            this->residual_[scvIdx][eqIdx] = volVars.massFraction(compIdx);
                        Valgrind::CheckDefined(this->residual_[scvIdx][eqIdx]);
                    }
                }
			}
			// set temperature
			if (bcTypes.isCouplingOutflow(energyEqIdx))
			{
				this->residual_[scvIdx][energyEqIdx] = volVars.temperature();
				Valgrind::CheckDefined(this->residual_[scvIdx][energyEqIdx]);
			}
		}
		
		void evalBeaversJoseph_(const IntersectionIterator &isIt,
								const int scvIdx,
								const int boundaryFaceIdx,
								const FluxVariables &boundaryVars) //TODO: required
		{
			const BoundaryTypes &bcTypes = this->bcTypes_(scvIdx);
			
			Scalar beaversJosephCoeff = this->problem_().beaversJosephCoeff(this->element_(),
																			this->fvGeometry_(),
																			*isIt,
																			scvIdx,
																			boundaryFaceIdx);
			
			// only enter here, if we are on a coupling boundary (see top)
			// and the BJ coefficient is not zero
			if (beaversJosephCoeff)
			{
				const Scalar Kxx = this->problem_().permeability(this->element_(),
																 this->fvGeometry_(),
																 scvIdx);
				beaversJosephCoeff /= sqrt(Kxx);
				const DimVector& elementUnitNormal = isIt->centerUnitOuterNormal();
				
				// implementation as NEUMANN condition /////////////////////////////////////////////
				// (v.n)n
				if (bcTypes.isCouplingOutflow(momentumXIdx))
				{
					const Scalar normalComp = boundaryVars.velocity()*elementUnitNormal;
					DimVector normalV = elementUnitNormal;
					normalV *= normalComp; // v*n*n
					
					// v_tau = v - (v.n)n
					const DimVector tangentialV = boundaryVars.velocity() - normalV;
					const Scalar boundaryFaceArea = boundaryVars.face().area;
					
					for (int dimIdx=0; dimIdx < dim; ++dimIdx)
						this->residual_[scvIdx][dimIdx] += beaversJosephCoeff *
						boundaryFaceArea *
						tangentialV[dimIdx] *
						(boundaryVars.dynamicViscosity() + boundaryVars.dynamicEddyViscosity());
					
					////////////////////////////////////////////////////////////////////////////////////
					//normal component has only to be set if no coupling conditions are defined
					//set NEUMANN flux (set equal to pressure in problem)
					//                PrimaryVariables priVars(0.0);
					//                this->problem_().neumann(priVars, this->element_(), this->fvGeometry_(),
					//                                  *isIt, scvIdx, boundaryFaceIdx);
					//                for (int dimIdx=0; dimIdx < dim; ++dimIdx)
					//                    this->residual_[scvIdx][dimIdx] += priVars[dimIdx]*
					//                                                       boundaryFaceArea;
				}
				if (bcTypes.isCouplingInflow(momentumXIdx)) //TODO: 3D
				{
					///////////////////////////////////////////////////////////////////////////////////////////
					//IMPLEMENTATION AS DIRICHLET CONDITION
					//tangential componment: vx = sqrt K /alpha * (grad v n(unity))t
					DimVector tangentialVelGrad(0);
					boundaryVars.velocityGrad().umv(elementUnitNormal, tangentialVelGrad);
					tangentialVelGrad /= -beaversJosephCoeff; // was - before
					this->residual_[scvIdx][momentumXIdx] =
					tangentialVelGrad[momentumXIdx] - this->curPriVars_(scvIdx)[momentumXIdx];
					
					////////////////////////////////////////////////////////////////////////////////////
					//for testing Beavers and Joseph without adjacent porous medium set vy to 0
					//                Scalar normalVel(0);
					//                this->residual_[scvIdx][momentumYIdx] = this->curPrimaryVars_(scvIdx)[momentumYIdx] - normalVel;
					////////////////////////////////////////////////////////////////////////////////////
				}
			}
		}
		
		// return true, if at least one equation on the boundary has a  coupling condition
		bool boundaryHasCoupling_(const BoundaryTypes& bcTypes) const
		{
			for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
				if (bcTypes.isCouplingInflow(eqIdx) || bcTypes.isCouplingOutflow(eqIdx))
					return true;
			return false;
		}
		
	    // return true, if at least one equation on the boundary has a mortar coupling condition
	    bool boundaryHasMortarCoupling_(const BoundaryTypes& bcTypes) const
	    {
	        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
	            if (bcTypes.isMortarCoupling(eqIdx))
	                return true;
	        return false;
	    }

	private:
		Implementation *asImp_()
		{ return static_cast<Implementation *>(this); }
		const Implementation *asImp_() const
		{ return static_cast<const Implementation *>(this); }
		
	};
	
} // Dumux

#endif // DUMUX_STOKESNCNI_COUPLING_LOCAL_RESIDUAL_HH
