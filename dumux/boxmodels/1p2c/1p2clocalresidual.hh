/*****************************************************************************
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
/*!
 * \file
 *
 * \brief Element-wise calculation the local Jacobian for the single-phase,
 *        two-component model in the BOX scheme.
 */

#ifndef DUMUX_ONEP_TWOC_LOCAL_RESIDUAL_HH
#define DUMUX_ONEP_TWOC_LOCAL_RESIDUAL_HH
#define VELOCITY_OUTPUT 1 //1 turns velocity output on, 0 turns it off

#include <dumux/boxmodels/common/boxmodel.hh>

#include <dumux/boxmodels/1p2c/1p2cproperties.hh>
#include <dumux/boxmodels/1p2c/1p2cvolumevariables.hh>
#include <dumux/boxmodels/1p2c/1p2cfluxvariables.hh>
#include <dumux/boxmodels/1p2c/1p2cboundaryvariables.hh>

#include <dune/common/collectivecommunication.hh>
#include <vector>
#include <iostream>

namespace Dumux
{
/*!
 *
 * \ingroup OnePTwoCBoxModel
 *
 * \brief Calculate the local Jacobian for the single-phase,
 *        two-component model in the BOX scheme.
 *
 *  This class is used to fill the gaps in BoxLocalResidual for the 1P-2C flow.
 */
template<class TypeTag>
class OnePTwoCLocalResidual : public BoxLocalResidual<TypeTag>
{
protected:
    typedef OnePTwoCLocalResidual<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalResidual)) Implementation;
    typedef BoxLocalResidual<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVariables)) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryVariables)) BoundaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePTwoCIndices)) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

    enum
    {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),

        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,
        x1Idx = Indices::x1Idx,

        phaseIdx = Indices::phaseIdx,
        comp0Idx = Indices::comp0Idx,
        comp1Idx = Indices::comp1Idx,

        // indices of the equations
        contiEqIdx = Indices::contiEqIdx,
        transEqIdx = Indices::transEqIdx,
    };

    static const bool useMoles = GET_PROP_VALUE(TypeTag, PTAG(UseMoles));

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef Dune::FieldVector<Scalar, dimWorld> Vector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> Tensor;

public:
    /*!
     * \brief Constructor. Sets the upwind weight.
     */
    OnePTwoCLocalResidual()
    {
        // retrieve the upwind weight for the mobility. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        upwindAlpha_ = GET_PARAM(TypeTag, Scalar, UpwindAlpha);
    };

    /*!
     * \brief Evaluate the amount of all conservation quantities
     *        (e.g. phase mass) within a finite volume.
     *
     *        \param result The mass of the component within the sub-control volume
     *        \param scvIdx The index of the considered face of the sub control volume
     *        \param usePrevSol Evaluate function with solution of current or previous time step
     */
    void computeStorage(PrimaryVariables &result, int scvIdx, bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const ElementVolumeVariables &elemVolVars = usePrevSol ? this->prevVolVars_() : this->curVolVars_();
        const VolumeVariables &volVars = elemVolVars[scvIdx];

        result = 0;
        if(!useMoles)
        {
            // storage term of continuity equation - massfractions
            result[contiEqIdx] +=
                volVars.density()*volVars.porosity();
            //storage term of the transport equation - massfractions
            result[transEqIdx] +=
                volVars.density() * volVars.massFrac(comp1Idx) * volVars.porosity();
        }
        else
        {
            // storage term of continuity equation- molefractions
            //careful: molarDensity changes with moleFrac!
            result[contiEqIdx] += volVars.molarDensity()*volVars.porosity();
            // storage term of the transport equation - molefractions
            result[transEqIdx] +=
                volVars.molarDensity()*volVars.moleFrac(comp1Idx) *
                volVars.porosity();
        }

    }

    /*!
     * \brief Evaluates the mass flux over a face of a subcontrol
     *        volume.
     *
     *        \param flux The flux over the SCV (sub-control-volume) face for each component
     *        \param faceId The index of the considered face of the sub control volume
     */
    void computeFlux(PrimaryVariables &flux, int faceId) const
    {
        flux = 0;
        FluxVariables fluxVars(this->problem_(),
                               this->elem_(),
                               this->fvElemGeom_(),
                               faceId,
                               this->curVolVars_());

        asImp_()->computeAdvectiveFlux(flux, fluxVars);
        asImp_()->computeDiffusiveFlux(flux, fluxVars);
    }

    /*!
         * \brief Evaluates the advective mass flux of all components over
         *        a face of a subcontrol volume.
         *
         * \param flux The advective flux over the sub-control-volume face for each component
         * \param vars The flux variables at the current SCV
         */
    void computeAdvectiveFlux(PrimaryVariables &flux, const FluxVariables &fluxVars) const
    {
        ////////
        // advective fluxes of all components in all phases
        ////////

        // data attached to upstream and the downstream vertices
       // of the current phase
       const VolumeVariables &up =
           this->curVolVars_(fluxVars.upstreamIdx());
       const VolumeVariables &dn =
           this->curVolVars_(fluxVars.downstreamIdx());

        if(!useMoles)
        {
            // total mass flux - massfraction
            //KmvpNormal is the Darcy velocity multiplied with the normal vector, calculated in 1p2cfluxvariables.hh
            flux[contiEqIdx] +=
               fluxVars.KmvpNormal() *
               ((     upwindAlpha_)*up.density()/up.viscosity()
                +
                ((1 - upwindAlpha_)*dn.density()/dn.viscosity()));

            // advective flux of the second component - massfraction
            flux[transEqIdx] +=
               fluxVars.KmvpNormal() *
               ((    upwindAlpha_)*up.density() * up.massFrac(comp1Idx)/up.viscosity()
                +
                (1 - upwindAlpha_)*dn.density()*dn.massFrac(comp1Idx)/dn.viscosity());
        }
        else
        {
            // total mass flux - molefraction
            //KmvpNormal is the Darcy velocity multiplied with the normal vector, calculated in 1p2cfluxvariables.hh
            flux[contiEqIdx] +=
               fluxVars.KmvpNormal() *
               ((     upwindAlpha_)*up.molarDensity()/up.viscosity()
                +
                ((1 - upwindAlpha_)*dn.molarDensity()/dn.viscosity()));

            // advective flux of the second component -molefraction
            flux[transEqIdx] +=
               fluxVars.KmvpNormal() *
               ((    upwindAlpha_)*up.molarDensity() * up.moleFrac(comp1Idx)/up.viscosity()
                +
                (1 - upwindAlpha_)*dn.molarDensity() * dn.moleFrac(comp1Idx)/dn.viscosity());
        }

    }

    /*!
         * \brief Adds the diffusive mass flux of all components over
         *        a face of a subcontrol volume.
         *
         * \param flux The diffusive flux over the sub-control-volume face for each component
         * \param vars The flux variables at the current SCV
         */
    void computeDiffusiveFlux(PrimaryVariables &flux, const FluxVariables &fluxVars) const
    {
        Scalar tmp(0);

        // diffusive flux of second component
        if(!useMoles)
        {
            // diffusive flux of the second component - massfraction
            tmp = - fluxVars.porousDiffCoeff() * fluxVars.densityAtIP()*
              (fluxVars.massFracGrad(comp1Idx) * fluxVars.face().normal);

            flux[transEqIdx] += tmp;// * FluidSystem::molarMass(comp1Idx);
        }
        else
        {
            // diffusive flux of the second component - molefraction
            tmp = - fluxVars.porousDiffCoeff() * fluxVars.molarDensityAtIP()*
               (fluxVars.moleFracGrad(comp1Idx) * fluxVars.face().normal);

            // dispersive flux of second component - molefraction
//            Vector normalDisp;
//            fluxVars.dispersionTensor().mv(fluxVars.face().normal, normalDisp);
//            tmp -= fluxVars.molarDensityAtIP()*
//                (normalDisp * fluxVars.moleFracGrad(comp1Idx));

           flux[transEqIdx] += tmp;
        }
    }
    /*!
     * \brief Calculate the source term of the equation
     *        \param q The source/sink in the SCV for each component
     *        \param localVertexIdx The index of the vertex of the sub control volume
     *
     */
    void computeSource(PrimaryVariables &q, int localVertexIdx)
    {
        this->problem_().source(q,
                                this->elem_(),
                                this->fvElemGeom_(),
                                localVertexIdx);
    }

     void evalBoundary_()
    {
        ParentType::evalBoundary_();

        if (this->bcTypes_().hasOutflow())
       {
           evalOutflow_();
       }

        typedef Dune::GenericReferenceElements<Scalar, dim> ReferenceElements;
        typedef Dune::GenericReferenceElement<Scalar, dim> ReferenceElement;
        const ReferenceElement &refElem = ReferenceElements::general(this->elem_().geometry().type());

        // loop over vertices of the element
        for (int idx = 0; idx < this->fvElemGeom_().numVertices; idx++)
        {
            // consider only SCVs on the boundary
            if (this->fvElemGeom_().subContVol[idx].inner)
                continue;
            int numberOfOuterFaces = 0;

            // evaluate boundary conditions for the intersections of
            // the current element
            IntersectionIterator isIt = this->gridView_().ibegin(this->elem_());
            const IntersectionIterator &endIt = this->gridView_().iend(this->elem_());
            for (; isIt != endIt; ++isIt)
            {
                // handle only intersections on the boundary
                if (!isIt->boundary())
                    continue;

                // assemble the boundary for all vertices of the current face
                const int faceIdx = isIt->indexInInside();
                const int numFaceVertices = refElem.size(faceIdx, 1, dim);

                // loop over the single vertices on the current face
                for (int faceVertIdx = 0; faceVertIdx < numFaceVertices; ++faceVertIdx)
                {
                    const int elemVertIdx = refElem.subEntity(faceIdx, 1, faceVertIdx, dim);
                    // only evaluate, if we consider the same face vertex as in the outer
                    // loop over the element vertices
                    if (elemVertIdx != idx)
                        continue;

                    if (boundaryHasCoupling_(this->bcTypes_(idx)))
                    {
                        const int boundaryFaceIdx = this->fvElemGeom_().boundaryFaceIndex(faceIdx, faceVertIdx);

                        asImp_()->evalCouplingVertex_(isIt, elemVertIdx, boundaryFaceIdx);
                    }

                    // count the number of outer faces to determine, if we are on
                    // a corner point and if an interpolation should be done
                    numberOfOuterFaces++;
                }
            } // end loop over intersections

            if (numberOfOuterFaces == 2 && this->fvElemGeom_().numVertices == 4)
                // do interpolation at corners of the grid
                interpolateCornerPoints_(idx);
        } // end loop over vertices
    }

protected:
     /*!
         * \brief Add all Outflow boundary conditions to the local
         *        residual.
         */
        void evalOutflow_()
        {
            Dune::GeometryType geoType = this->elem_().geometry().type();

            typedef typename Dune::GenericReferenceElements<Scalar, dim> ReferenceElements;
            typedef typename Dune::GenericReferenceElement<Scalar, dim> ReferenceElement;
            const ReferenceElement &refElem = ReferenceElements::general(geoType);

            IntersectionIterator isIt = this->gridView_().ibegin(this->elem_());
            const IntersectionIterator &endIt = this->gridView_().iend(this->elem_());
            for (; isIt != endIt; ++isIt)
            {
                // handle only faces on the boundary
                if (!isIt->boundary())
                    continue;

                // Assemble the boundary for all vertices of the current
                // face
                int faceIdx = isIt->indexInInside();
                int numFaceVerts = refElem.size(faceIdx, 1, dim);
                for (int faceVertIdx = 0;
                     faceVertIdx < numFaceVerts;
                     ++faceVertIdx)
                {
                    int elemVertIdx = refElem.subEntity(faceIdx,
                                                        1,
                                                        faceVertIdx,
                                                        dim);

                    int boundaryFaceIdx =
                        this->fvElemGeom_().boundaryFaceIndex(faceIdx, faceVertIdx);

                    // add the residual of all vertices of the boundary
                    // segment
                    evalOutflowSegment_(isIt,
                                        elemVertIdx,
                                        boundaryFaceIdx);
                }
            }
        }

        /*!
        * \brief Add Outflow boundary conditions for a single sub-control
        *        volume face to the local residual.
        */
       void evalOutflowSegment_(const IntersectionIterator &isIt,
                                int scvIdx,
                                int boundaryFaceIdx)
       {
           // temporary vector to store the neumann boundary fluxes
           PrimaryVariables flux(0.0);
           const BoundaryTypes &bcTypes = this->bcTypes_(scvIdx);

           // deal with neumann boundaries
           if (bcTypes.hasOutflow())
           {
               const BoundaryVariables boundaryVars(this->problem_(),
                                                    this->elem_(),
                                                    this->fvElemGeom_(),
                                                    boundaryFaceIdx,
                                                    this->curVolVars_(),
                                                    scvIdx);

               const VolumeVariables& vertVars = this->curVolVars_()[scvIdx];

               // mass balance
               if (bcTypes.isOutflow(contiEqIdx))
               {
                   if(!useMoles) //use massfractions
                   {
                       flux[contiEqIdx] += boundaryVars.KmvpNormal()*vertVars.density()/vertVars.viscosity();
                   }
                   else //use molefractions
                   {
                       flux[contiEqIdx] += boundaryVars.KmvpNormal()*vertVars.molarDensity()/vertVars.viscosity();
                   }
               }

               // component transport
               if (bcTypes.isOutflow(transEqIdx))
               {
                   if(!useMoles)//use massfractions
                   {
                       // advective flux
                       flux[transEqIdx]+= boundaryVars.KmvpNormal()*vertVars.density()/vertVars.viscosity()
                                        *vertVars.fluidState().massFrac(phaseIdx, comp1Idx);

                       // diffusive flux of comp1 component in phase0
                       Scalar tmp = -boundaryVars.porousDiffCoeff()*boundaryVars.densityAtIP()
                                        *(boundaryVars.massFracGrad(comp1Idx)*boundaryVars.boundaryFace().normal);
                       flux[transEqIdx] += tmp;//* FluidSystem::molarMass(comp1Idx);
                   }
                   else //use molefractions
                   {
                       // advective flux
                       flux[transEqIdx]+= boundaryVars.KmvpNormal()*vertVars.molarDensity()/vertVars.viscosity()
                                       *vertVars.fluidState().moleFrac(phaseIdx, comp1Idx);

                       // diffusive flux of comp1 component in phase0
                       Scalar tmp = -boundaryVars.porousDiffCoeff()*boundaryVars.molarDensityAtIP()
                                        *(boundaryVars.moleFracGrad(comp1Idx)*boundaryVars.boundaryFace().normal);
                       flux[transEqIdx] += tmp;
                   }
               }

               this->residual_[scvIdx] += flux;
           }
       }

     void evalCouplingVertex_(const IntersectionIterator &isIt,
                            const int scvIdx,
                            const int boundaryFaceIdx)
     {
         const VolumeVariables &volVars = this->curVolVars_()[scvIdx];
         // set pressure as part of the momentum coupling
         if (this->bcTypes_(scvIdx).isCouplingOutflow(contiEqIdx))
         {
             this->residual_[scvIdx][contiEqIdx] = volVars.pressure();
         }

         if (this->bcTypes_(scvIdx).isCouplingOutflow(transEqIdx))
         {
             if(!useMoles)
                 this->residual_[scvIdx][transEqIdx] = volVars.fluidState().massFrac(phaseIdx, comp1Idx);
             else
                 this->residual_[scvIdx][transEqIdx] = volVars.fluidState().moleFrac(phaseIdx, comp1Idx);
         }
     }

     /*!
          * \brief Interpolate the corner points of the grid.
          */
         void interpolateCornerPoints_(const int idx)
         {
     //        for (int eqnIdx=0; eqnIdx<numEq; ++eqnIdx)
     //        {
     //            //do not interpolate in case of Beavers-Joseph or Dirichlet condition for coupling interface
     //            if (boundaryHasCoupling_(this->bcTypes_(idx)))
     //                continue;
     //
     //            this->residual_[idx][eqnIdx] = this->curPrimaryVars_(0)[eqnIdx]+this->curPrimaryVars_(3)[eqnIdx]
     //                                                -this->curPrimaryVars_(1)[eqnIdx]-this->curPrimaryVars_(2)[eqnIdx];
     //        }
         }

     /*!
          * \brief Check if one of the boundary conditions is coupling.
          */
     bool boundaryHasCoupling_(const BoundaryTypes& bcTypes) const
     {
         bool hasCoupling = false;
         for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
             if (bcTypes.isCouplingInflow(eqIdx) || bcTypes.isCouplingOutflow(eqIdx))
                 hasCoupling = true;

         return hasCoupling;
     }

    Implementation *asImp_()
    { return static_cast<Implementation *> (this); }
    const Implementation *asImp_() const
    { return static_cast<const Implementation *> (this); }

private:
    Scalar upwindAlpha_;
};

}

#endif
