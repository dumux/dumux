/*****************************************************************************
 *   Copyright (C) 2009-2011 by Katherina Baber, Klaus Mosthaf               *
 *   Copyright (C) 2008-2009 by Bernd Flemisch, Andreas Lauser               *
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
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the stokes box model.
 */

#ifndef DUMUX_STOKES_LOCAL_RESIDUAL_BASE_HH
#define DUMUX_STOKES_LOCAL_RESIDUAL_BASE_HH

#include <dumux/boxmodels/common/boxmodel.hh>
//#include <dumux/boxmodels/common/boxcouplinglocalresidual.hh>

#include "stokesproperties.hh"
#include "stokesvolumevariables.hh"
#include "stokesfluxvariables.hh"

#include <dune/common/collectivecommunication.hh>
#include <vector>
#include <iostream>

#include<dune/grid/common/grid.hh>

namespace Dumux
{
/*!
 * \ingroup BoxStokesModel
 * \ingroup BoxLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the stokes box model.
 *
 * This class is also used for the non-isothermal and the stokes transport
 * model, which means that it uses static polymorphism.
 */
template<class TypeTag>
class StokesLocalResidual : public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
protected:
    typedef typename GET_PROP_TYPE(TypeTag, BaseLocalResidual) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, StokesIndices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum {
        dim = GridView::dimension,
        numEq = GET_PROP_VALUE(TypeTag, NumEq)
    };
    enum {
        massBalanceIdx = Indices::massBalanceIdx, //!< Index of the mass balance

        momentumXIdx = Indices::momentumXIdx, //!< Index of the x-component of the momentum balance
        momentumYIdx = Indices::momentumYIdx, //!< Index of the y-component of the momentum balance
        momentumZIdx = Indices::momentumZIdx, //!< Index of the z-component of the momentum balance
        lastMomentumIdx = Indices::lastMomentumIdx //!< Index of the last component of the momentum balance
    };
    enum {
        dimXIdx = Indices::dimXIdx, //!< Index for the first component of a vector
        dimYIdx = Indices::dimYIdx, //!< Index for the second component of a vector
        dimZIdx = Indices::dimZIdx //!< Index for the third component of a vector
    };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::GenericReferenceElements<Scalar, dim> ReferenceElements;
    typedef Dune::GenericReferenceElement<Scalar, dim> ReferenceElement;

    typedef Dune::FieldVector<CoordScalar, dim> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> ScalarGradient;
    typedef Dune::FieldVector<Scalar, dim> FieldVector;

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    typedef typename GridView::Intersection Intersection;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

public:
    /*!
     * \brief Constructor. Sets the upwind weight.
     */
    StokesLocalResidual()
    {
        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        massUpwindWeight_ = GET_PARAM(TypeTag, Scalar, MassUpwindWeight);
        stabilizationAlpha_ = GET_PARAM(TypeTag, Scalar, StabilizationAlpha);
        stabilizationBeta_ = GET_PARAM(TypeTag, Scalar, StabilizationBeta);
    };

    /*!
     * \brief Evaluate the amount all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     *  \param result The mass of the component within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    void computeStorage(PrimaryVariables &result, int scvIdx, bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit Euler method.
        const ElementVolumeVariables &elemVolVars = usePrevSol ? this->prevVolVars_()
                : this->curVolVars_();
        const VolumeVariables &vertexData = elemVolVars[scvIdx];

        result = 0.0;

        // mass balance
        result[massBalanceIdx] = vertexData.density();

        // momentum balance
        for (int momentumIdx = momentumXIdx; momentumIdx <= lastMomentumIdx; ++momentumIdx)
            result[momentumIdx] = vertexData.density()
                * vertexData.velocity()[momentumIdx-momentumXIdx];
    }

    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a sub-control volume.
     *
     * \param flux The flux over the SCV (sub-control-volume) face
     * \param faceIdx The index of the SCV face
     * \param onBoundary Indicates, if the flux is evaluated on a boundary face
     */
    void computeFlux(PrimaryVariables &flux, int faceIdx, bool onBoundary=false) const
    {
        const FluxVariables fluxVars(this->problem_(),
                      this->elem_(),
                      this->fvElemGeom_(),
                      faceIdx,
                      this->curVolVars_(),
                      onBoundary);
        flux = 0.0;

        asImp_()->computeAdvectiveFlux(flux, fluxVars);
        Valgrind::CheckDefined(flux);
        asImp_()->computeDiffusiveFlux(flux, fluxVars);
        Valgrind::CheckDefined(flux);
    }

    /*!
     * \brief Evaluates the advective fluxes over
     *        a face of a subcontrol volume.
     *
     * \param flux The advective flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current SCV
     */
    void computeAdvectiveFlux(PrimaryVariables &flux,
                              const FluxVariables &fluxVars) const
    {
        // if the momentum balance has a dirichlet b.c., the mass balance
        // is replaced, thus we do not need to calculate outflow fluxes here
        if (fluxVars.onBoundary())
            if (momentumBalanceDirichlet_(this->bcTypes_(fluxVars.upstreamIdx())))
                return;

        // data attached to upstream and the downstream vertices
        const VolumeVariables &up = this->curVolVars_(fluxVars.upstreamIdx());
        const VolumeVariables &dn = this->curVolVars_(fluxVars.downstreamIdx());

        // mass balance with upwinded density
        FieldVector massBalanceResidual = fluxVars.velocityAtIP();
        if (massUpwindWeight_ == 1.0) // fully upwind
            massBalanceResidual *= up.density();
        else
            massBalanceResidual *= massUpwindWeight_ * up.density() +
                    (1.-massUpwindWeight_) * dn.density();

        if (!fluxVars.onBoundary())
        {
            // stabilization of the mass balance
            // with 0.5*alpha*(V_i + V_j)*grad P
            FieldVector stabilizationTerm = fluxVars.pressureGradAtIP();
            stabilizationTerm *= stabilizationAlpha_*
                    fluxVars.averageSCVVolume();
            massBalanceResidual += stabilizationTerm;
        }

        flux[massBalanceIdx] +=
                massBalanceResidual*fluxVars.face().normal;

        // momentum balance - pressure is evaluated as volume term
        // at the center of the SCV in computeSource
        // viscosity is upwinded

        // compute symmetrized gradient for the momentum flux:
        // mu (grad v + (grad v)^t)
        Dune::FieldMatrix<Scalar, dim, dim> symmVelGrad = fluxVars.velocityGradAtIP();
        for (int i=0; i<dim; ++i)
            for (int j=0; j<dim; ++j)
                    symmVelGrad[i][j] += fluxVars.velocityGradAtIP()[j][i];

        FieldVector velGradComp(0.);
        for (int velIdx = 0; velIdx < dim; ++velIdx)
        {
            velGradComp = symmVelGrad[velIdx];

            // TODO: dilatation term has to be accounted for in outflow, coupling, neumann
//            velGradComp[velIdx] += 2./3*fluxVars.velocityDivAtIP;

            if (massUpwindWeight_ == 1.0) // fully upwind
                velGradComp *= up.viscosity();
            else
                velGradComp *= massUpwindWeight_ * up.viscosity() +
                        (1.0 - massUpwindWeight_) * dn.viscosity();

            flux[momentumXIdx + velIdx] -=
                    velGradComp*fluxVars.face().normal;

// gravity is accounted for in computeSource; alternatively:
//            Scalar gravityTerm = fluxVars.densityAtIP *
//                    this->problem_().gravity()[dim-1] *
//                    fluxVars.face().ipGlobal[dim-1]*
//                    fluxVars.face().normal[velIdx];
//            flux[momentumXIdx + velIdx] -=
//                    gravityTerm;

        }
    }

    /*!
     * \brief Adds the diffusive flux to the flux vector over
     *        the face of a sub-control volume.
     *
     * It doesn't do anything in the Stokes model but is used by the
     * transport and non-isothermal models to calculate diffusive and
     * conductive fluxes
     *
     * \param flux The diffusive flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current sub control volume face
     */
    void computeDiffusiveFlux(PrimaryVariables &flux,
                              const FluxVariables &fluxVars) const
    { }

    /*!
     * \brief Calculate the source term of the equation,
     *        compute the pressure gradient at the center of a SCV
     *        and evaluate the gravity term
     *
     * \param q The source/sink in the sub control volume for each component
     * \param localVertexIdx The index of the sub control volume
     */
    void computeSource(PrimaryVariables &q, int localVertexIdx)
    {
        const ElementVolumeVariables &elemVolVars = this->curVolVars_();
        const VolumeVariables &vertexData = elemVolVars[localVertexIdx];

        // retrieve the source term intrinsic to the problem
        this->problem_().source(q,
                                this->elem_(),
                                this->fvElemGeom_(),
                                localVertexIdx);

        // ATTENTION: The source term of the mass balance has to be chosen as
        // div (q_momentum) in the problem file
        const Scalar alphaH2 = stabilizationAlpha_*
                this->fvElemGeom_().subContVol[localVertexIdx].volume;
        q[massBalanceIdx] *= alphaH2; // stabilization of the source term

        // pressure gradient at the center of the SCV,
        // the pressure is discretized as volume term,
        // while -mu grad v is calculated in computeFlux
        ScalarGradient pressureGradAtSCVCenter(0.0);
        ScalarGradient grad(0.0);

        for (int vertexIdx = 0; vertexIdx < this->fvElemGeom_().numVertices; vertexIdx++)
        {
            grad = this->fvElemGeom_().subContVol[localVertexIdx].gradCenter[vertexIdx];
            Valgrind::CheckDefined(grad);
            grad *= elemVolVars[vertexIdx].pressure();

            pressureGradAtSCVCenter += grad;
        }

        // add the component of the pressure gradient to the respective part
        // of the momentum equation and take the gravity term into account
        // signs are inverted, since q is subtracted
        for (int dimIdx=0; dimIdx<dim; ++dimIdx)
        {
            q[momentumXIdx + dimIdx] -= pressureGradAtSCVCenter[dimIdx];
            q[momentumXIdx + dimIdx] += vertexData.density()*this->problem_().gravity()[dimIdx];
        }
     }

    // the stokes model needs a modified treatment of the BCs
    void evalBoundary_()
    {
        assert(this->residual_.size() == this->fvElemGeom_().numVertices);
        const ReferenceElement &refElem = ReferenceElements::general(this->elem_().geometry().type());

        // loop over vertices of the element
        for (int vertexIdx = 0; vertexIdx < this->fvElemGeom_().numVertices; vertexIdx++)
        {
            // consider only SCVs on the boundary
            if (this->fvElemGeom_().subContVol[vertexIdx].inner)
                continue;

            // important at corners of the grid
            FieldVector momentumResidual(0.0);
            FieldVector averagedNormal(0.0);
            int numberOfOuterFaces = 0;
            // evaluate boundary conditions for the intersections of
            // the current element
            const BoundaryTypes &bcTypes = this->bcTypes_(vertexIdx);
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
                    // only evaluate, if we consider the same face vertex as in the outer
                    // loop over the element vertices
                    if (refElem.subEntity(faceIdx, 1, faceVertIdx, dim)
                            != vertexIdx)
                        continue;

                    const int boundaryFaceIdx = this->fvElemGeom_().boundaryFaceIndex(faceIdx, faceVertIdx);
                    const FluxVariables boundaryVars(this->problem_(),
                                                     this->elem_(),
                                                     this->fvElemGeom_(),
                                                     boundaryFaceIdx,
                                                     this->curVolVars_(),
                                                     true);

                    // the computed residual of the momentum equations is stored
                    // into momentumResidual for the replacement of the mass balance
                    // in case of Dirichlet conditions for the momentum balance;
                    // the fluxes at the boundary are added in the second step
                    if (momentumBalanceDirichlet_(bcTypes))
                    {
                        FieldVector muGradVelNormal(0.);
                        const FieldVector &boundaryFaceNormal =
                                boundaryVars.face().normal;

                        boundaryVars.velocityGradAtIP().umv(boundaryFaceNormal, muGradVelNormal);
                        muGradVelNormal *= boundaryVars.viscosityAtIP();

                        for (int i=0; i < this->residual_.size(); i++)
                            Valgrind::CheckDefined(this->residual_[i]);
                        for (int dimIdx=0; dimIdx < dim; ++dimIdx)
                            momentumResidual[dimIdx] = this->residual_[vertexIdx][momentumXIdx+dimIdx];

                        //Sign is right!!!: boundary flux: -mu grad v n
                        //but to compensate outernormal -> residual - (-mu grad v n)
                        momentumResidual += muGradVelNormal;
                        averagedNormal += boundaryFaceNormal;
                    }

                    // evaluate fluxes at a single boundary segment
                    asImp_()->evalNeumannSegment_(isIt, vertexIdx, boundaryFaceIdx, boundaryVars);
                    asImp_()->evalOutflowSegment_(isIt, vertexIdx, boundaryFaceIdx, boundaryVars);

                    // count the number of outer faces to determine, if we are on
                    // a corner point and if an interpolation should be done
                    numberOfOuterFaces++;
                } // end loop over face vertices
            } // end loop over intersections

            // replace defect at the corner points of the grid
            // by the interpolation of the primary variables
            if(!bcTypes.isDirichlet(massBalanceIdx))
            {
                if (momentumBalanceDirichlet_(bcTypes))
                    replaceMassbalanceResidual_(momentumResidual, averagedNormal, vertexIdx);
                else // de-stabilize (remove alpha*grad p - alpha div f
                    // from computeFlux on the boundary)
                    removeStabilizationAtBoundary_(vertexIdx);
            }
            if (numberOfOuterFaces == 2)
                interpolateCornerPoints_(bcTypes, vertexIdx);
        } // end loop over element vertices

        // evaluate the dirichlet conditions of the element
        if (this->bcTypes_().hasDirichlet())
            asImp_()->evalDirichlet_();
    }

protected:
    /*!
     * \brief Add Neumann boundary conditions for a single sub-control
     *        volume face to the local residual.
     */
    void evalNeumannSegment_(const IntersectionIterator &isIt,
                             const int scvIdx,
                             const int boundaryFaceIdx,
                             const FluxVariables &boundaryVars)
    {
        const BoundaryTypes &bcTypes = this->bcTypes_(scvIdx);

        // deal with neumann boundaries
        if (bcTypes.hasNeumann())
        {
            // call evalNeumannSegment_() of the base class first
            ParentType::evalNeumannSegment_(isIt, scvIdx, boundaryFaceIdx);

            // temporary vector to store the neumann boundary fluxes
            PrimaryVariables values(0.0);
            if (momentumBalanceHasNeumann_(bcTypes))
            {
                // Neumann BC of momentum equation needs special treatment
                // mathematically Neumann BC: p n - mu grad v n = q
                // boundary terms: -mu grad v n
                // implement q * A (from evalBoundarySegment) - p n(unity) A
                FieldVector pressureCorrection(boundaryVars.face().normal);
                pressureCorrection *= this->curVolVars_(scvIdx).pressure();
                for (int momentumIdx = momentumXIdx; momentumIdx <= lastMomentumIdx; momentumIdx++)
                	if(bcTypes.isNeumann(momentumIdx))
                   		this->residual_[scvIdx][momentumIdx] += pressureCorrection[momentumIdx];

                // beta-stabilization at the boundary
                // in case of neumann conditions for the momentum equation;
                // calculate  mu grad v t t
                // center in the face of the reference element
                FieldVector tangent;
                const FieldVector& elementUnitNormal = isIt->centerUnitOuterNormal();
                tangent[0] = elementUnitNormal[1];  //TODO: 3D
                tangent[1] = -elementUnitNormal[0];
                FieldVector tangentialVelGrad;
                boundaryVars.velocityGradAtIP().mv(tangent, tangentialVelGrad);
                tangentialVelGrad *= boundaryVars.viscosityAtIP();

                this->residual_[scvIdx][massBalanceIdx] -= stabilizationBeta_*0.5*
                                                    this->curVolVars_(scvIdx).pressure();
                this->residual_[scvIdx][massBalanceIdx] -= stabilizationBeta_*0.5*
                                                    (tangentialVelGrad*tangent);

                for (int momentumIdx = momentumXIdx; momentumIdx <= lastMomentumIdx; momentumIdx++)
                    this->residual_[scvIdx][massBalanceIdx] -= stabilizationBeta_*0.5
                        * values[momentumIdx]*elementUnitNormal[momentumIdx-momentumXIdx];
                Valgrind::CheckDefined(this->residual_);
            }
        }
    }

    void evalOutflowSegment_(const IntersectionIterator &isIt,
                             const int scvIdx,
                             const int boundaryFaceIdx,
                             const FluxVariables &boundaryVars)
    {
        // temporary vector to store the neumann boundary fluxes
        const BoundaryTypes &bcTypes = this->bcTypes_(scvIdx);

        // deal with outflow boundaries
        if (bcTypes.hasOutflow())
        {
            PrimaryVariables values(0.0);

            asImp_()->computeFlux(values, boundaryFaceIdx, /*onBoundary=*/true);
            Valgrind::CheckDefined(values);

            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            {
                if (!bcTypes.isOutflow(eqIdx) )
                    continue;
                // do not calculate outflow for the mass balance
                // if the momentum balance is dirichlet -
                // it is replaced in that case
                if (eqIdx==massBalanceIdx && momentumBalanceDirichlet_(bcTypes))
                    continue;
                // deduce outflow
                this->residual_[scvIdx][eqIdx] += values[eqIdx];
            }

            // beta-stabilization at the boundary in case of outflow condition
            // for the momentum balance
            if(momentumBalanceOutflow_(bcTypes) && stabilizationBeta_ != 0)
            {
                // calculate  mu grad v t t for beta-stabilization
                // center in the face of the reference element
                FieldVector tangent;
                const FieldVector& elementUnitNormal = isIt->centerUnitOuterNormal();
                tangent[0] = elementUnitNormal[1];
                tangent[1] = -elementUnitNormal[0];
                FieldVector tangentialVelGrad;
                boundaryVars.velocityGradAtIP().mv(tangent, tangentialVelGrad);

                this->residual_[scvIdx][massBalanceIdx] -= 0.5*stabilizationBeta_
                                                        * boundaryVars.viscosityAtIP()
                                                        * (tangentialVelGrad*tangent);
            }
        }
    }

    void removeStabilizationAtBoundary_(const int vertexIdx)
    {
        // loop over the edges of the element
        for (int faceIdx = 0; faceIdx < this->fvElemGeom_().numEdges; faceIdx++)
        {
            const FluxVariables fluxVars(this->problem_(),
                                         this->elem_(),
                                         this->fvElemGeom_(),
                                         faceIdx,
                                         this->curVolVars_());

            const int i = this->fvElemGeom_().subContVolFace[faceIdx].i;
            const int j = this->fvElemGeom_().subContVolFace[faceIdx].j;

            if (i != vertexIdx && j != vertexIdx)
                continue;

            const Scalar alphaH2 = stabilizationAlpha_*
                    fluxVars.averageSCVVolume();
            Scalar stabilizationTerm = fluxVars.pressureGradAtIP() *
                  this->fvElemGeom_().subContVolFace[faceIdx].normal;

            stabilizationTerm *= alphaH2;

            if (vertexIdx == i)
                this->residual_[i][massBalanceIdx] += stabilizationTerm;
            if (vertexIdx == j)
                this->residual_[j][massBalanceIdx] -= stabilizationTerm;
        }

        //destabilize source term
        PrimaryVariables q(0.0);
        this->problem_().source(q,
                                this->elem_(),
                                this->fvElemGeom_(),
                                vertexIdx);
        const Scalar alphaH2 = stabilizationAlpha_*this->fvElemGeom_().subContVol[vertexIdx].volume;
        this->residual_[vertexIdx][massBalanceIdx] += alphaH2*q[massBalanceIdx]*
                                              this->fvElemGeom_().subContVol[vertexIdx].volume;
    }

    void interpolateCornerPoints_(const BoundaryTypes &bcTypes, const int vertexIdx)
    {
        // TODO: 3D
        if (bcTypes.isCouplingInflow(massBalanceIdx) || bcTypes.isCouplingOutflow(massBalanceIdx))
        {
            if (vertexIdx == 0 || vertexIdx == 3)
                this->residual_[vertexIdx][massBalanceIdx] =
                        this->curPrimaryVars_(0)[massBalanceIdx]-this->curPrimaryVars_(3)[massBalanceIdx];
            if (vertexIdx == 1 || vertexIdx == 2)
                this->residual_[vertexIdx][massBalanceIdx] =
                        this->curPrimaryVars_(1)[massBalanceIdx]-this->curPrimaryVars_(2)[massBalanceIdx];
        }
        else
        {
            if (!bcTypes.isDirichlet(massBalanceIdx)) // do nothing in case of dirichlet
                this->residual_[vertexIdx][massBalanceIdx] =
                        this->curPrimaryVars_(0)[massBalanceIdx]+this->curPrimaryVars_(3)[massBalanceIdx]-
                        this->curPrimaryVars_(1)[massBalanceIdx]-this->curPrimaryVars_(2)[massBalanceIdx];
        }
    }

    void replaceMassbalanceResidual_(const FieldVector& momentumResidual,
                                     FieldVector& averagedNormal,
                                     const int vertexIdx)
    {
        assert(averagedNormal.two_norm() != 0.0);

        // divide averagedNormal by its length
        averagedNormal /= averagedNormal.two_norm();
        // replace the mass balance by the sum of the momentum balances
        this->residual_[vertexIdx][massBalanceIdx] = momentumResidual*averagedNormal;
    }

    // returns true, if all conditions for the momentum balance are dirichlet
    bool momentumBalanceDirichlet_(const BoundaryTypes& bcTypes) const
    {
        for (int momentumIdx=momentumXIdx; momentumIdx<=lastMomentumIdx; ++momentumIdx)
            if (!bcTypes.isDirichlet(momentumIdx))
                return false;
        return true;
    }

    // returns true, if one condition of the momentum balance is neumann
    bool momentumBalanceHasNeumann_(const BoundaryTypes& bcTypes) const
    {
        for (int momentumIdx=momentumXIdx; momentumIdx<=lastMomentumIdx; ++momentumIdx)
            if (bcTypes.isNeumann(momentumIdx))
                return true;
        return false;
    }

    // returns true, if all conditions for the momentum balance are outlow
    bool momentumBalanceOutflow_(const BoundaryTypes& bcTypes) const
    {
        for (int momentumIdx=momentumXIdx; momentumIdx<=lastMomentumIdx; ++momentumIdx)
            if (!bcTypes.isOutflow(momentumIdx))
                return false;
        return true;
    }

    Scalar massUpwindWeight_;
    Scalar stabilizationAlpha_;
    Scalar stabilizationBeta_;

    Implementation *asImp_()
    { return static_cast<Implementation *>(this); }
    const Implementation *asImp_() const
    { return static_cast<const Implementation *>(this); }
};

}

#endif
