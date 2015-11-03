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
 *        using the Stokes box model.
 */

#ifndef DUMUX_STOKES_LOCAL_RESIDUAL_BASE_HH
#define DUMUX_STOKES_LOCAL_RESIDUAL_BASE_HH

#include <dune/common/float_cmp.hh>
#include <dune/grid/common/grid.hh>

#include <dumux/implicit/common/implicitmodel.hh>
#include "stokesproperties.hh"
#include "stokesvolumevariables.hh"
#include "stokesfluxvariables.hh"

namespace Dumux
{
/*!
 * \ingroup BoxStokesModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the local Jacobian matrix for problems
 *        using the Stokes box model.
 *
 * This class is also used for the non-isothermal and the two-component Stokes
 * model (static polymorphism).
 */
template<class TypeTag>
class StokesLocalResidual : public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
protected:
    typedef typename GET_PROP_TYPE(TypeTag, BaseLocalResidual) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum {
        dim = GridView::dimension,
        numEq = GET_PROP_VALUE(TypeTag, NumEq)
    };
    enum {
        massBalanceIdx = Indices::massBalanceIdx, //!< Index of the mass balance
        momentumXIdx = Indices::momentumXIdx, //!< Index of the x-component of the momentum balance
        lastMomentumIdx = Indices::lastMomentumIdx //!< Index of the last component of the momentum balance
    };
    enum { pressureIdx = Indices::pressureIdx }; //!< Index of the pressure in a solution vector

    typedef typename Dune::ReferenceElements<Scalar, dim> ReferenceElements;
    typedef typename Dune::ReferenceElement<Scalar, dim> ReferenceElement;

    typedef Dune::FieldVector<Scalar, dim> DimVector;

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    static const bool enableUnsymmetrizedVelocityGradient = GET_PROP_VALUE(TypeTag, EnableUnsymmetrizedVelocityGradient);
    static const bool calculateNavierStokes = GET_PROP_VALUE(TypeTag, EnableNavierStokes);
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

 public:
    /*!
     * \brief Constructor. Sets the upwind weight and the stabilization parameters.
     */
    StokesLocalResidual()
    {
        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        massUpwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);
        stabilizationAlpha_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Stokes, StabilizationAlpha);
        stabilizationBeta_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Stokes, StabilizationBeta);
    }

    /*!
     * \brief Evaluates the amount of all conservation quantities
     *        (mass and momentum) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     *  \param storage The mass of the component within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    void computeStorage(PrimaryVariables &storage, const int scvIdx, const bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit Euler method.
        const ElementVolumeVariables &elemVolVars = usePrevSol ? this->prevVolVars_()
            : this->curVolVars_();
        const VolumeVariables &volVars = elemVolVars[scvIdx];

        storage = 0.0;

        if(useMoles)
            // mass balance mole fraction based
            storage[massBalanceIdx] = volVars.molarDensity();
        else
            // mass balance mass fraction based
            storage[massBalanceIdx] = volVars.density();

        // momentum balance
        for (int momentumIdx = momentumXIdx; momentumIdx <= lastMomentumIdx; ++momentumIdx)
            storage[momentumIdx] = volVars.density()
                * volVars.velocity()[momentumIdx-momentumXIdx];
    }

    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a sub-control volume. The face may be within
     *        an element (SCV face) or on the boundary. The advective and
     *        the diffusive fluxes are computed.
     *
     * \param flux The flux over the SCV (sub-control-volume) face
     * \param fIdx The index of the SCV face (may also be a boundary face)
     * \param onBoundary Indicates, if the flux is evaluated on a boundary face. If it is true,
     *        the created fluxVars object contains boundary variables evaluated at the IP of the
     *        boundary face
     */
    void computeFlux(PrimaryVariables &flux, const int fIdx, const bool onBoundary=false) const
    {
        const FluxVariables fluxVars(this->problem_(),
                                     this->element_(),
                                     this->fvGeometry_(),
                                     fIdx,
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
     *        a face of a sub-control volume.
     *
     * \param flux The advective flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current SCV/boundary face
     */
    void computeAdvectiveFlux(PrimaryVariables &flux,
                              const FluxVariables &fluxVars) const
    {
        // if the momentum balance has a dirichlet b.c., the mass balance
        // is replaced, thus we do not need to calculate outflow fluxes here
        if (fluxVars.onBoundary() &&
            momentumBalanceDirichlet_(this->bcTypes_(fluxVars.upstreamIdx())))
        {
            return;
        }

        // data attached to upstream and the downstream vertices
        const VolumeVariables &up = this->curVolVars_(fluxVars.upstreamIdx());
        const VolumeVariables &dn = this->curVolVars_(fluxVars.downstreamIdx());

        DimVector massBalanceResidual = fluxVars.velocity();

        if(useMoles)
        {
            massBalanceResidual *= (massUpwindWeight_ * up.molarDensity()
                                    + (1.-massUpwindWeight_) * dn.molarDensity());
        }
        else
        {
            massBalanceResidual *= (massUpwindWeight_ * up.density()
                                    + (1.-massUpwindWeight_) * dn.density());
        }

        if (!fluxVars.onBoundary())
        {
            // stabilization of the mass balance
            // with 0.5*alpha*(V_i + V_j)*grad P
            DimVector stabilizationTerm = fluxVars.pressureGrad();
            stabilizationTerm *= stabilizationAlpha_*
                fluxVars.averageSCVVolume();
            massBalanceResidual += stabilizationTerm;
        }

        flux[massBalanceIdx] +=
            massBalanceResidual*fluxVars.face().normal;

        // momentum balance - pressure is evaluated as volume term
        // at the center of the SCV in computeSource
        // dynamic viscosity is upwinded

        Dune::FieldMatrix<Scalar, dim, dim> velGrad = fluxVars.velocityGrad();
        if (enableUnsymmetrizedVelocityGradient)
        {
            // nothing has to be done in this case:
            // grad v
        }
        else
        {
            // compute symmetrized gradient for the momentum flux:
            // grad v + (grad v)^T
            for (int i=0; i<dim; ++i)
                for (int j=0; j<dim; ++j)
                    velGrad[i][j] += fluxVars.velocityGrad()[j][i];
        }

        DimVector velGradComp(0.);
        for (int velIdx = 0; velIdx < dim; ++velIdx)
        {
            velGradComp = velGrad[velIdx];

            // TODO: dilatation term has to be accounted for in outflow, coupling, neumann
            //            velGradComp[velIdx] += 2./3*fluxVars.velocityDiv;

            velGradComp *= fluxVars.dynamicViscosity() + fluxVars.dynamicEddyViscosity();

            flux[momentumXIdx + velIdx] -=
                velGradComp*fluxVars.face().normal;

            // gravity is accounted for in computeSource; alternatively:
            //            Scalar gravityTerm = fluxVars.density *
            //                    this->problem_().gravity()[dim-1] *
            //                    fluxVars.face().ipGlobal[dim-1]*
            //                    fluxVars.face().normal[velIdx];
            //            flux[momentumXIdx + velIdx] -=
            //                    gravityTerm;

        }

        // this term changes the Stokes equation to the Navier-Stokes equation
        // rho v (v*n)
        // rho and first v are upwinded, second v is evaluated at the face
        if (calculateNavierStokes)
        {
            for (int dimIndex = 0; dimIndex < dim; ++dimIndex)
                flux[momentumXIdx + dimIndex] +=
                        up.density() * up.velocity()[dimIndex] * fluxVars.normalVelocity();
        }
    }

    /*!
     * \brief Adds the diffusive flux to the flux vector over
     *        a SCV face or a boundary face.
     *
     * It doesn't do anything in the Stokes model but is used by the
     * transport and non-isothermal models to calculate diffusive and
     * conductive fluxes.
     *
     * \param flux The diffusive flux over the SCV face or boundary face for each component
     * \param fluxVars The flux variables at the current SCV/boundary face
     */
    void computeDiffusiveFlux(PrimaryVariables &flux,
                              const FluxVariables &fluxVars) const
    { }

    /*!
     * \brief Calculate the source term of all equations.
     *        The pressure gradient at the center of a SCV is computed
     *        and the gravity term evaluated.
     *
     * \param source The source/sink in the sub control volume for each component
     * \param scvIdx The local index of the sub-control volume
     */
    void computeSource(PrimaryVariables &source, const int scvIdx)
    {
        const ElementVolumeVariables &elemVolVars = this->curVolVars_();
        const VolumeVariables &volVars = elemVolVars[scvIdx];

        // retrieve the source term intrinsic to the problem
        this->problem_().solDependentSource(source,
                                     this->element_(),
                                     this->fvGeometry_(),
                                     scvIdx,
                                     this->curVolVars_());

        // ATTENTION: The source term of the mass balance has to be chosen as
        // div (q_momentum) in the problem file
        const Scalar alphaH2 = stabilizationAlpha_*
            this->fvGeometry_().subContVol[scvIdx].volume;
        source[massBalanceIdx] *= alphaH2; // stabilization of the source term

        // pressure gradient at the center of the SCV,
        // the pressure is discretized as volume term,
        // while -mu grad v is calculated in computeFlux
        DimVector pressureGradAtSCVCenter(0.0);
        DimVector grad(0.0);

        for (int scvIdx2 = 0; scvIdx2 < this->fvGeometry_().numScv; scvIdx2++)
        {
            grad = this->fvGeometry_().subContVol[scvIdx].gradCenter[scvIdx2];
            Valgrind::CheckDefined(grad);
            grad *= elemVolVars[scvIdx2].pressure();

            pressureGradAtSCVCenter += grad;
        }

        // add the component of the pressure gradient to the respective part
        // of the momentum equation and take the gravity term into account
        // signs are inverted, since q is subtracted
        for (int dimIdx = 0; dimIdx < dim; ++dimIdx)
        {
            source[momentumXIdx + dimIdx] -= pressureGradAtSCVCenter[dimIdx];
            source[momentumXIdx + dimIdx] += volVars.density()*this->problem_().gravity()[dimIdx];
        }
    }

    /*!
     * \brief The Stokes model needs a modified treatment of the boundary conditions as
     *        the common box models
     */
    void evalBoundary_()
    {
        assert(this->residual_.size() == this->fvGeometry_().numScv);
        const ReferenceElement &refElement = ReferenceElements::general(this->element_().geometry().type());

        // loop over sub-control volums of the element
        for (int scvIdx = 0; scvIdx < this->fvGeometry_().numScv; scvIdx++)
        {
            // consider only SCVs on the boundary
            if (this->fvGeometry_().subContVol[scvIdx].inner)
                continue;

            // important at corners of the grid
            DimVector momentumResidual(0.0);
            DimVector averagedNormal(0.0);
            int numberOfOuterFaces = 0;
            // evaluate boundary conditions for the intersections of
            // the current element
            const BoundaryTypes &bcTypes = this->bcTypes_(scvIdx);
            for (const auto& intersection : Dune::intersections(this->gridView_(), this->element_()))
            {
                // handle only intersections on the boundary
                if (!intersection.boundary())
                    continue;

                // assemble the boundary for all vertices of the current face
                const int fIdx = intersection.indexInInside();
                const int numFaceVertices = refElement.size(fIdx, 1, dim);

                // loop over the single vertices on the current face
                for (int faceVertexIdx = 0; faceVertexIdx < numFaceVertices; ++faceVertexIdx)
                {
                    // only evaluate, if we consider the same face vertex as in the outer
                    // loop over the element vertices
                    if (refElement.subEntity(fIdx, 1, faceVertexIdx, dim)
                        != scvIdx)
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
                    if (momentumBalanceDirichlet_(bcTypes))
                    {
                        const DimVector &boundaryFaceNormal =
                            boundaryVars.face().normal;

                        Dune::FieldMatrix<Scalar, dim, dim> velGrad = boundaryVars.velocityGrad();

                        DimVector muGradVelNormal(0.0);
                        velGrad.umv(boundaryFaceNormal, muGradVelNormal);
                        muGradVelNormal *= boundaryVars.dynamicViscosity() + boundaryVars.dynamicEddyViscosity();

                        for (unsigned int i=0; i < this->residual_.size(); i++)
                            Valgrind::CheckDefined(this->residual_[i]);
                        for (int dimIdx = 0; dimIdx < dim; ++dimIdx)
                            momentumResidual[dimIdx] = this->residual_[scvIdx][momentumXIdx+dimIdx];

                        //Sign is right!!!: boundary flux: -mu grad v n
                        //but to compensate outernormal -> residual - (-mu grad v n)
                        momentumResidual += muGradVelNormal;
                        averagedNormal += boundaryFaceNormal;
                    }

                    // evaluate fluxes at a single boundary segment
                    asImp_()->evalNeumannSegment_(&intersection, scvIdx, boundaryFaceIdx, boundaryVars);
                    asImp_()->evalOutflowSegment_(&intersection, scvIdx, boundaryFaceIdx, boundaryVars);

                    // count the number of outer faces to determine, if we are on
                    // a corner point and if an interpolation should be done
                    numberOfOuterFaces++;
                } // end loop over face vertices
            } // end loop over intersections

            if(!bcTypes.isDirichlet(massBalanceIdx))
            {
                // replace the mass balance by the sum of the residua of the momentum balance
                if (momentumBalanceDirichlet_(bcTypes))
                    replaceMassbalanceResidual_(momentumResidual, averagedNormal, scvIdx);
                else // de-stabilize (remove alpha*grad p - alpha div f
                     // from computeFlux on the boundary)
                    removeStabilizationAtBoundary_(scvIdx);
            }
            if (numberOfOuterFaces == 2)
                interpolateCornerPoints_(bcTypes, scvIdx);
        } // end loop over element vertices

        // evaluate the dirichlet conditions of the element
        if (this->bcTypes_().hasDirichlet())
            asImp_()->evalDirichlet_();
    }

protected:
    /*!
     * \brief Evaluate and add Neumann boundary conditions for a single sub-control
     *        volume face to the local residual.
     */
    template <class IntersectionIterator>
    void evalNeumannSegment_(const IntersectionIterator &isIt,
                             const int scvIdx,
                             const int boundaryFaceIdx,
                             const FluxVariables &boundaryVars)
    {
        const BoundaryTypes &bcTypes = this->bcTypes_(scvIdx);

        if (bcTypes.hasNeumann())
        {
            // call evalNeumannSegment_() of the base class first
            ParentType::evalNeumannSegment_(isIt, scvIdx, boundaryFaceIdx);

            // temporary vector to store the Neumann boundary fluxes
            PrimaryVariables values(0.0);
            if (momentumBalanceHasNeumann_(bcTypes))
            {
                // Neumann BC of momentum equation needs special treatment
                // mathematically Neumann BC: p n - mu grad v n = q
                // boundary terms: -mu grad v n
                // implement q * A (from evalBoundarySegment) - p n(unity) A
                DimVector pressureCorrection(boundaryVars.face().normal);
                pressureCorrection *= this->curVolVars_(scvIdx).pressure();
                for (int momentumIdx = momentumXIdx; momentumIdx <= lastMomentumIdx; momentumIdx++)
                    if(bcTypes.isNeumann(momentumIdx))
                        this->residual_[scvIdx][momentumIdx] += pressureCorrection[momentumIdx];

                // beta-stabilization at the boundary
                // in case of neumann conditions for the momentum equation;
                // calculate  mu grad v t t
                // center in the face of the reference element
                DimVector tangent;
                if (Dune::FloatCmp::ne<Scalar, Dune::FloatCmp::absolute>(stabilizationBeta_, 0.0, 1.0e-30))
                {
                    if (dim == 3)
                        DUNE_THROW(Dune::NotImplemented, "The beta-stabilization is only implemented for 2D.");

                    const DimVector& elementUnitNormal = isIt->centerUnitOuterNormal();

                    tangent[0] = elementUnitNormal[1];  //TODO: 3D
                    tangent[1] = -elementUnitNormal[0];

                    Dune::FieldMatrix<Scalar, dim, dim> velGrad = boundaryVars.velocityGrad();

                    DimVector muGradVelTangential(0.0);
                    velGrad.mv(tangent, muGradVelTangential);
                    muGradVelTangential *= boundaryVars.dynamicViscosity() + boundaryVars.dynamicEddyViscosity();

                    this->residual_[scvIdx][massBalanceIdx] -= stabilizationBeta_*0.5*
                        this->curVolVars_(scvIdx).pressure();
                    this->residual_[scvIdx][massBalanceIdx] -= stabilizationBeta_*0.5*
                        (muGradVelTangential*tangent);

                    for (int momentumIdx = momentumXIdx; momentumIdx <= lastMomentumIdx; momentumIdx++)
                        this->residual_[scvIdx][massBalanceIdx] -= stabilizationBeta_*0.5
                            * values[momentumIdx]*elementUnitNormal[momentumIdx-momentumXIdx];
                }
                Valgrind::CheckDefined(this->residual_);
            }
        }
    }

    /*!
     * \brief Evaluate outflow boundary conditions for a single SCV face on the boundary.
     */
    template <class IntersectionIterator>
    void evalOutflowSegment_(const IntersectionIterator &isIt,
                             const int scvIdx,
                             const int boundaryFaceIdx,
                             const FluxVariables &boundaryVars)
    {
        const BoundaryTypes &bcTypes = this->bcTypes_(scvIdx);

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
            if(momentumBalanceOutflow_(bcTypes) && Dune::FloatCmp::ne<Scalar, Dune::FloatCmp::absolute>(stabilizationBeta_, 0.0, 1.0e-30))
            {
                // calculate  mu grad v t t for beta-stabilization
                // center in the face of the reference element
                DimVector tangent;
                const DimVector& elementUnitNormal = isIt->centerUnitOuterNormal();
                tangent[0] = elementUnitNormal[1];
                tangent[1] = -elementUnitNormal[0];

                Dune::FieldMatrix<Scalar, dim, dim> velGrad = boundaryVars.velocityGrad();
                if (enableUnsymmetrizedVelocityGradient)
                {
                    // nothing has to be done in this case:
                    // grad v
                }
                else
                {
                    // compute symmetrized gradient for the momentum flux:
                    // mu (grad v + (grad v)^t)
                    for (int i=0; i<dim; ++i)
                        for (int j=0; j<dim; ++j)
                            velGrad[i][j] += boundaryVars.velocityGrad()[j][i];
                }

                DimVector muGradVelTangential(0.0);
                velGrad.umv(tangent, muGradVelTangential);
                muGradVelTangential *= boundaryVars.dynamicViscosity() + boundaryVars.dynamicEddyViscosity();

                this->residual_[scvIdx][massBalanceIdx] -= 0.5*stabilizationBeta_
                    * (muGradVelTangential*tangent);
            }
        }
    }

    /*!
     * \brief Remove the alpha stabilization at boundaries.
     */
    void removeStabilizationAtBoundary_(const int scvIdx)
    {
        if (Dune::FloatCmp::ne<Scalar, Dune::FloatCmp::absolute>(stabilizationAlpha_, 0.0, 1.0e-30))
        {
            // loop over the edges of the element
            for (int fIdx = 0; fIdx < this->fvGeometry_().numScvf; fIdx++)
            {
                const FluxVariables fluxVars(this->problem_(),
                                             this->element_(),
                                             this->fvGeometry_(),
                                             fIdx,
                                             this->curVolVars_());

                const int i = this->fvGeometry_().subContVolFace[fIdx].i;
                const int j = this->fvGeometry_().subContVolFace[fIdx].j;

                if (i != scvIdx && j != scvIdx)
                    continue;

                const Scalar alphaH2 = stabilizationAlpha_*
                    fluxVars.averageSCVVolume();
                Scalar stabilizationTerm = fluxVars.pressureGrad() *
                    this->fvGeometry_().subContVolFace[fIdx].normal;

                stabilizationTerm *= alphaH2;

                if (scvIdx == i)
                    this->residual_[i][massBalanceIdx] += stabilizationTerm;
                if (scvIdx == j)
                    this->residual_[j][massBalanceIdx] -= stabilizationTerm;
            }

            //destabilize source term
            PrimaryVariables source(0.0);
            this->problem_().solDependentSource(source,
                                     this->element_(),
                                     this->fvGeometry_(),
                                     scvIdx,
                                     this->curVolVars_());
            const Scalar alphaH2 = stabilizationAlpha_ * this->fvGeometry_().subContVol[scvIdx].volume;
            this->residual_[scvIdx][massBalanceIdx] += alphaH2 * source[massBalanceIdx] *
                this->fvGeometry_().subContVol[scvIdx].volume;
        }
    }

    /*!
     * \brief Interpolate the pressure at corner points of the grid, thus taking the degree of freedom there.
     *        This is required due to stability reasons.
     */
    void interpolateCornerPoints_(const BoundaryTypes &bcTypes, const int scvIdx)
    {
        if (dim == 3)
            DUNE_THROW(Dune::NotImplemented, "The function interpolateCornerPoints_() is only implemented for 2D.");

        if (bcTypes.isCouplingInflow(massBalanceIdx) || bcTypes.isCouplingOutflow(massBalanceIdx))
        {
            if (scvIdx == 0 || scvIdx == 3)
                this->residual_[scvIdx][massBalanceIdx] =
                    this->curPriVars_(0)[pressureIdx] - this->curPriVars_(3)[pressureIdx];
            if (scvIdx == 1 || scvIdx == 2)
                this->residual_[scvIdx][massBalanceIdx] =
                    this->curPriVars_(1)[pressureIdx] - this->curPriVars_(2)[pressureIdx];
        }
        else
        {
            if (!bcTypes.isDirichlet(massBalanceIdx)) // do nothing in case of dirichlet
                this->residual_[scvIdx][massBalanceIdx] =
                    this->curPriVars_(0)[pressureIdx] + this->curPriVars_(3)[pressureIdx]-
                    this->curPriVars_(1)[pressureIdx] - this->curPriVars_(2)[pressureIdx];
        }
    }

    /*!
     * \brief Replace the local residual of the mass balance equation by
     *        the sum of the residuals of the momentum balance equation.
     */
    void replaceMassbalanceResidual_(const DimVector& momentumResidual,
                                     DimVector& averagedNormal,
                                     const int scvIdx)
    {
        assert( (Dune::FloatCmp::ne<Scalar, Dune::FloatCmp::absolute>(averagedNormal.two_norm(), 0.0, 1.0e-30)) );

        // divide averagedNormal by its length
        averagedNormal /= averagedNormal.two_norm();
        // replace the mass balance by the sum of the residuals of the momentum balances
        this->residual_[scvIdx][massBalanceIdx] = momentumResidual*averagedNormal;
    }

    /*!
     * \brief Returns true, if all boundary conditions for the momentum balance
     *        at the considered vertex are Dirichlet.
     */
    bool momentumBalanceDirichlet_(const BoundaryTypes& bcTypes) const
    {
        for (int momentumIdx=momentumXIdx; momentumIdx<=lastMomentumIdx; ++momentumIdx)
            if (!bcTypes.isDirichlet(momentumIdx))
                return false;
        return true;
    }

    /*!
     * \brief Returns true, if at least one boundary condition of the momentum balance is Neumann.
     */
    bool momentumBalanceHasNeumann_(const BoundaryTypes& bcTypes) const
    {
        for (int momentumIdx=momentumXIdx; momentumIdx<=lastMomentumIdx; ++momentumIdx)
            if (bcTypes.isNeumann(momentumIdx))
                return true;
        return false;
    }

    /*!
     * \brief Returns true, if all boundary conditions for the momentum balance are outflow.
     */
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
