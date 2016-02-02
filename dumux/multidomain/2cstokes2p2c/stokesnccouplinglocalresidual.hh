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
 *        using the coupled compositional Stokes box model.
 */

#ifndef DUMUX_STOKESNC_COUPLING_LOCAL_RESIDUAL_HH
#define DUMUX_STOKESNC_COUPLING_LOCAL_RESIDUAL_HH

#include <dune/common/deprecated.hh>

#include <dumux/freeflow/stokesnc/localresidual.hh>
#include <dumux/freeflow/stokesnc/model.hh>

namespace Dumux
{
/*!
 * \ingroup ImplicitLocalResidual
 * \ingroup TwoPTwoCStokesTwoCModel
 * \ingroup TwoPTwoCZeroEqTwoCModel
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the coupled compositional Stokes box model.
 *        It is derived from the compositional Stokes box model.
 */
template<class TypeTag>
class StokesncCouplingLocalResidual : public StokesncLocalResidual<TypeTag>
{
    typedef StokesncLocalResidual<TypeTag> ParentType;

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
        transportEqIdx = Indices::transportEqIdx, //!< Index of the transport equation
        conti0EqIdx = Indices::conti0EqIdx
    };
    enum {
        //indices of phase and transported component
        phaseIdx = Indices::phaseIdx,
        transportCompIdx = Indices::transportCompIdx
    };
    enum {
        dimXIdx = Indices::dimXIdx, //!< Index for the first component of a vector
        dimYIdx = Indices::dimYIdx, //!< Index for the second component of a vector
        dimZIdx = Indices::dimZIdx //!< Index for the third component of a vector
    };

    typedef typename GridView::ctype CoordScalar;

    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> DimVector;

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

public:
    /*!
     * \brief Implementation of the boundary evaluation for the Stokes model
     *
     * Evaluate one part of the Dirichlet-like coupling conditions for a single
     * sub-control volume face; rest is done in the local coupling operator
     */
    void evalBoundary_()
    {
        ParentType::evalBoundary_();

        typedef Dune::ReferenceElements<Scalar, dim> ReferenceElements;
        typedef Dune::ReferenceElement<Scalar, dim> ReferenceElement;
        const ReferenceElement &refElement = ReferenceElements::general(this->element_().geometry().type());

        // loop over vertices of the element
        for (int scvIdx = 0; scvIdx < this->fvGeometry_().numScv; scvIdx++)
        {
            // consider only SCVs on the boundary
            if (this->fvGeometry_().subContVol[scvIdx].inner)
                continue;

            const BoundaryTypes &bcTypes = this->bcTypes_(scvIdx);

            // evaluate boundary conditions for the intersections of the current element
            for (const auto& intersection : intersections(this->gridView_(), this->element_()))
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
                    if (refElement.subEntity(fIdx, 1, faceVertexIdx, dim) != scvIdx)
                        continue;

                    const int boundaryFaceIdx = this->fvGeometry_().boundaryFaceIndex(fIdx, faceVertexIdx);
                    const FluxVariables boundaryVars(this->problem_(),
                                                     this->element_(),
                                                     this->fvGeometry_(),
                                                     boundaryFaceIdx,
                                                     this->curVolVars_(),
                                                     true);
                    const VolumeVariables &volVars = this->curVolVars_()[scvIdx];

                    // set velocity normal to the interface
                    if (bcTypes.isCouplingDirichlet(momentumYIdx))
                    {
                        this->residual_[scvIdx][momentumYIdx] = volVars.velocity()
                                                                * boundaryVars.face().normal
                                                                / boundaryVars.face().normal.two_norm();
                        Valgrind::CheckDefined(this->residual_[scvIdx][momentumYIdx]);
                    }

                    // add pressure correction - required for pressure coupling,
                    // if p.n comes from the pm
                    if (bcTypes.isCouplingNeumann(momentumYIdx) || bcTypes.isCouplingMortar(momentumYIdx))
                    {
                        DimVector pressureCorrection(intersection.centerUnitOuterNormal());
                        pressureCorrection *= volVars.pressure();
                        this->residual_[scvIdx][momentumYIdx] += pressureCorrection[momentumYIdx]
                                                                 * boundaryVars.face().area;
                        Valgrind::CheckDefined(this->residual_[scvIdx][momentumYIdx]);
                    }

                    // set mole or mass fraction for the transported components
                    for (int compIdx = 0; compIdx < numComponents; compIdx++)
                    {
                        int eqIdx = conti0EqIdx + compIdx;
                        if (eqIdx != massBalanceIdx)
                        {
                            if (bcTypes.isCouplingDirichlet(eqIdx))
                            {
                                if(useMoles)
                                    this->residual_[scvIdx][eqIdx] = volVars.moleFraction(compIdx);
                                else
                                    this->residual_[scvIdx][eqIdx] = volVars.massFraction(compIdx);
                            }
                        }
                    }
                }
            }
        }
    }

    /*!
     * \brief Removes the stabilization for the Stokes model.
     */
    void evalBoundaryPDELab_()
    {
        // loop over vertices of the element
        for (int scvIdx = 0; scvIdx < this->fvGeometry_().numScv; scvIdx++)
        {
            // consider only SCVs on the boundary
            if (this->fvGeometry_().subContVol[scvIdx].inner)
                continue;

            this->removeStabilizationAtBoundary_(scvIdx);
        }
    }

protected:
    /*!
     * \brief Evaluate one part of the Dirichlet-like coupling conditions for a single
     *        sub-control volume face; rest is done in the local coupling operator
     */
    template <class IntersectionIterator>
    DUNE_DEPRECATED_MSG("evalCouplingVertex_ is deprecated. Its functionality is now included in evalBoundary_.")
    void evalCouplingVertex_(const IntersectionIterator &isIt,
                             const int scvIdx,
                             const int boundaryFaceIdx,
                             const FluxVariables &boundaryVars)
    {
        const VolumeVariables &volVars = this->curVolVars_()[scvIdx];
        const BoundaryTypes &bcTypes = this->bcTypes_(scvIdx);

        // set velocity normal to the interface
        if (bcTypes.isCouplingNeumann(momentumYIdx))
            this->residual_[scvIdx][momentumYIdx] =
                    volVars.velocity() *
                    boundaryVars.face().normal /
                    boundaryVars.face().normal.two_norm();
        Valgrind::CheckDefined(this->residual_[scvIdx][momentumYIdx]);

        // add pressure correction - required for pressure coupling,
        // if p.n comes from the pm
        if (bcTypes.isCouplingDirichlet(momentumYIdx) || bcTypes.isCouplingMortar(momentumYIdx))
        {
            DimVector pressureCorrection(isIt->centerUnitOuterNormal());
            pressureCorrection *= volVars.pressure();
            this->residual_[scvIdx][momentumYIdx] += pressureCorrection[momentumYIdx]*
                    boundaryVars.face().area;
        }

        // set mole or mass fraction for the transported components
        for (int compIdx = 0; compIdx < numComponents; compIdx++)
        {
            int eqIdx = dim + compIdx;
            if (eqIdx != massBalanceIdx) {
                if (bcTypes.isCouplingDirichlet(eqIdx))
                {
                    if(useMoles)
                        this->residual_[scvIdx][eqIdx] = volVars.moleFraction(compIdx);
                    else
                        this->residual_[scvIdx][eqIdx] = volVars.massFraction(compIdx);
                    Valgrind::CheckDefined(this->residual_[scvIdx][eqIdx]);
                }
            }
        }
    }

    template <class IntersectionIterator>
    DUNE_DEPRECATED_MSG("evalBeaversJoseph_ is deprecated. Its functionality is now included in the LOP function evalCoupling21().")
    void evalBeaversJoseph_(const IntersectionIterator &isIt,
                            const int scvIdx,
                            const int boundaryFaceIdx,
                            const FluxVariables &boundaryVars)
    {
        const BoundaryTypes &bcTypes = this->bcTypes_(scvIdx);

        const GlobalPosition &globalPos = this->fvGeometry_().boundaryFace[boundaryFaceIdx].ipGlobal;
        Scalar beaversJosephCoeff = this->problem_().beaversJosephCoeffAtPos(globalPos);

        const Scalar Kxx = this->problem_().permeability(this->element_(),
                                                         this->fvGeometry_(),
                                                         scvIdx);
        beaversJosephCoeff /= std::sqrt(Kxx);
        const DimVector& elementUnitNormal = isIt->centerUnitOuterNormal();

        // implementation as NEUMANN condition /////////////////////////////////////////////
        // (v.n)n
        if (bcTypes.isCouplingDirichlet(momentumXIdx))
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
                                        (boundaryVars.dynamicViscosity()
                                          + boundaryVars.dynamicEddyViscosity());

            ////////////////////////////////////////////////////////////////////////////////////
            //normal component has only to be set if no coupling conditions are defined
            //set NEUMANN flux (set equal to pressure in problem)
//             PrimaryVariables priVars(0.0);
//             this->problem_().neumann(priVars, this->element_(), this->fvGeometry_(),
//                             intersection, scvIdx, boundaryFaceIdx);
//             for (int dimIdx=0; dimIdx < dim; ++dimIdx)
//                 this->residual_[scvIdx][dimIdx] += priVars[dimIdx]*
//                                                    boundaryFaceArea;
        }
        if (bcTypes.isCouplingNeumann(momentumXIdx))
        {
            assert(beaversJosephCoeff > 0);
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

    //! \brief Return true, if at least one equation on the boundary has a  coupling condition
    DUNE_DEPRECATED_MSG("boundaryHasCoupling_() is deprecated")
    bool boundaryHasCoupling_(const BoundaryTypes& bcTypes) const
    {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            if (bcTypes.isCouplingInflow(eqIdx) || bcTypes.isCouplingOutflow(eqIdx))
                return true;
        return false;
    }

    //! \brief Return true, if at least one equation on the boundary has a mortar coupling condition
    DUNE_DEPRECATED_MSG("boundaryHasMortarCoupling_() is deprecated")
    bool boundaryHasMortarCoupling_(const BoundaryTypes& bcTypes) const
    {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            if (bcTypes.isMortarCoupling(eqIdx))
                return true;
        return false;
    }
};

} // namespace Dumux

#endif // DUMUX_STOKESNC_COUPLING_LOCAL_RESIDUAL_HH
