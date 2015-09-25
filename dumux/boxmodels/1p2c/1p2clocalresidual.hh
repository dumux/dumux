// $Id: 1p2clocalresidual.hh 3784 2010-06-24 13:43:57Z bernd $
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

#include <dumux/boxmodels/common/boxmodel.hh>

#include <dumux/boxmodels/1p2c/1p2cproperties.hh>
#include <dumux/boxmodels/1p2c/1p2cvolumevariables.hh>
#include <dumux/boxmodels/1p2c/1p2cfluxvariables.hh>

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
    typedef BoxLocalResidual<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVariables)) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePTwoCIndices)) Indices;
    enum
    {
        dimWorld = GridView::dimensionworld,

        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,
        x1Idx = Indices::x1Idx,

        // indices of the equations
        contiEqIdx = Indices::contiEqIdx,
        transEqIdx = Indices::transEqIdx,
    };

    static const Scalar upwindAlpha = GET_PROP_VALUE(TypeTag, PTAG(UpwindAlpha));

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef Dune::FieldVector<Scalar, dimWorld> Vector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> Tensor;

public:
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
        const VolumeVariables &volVars =
            usePrevSol ?
            this->prevVolVars_(scvIdx) :
            this->curVolVars_(scvIdx);

        // storage term of continuity equation
        result[contiEqIdx] =
            volVars.density()*volVars.porosity();

        // storage term of the transport equation
        result[transEqIdx] =
            volVars.concentration(1) *
            volVars.porosity();
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

        Vector tmpVec;

        fluxVars.intrinsicPermeability().mv(fluxVars.potentialGrad(), tmpVec);

        // "intrinsic" flux from cell i to cell j
        Scalar normalFlux = - (tmpVec*fluxVars.face().normal);
        const VolumeVariables &up = this->curVolVars_(fluxVars.upstreamIdx(normalFlux));
        const VolumeVariables &dn = this->curVolVars_(fluxVars.downstreamIdx(normalFlux));

        // total mass flux
        flux[contiEqIdx] =
            normalFlux *
            ((     upwindAlpha)*up.density()/up.viscosity()
             +
             ((1 - upwindAlpha)*dn.density()/dn.viscosity()));

        // advective flux of the second component
        flux[transEqIdx] +=
            normalFlux *
            ((    upwindAlpha)*up.concentration(1)/up.viscosity()
             +
             (1 - upwindAlpha)*dn.concentration(1)/dn.viscosity());

        // diffusive flux of second component
        flux[transEqIdx] -=
            fluxVars.porousDiffCoeff() *
            (fluxVars.concentrationGrad(1) * fluxVars.face().normal);

        // dispersive flux of second component
        Vector normalDisp;
        fluxVars.dispersionTensor().mv(fluxVars.face().normal, normalDisp);
        flux[transEqIdx] -=
            (normalDisp * fluxVars.concentrationGrad(1));
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
};

}

#endif
