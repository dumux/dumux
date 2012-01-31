// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Holger Class                                 *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Bernd Flemisch                               *
 *   Copyright (C) 2009-2010 by Andreas Lauser                               *
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
 *        using the two-phase two-component box model.
 */
#ifndef DUMUX_3P3C_LOCAL_RESIDUAL_HH
#define DUMUX_3P3C_LOCAL_RESIDUAL_HH

#include <dumux/boxmodels/common/boxmodel.hh>
#include <dumux/common/math.hh>

#include "3p3cproperties.hh"
#include "3p3cvolumevariables.hh"
#include "3p3cfluxvariables.hh"
#include "3p3cnewtoncontroller.hh"

#include <iostream>
#include <vector>

//#define VELOCITY_OUTPUT 1 // uncomment this line if an output of the velocity is needed

namespace Dumux
{
/*!
 * \ingroup ThreePThreeCModel
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the two-phase two-component box model.
 *
 * This class is used to fill the gaps in BoxLocalResidual for the 3P3C flow.
 */
template<class TypeTag>
class ThreePThreeCLocalResidual: public BoxLocalResidual<TypeTag>
{
protected:
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ThreePThreeCIndices) Indices;

    enum {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),

        conti0EqIdx = Indices::conti0EqIdx,
        contiWEqIdx = Indices::contiWEqIdx,
        contiCEqIdx = Indices::contiCEqIdx,
        contiAEqIdx = Indices::contiAEqIdx,

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,

        wCompIdx = Indices::wCompIdx,
        cCompIdx = Indices::cCompIdx,
        aCompIdx = Indices::aCompIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;

public:
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
        // using the implicit euler method.
        const ElementVolumeVariables &elemVolVars = 
            usePrevSol
            ? this->prevVolVars_()
            : this->curVolVars_();
        const VolumeVariables &volVars = elemVolVars[scvIdx];

        // compute storage term of all components within all phases
        result = 0;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                result[conti0EqIdx + compIdx] += 
                    volVars.porosity() 
                    * volVars.saturation(phaseIdx)
                    * volVars.molarDensity(phaseIdx)
                    * volVars.fluidState().moleFraction(phaseIdx, compIdx);
            }
        }
    }

    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a sub-control volume.
     *
     * \param flux The flux over the SCV (sub-control-volume) face for each component
     * \param faceIdx The index of the SCV face
     */
    void computeFlux(PrimaryVariables &flux, int faceIdx) const
    {
        FluxVariables vars(this->problem_(),
                           this->elem_(),
                           this->fvElemGeom_(),
                           faceIdx,
                           this->curVolVars_());

        flux = 0;
        asImp_()->computeAdvectiveFlux(flux, vars);
        asImp_()->computeDiffusiveFlux(flux, vars);
    }

    /*!
     * \brief Evaluates the advective mass flux of all components over
     *        a face of a subcontrol volume.
     *
     * \param flux The advective flux over the sub-control-volume face for each component
     * \param vars The flux variables at the current SCV
     */

    void computeAdvectiveFlux(PrimaryVariables &flux, const FluxVariables &vars) const
    {
        Scalar massUpwindWeight = GET_PARAM(TypeTag, Scalar, MassUpwindWeight);

        ////////
        // advective fluxes of all components in all phases
        ////////
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // data attached to upstream and the downstream vertices
            // of the current phase
            const VolumeVariables &up = this->curVolVars_(vars.upstreamIdx(phaseIdx));
            const VolumeVariables &dn = this->curVolVars_(vars.downstreamIdx(phaseIdx));

            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                // add advective flux of current component in current
                // phase
                // if alpha > 0 und alpha < 1 then both upstream and downstream
                // nodes need their contribution
                // if alpha == 1 (which is mostly the case) then, the downstream
                // node is not evaluated
                int eqIdx = conti0EqIdx + compIdx;
                flux[eqIdx] += vars.KmvpNormal(phaseIdx)
                    * (massUpwindWeight
                       * up.mobility(phaseIdx)
                       * up.fluidState().molarDensity(phaseIdx)
                       * up.fluidState().moleFraction(phaseIdx, compIdx)
                       +
                       (1.0 - massUpwindWeight)
                       * dn.mobility(phaseIdx)
                       * dn.fluidState().molarDensity(phaseIdx)
                       * dn.fluidState().moleFraction(phaseIdx, compIdx));
            }
        }
    }

    /*!
     * \brief Adds the diffusive mass flux of all components over
     *        a face of a subcontrol volume.
     *
     * \param flux The diffusive flux over the sub-control-volume face for each component
     * \param vars The flux variables at the current SCV
     */

    void computeDiffusiveFlux(PrimaryVariables &flux, const FluxVariables &vars) const
    {
#warning "Disabled for now"
        return;
        // TODO: reference!?  Dune::FieldMatrix<Scalar, numPhases, numComponents> averagedPorousDiffCoeffMatrix = vars.porousDiffCoeff();
        // add diffusive flux of gas component in liquid phase
        Scalar tmp = - vars.porousDiffCoeff()[wPhaseIdx][aCompIdx] * vars.molarDensityAtIP(wPhaseIdx);
        tmp *= (vars.molarAConcGrad(wPhaseIdx) * vars.face().normal);
        Scalar jAW = tmp;

        tmp = - vars.porousDiffCoeff()[wPhaseIdx][cCompIdx] * vars.molarDensityAtIP(wPhaseIdx);
        tmp *= (vars.molarCConcGrad(wPhaseIdx) * vars.face().normal);
        Scalar jCW = tmp;

        Scalar jWW = -(jAW+jCW);

        tmp = - vars.porousDiffCoeff()[gPhaseIdx][wCompIdx] * vars.molarDensityAtIP(gPhaseIdx);
        tmp *= (vars.molarWConcGrad(gPhaseIdx) * vars.face().normal);
        Scalar jWG = tmp;

        tmp = - vars.porousDiffCoeff()[gPhaseIdx][cCompIdx] * vars.molarDensityAtIP(gPhaseIdx);
        tmp *= (vars.molarCConcGrad(gPhaseIdx) * vars.face().normal);
        Scalar jCG = tmp;

        Scalar jAG = -(jWG+jCG);

        tmp = - vars.porousDiffCoeff()[nPhaseIdx][wCompIdx] * vars.molarDensityAtIP(nPhaseIdx);
        tmp *= (vars.molarWConcGrad(nPhaseIdx) * vars.face().normal);
        Scalar jWN = tmp;

        tmp = - vars.porousDiffCoeff()[nPhaseIdx][aCompIdx] * vars.molarDensityAtIP(nPhaseIdx);
        tmp *= (vars.molarAConcGrad(nPhaseIdx) * vars.face().normal);
        Scalar jAN = tmp;

        Scalar jCN = -(jAN+jWN);

        /*
          jWW = 0;
          jWG = 0;
          jWN = 0;
          jCW = 0;
          jCG = 0;
          jCN = 0;
          jAW = 0;
          jAG = 0;
          jAN = 0;
        */

        flux[contiWEqIdx] += jWW+jWG+jWN;
        flux[contiCEqIdx] += jCW+jCG+jCN;
        flux[contiAEqIdx] += jAW+jAG+jAN;
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param q The source/sink in the SCV for each component
     * \param localVertexIdx The index of the SCV
     */
    void computeSource(PrimaryVariables &q, int localVertexIdx)
    {
        this->problem_().boxSDSource(q,
                                     this->elem_(),
                                     this->fvElemGeom_(),
                                     localVertexIdx,
                                     this->curVolVars_());
    }

protected:
    Implementation *asImp_()
    {
        return static_cast<Implementation *> (this);
    }

    const Implementation *asImp_() const
    {
        return static_cast<const Implementation *> (this);
    }
};

} // end namepace

#endif
