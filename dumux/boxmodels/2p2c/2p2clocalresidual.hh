// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
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

#ifndef DUMUX_NEW_2P2C_LOCAL_RESIDUAL_BASE_HH
#define DUMUX_NEW_2P2C_LOCAL_RESIDUAL_BASE_HH

#include <dumux/boxmodels/common/boxmodel.hh>
#include <dumux/common/math.hh>

#include "2p2cproperties.hh"
#include "2p2cvolumevariables.hh"
#include "2p2cfluxvariables.hh"
#include "2p2cnewtoncontroller.hh"

#include <iostream>
#include <vector>

//#define VELOCITY_OUTPUT 1 // uncomment this line if an output of the velocity is needed

namespace Dumux
{
/*!
 * \ingroup TwoPTwoCModel
 * \ingroup BoxLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the two-phase two-component box model.
 *
 * This class is used to fill the gaps in BoxLocalResidual for the 2P-2C flow.
 */
template<class TypeTag>
class TwoPTwoCLocalResidual: public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
 protected:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryVariables) BoundaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    enum
    {
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };

    typedef typename GET_PROP_TYPE(TypeTag, TwoPTwoCIndices) Indices;
    enum
    {
        contiLEqIdx = Indices::contiLEqIdx,
        contiGEqIdx = Indices::contiGEqIdx,
        lPhaseIdx = Indices::lPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,
        lCompIdx = Indices::lCompIdx,
        gCompIdx = Indices::gCompIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    static constexpr unsigned int replaceCompEqIdx =
        GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx);

 public:
    /*!
     * \brief Constructor. Sets the upwind weight.
     */
    TwoPTwoCLocalResidual()
    {
        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        massUpwindWeight_ = GET_PARAM(TypeTag, Scalar, MassUpwindWeight);
    };

    /*!
     * \brief Evaluate the storage term of the current solution in a
     *        single phase.
     *
     * \param element The element
     * \param phaseIdx The index of the fluid phase
     */
    void evalPhaseStorage(const Element &element, int phaseIdx)
    {
        FVElementGeometry fvGeom;
        fvGeom.update(this->gridView_(), element);
        ElementBoundaryTypes bcTypes;
        bcTypes.update(this->problem_(), element, fvGeom);
        ElementVolumeVariables volVars;
        volVars.update(this->problem_(), element, fvGeom, false);

        this->storageTerm_.resize(fvGeom.numVertices);
        this->storageTerm_ = 0;

        this->elemPtr_ = &element;
        this->fvElemGeomPtr_ = &fvGeom;
        this->bcTypesPtr_ = &bcTypes;
        this->prevVolVarsPtr_ = 0;
        this->curVolVarsPtr_ = &volVars;
        evalPhaseStorage_(phaseIdx);
    }

    /*!
     * \brief Evaluate the outflow conditions at a boundary face.
     *        This is still a beta version and only available
     *        for the 2p2c and the 2p2cni model yet.
     *
     * \param isIt The intersection iterator
     * \param scvIdx The index of the sub-control volume containing the outflow face
     * \param boundaryFaceIdx  The index of the boundary face
     */
    void evalOutflowSegment(const IntersectionIterator &isIt,
                            int scvIdx,
                            int boundaryFaceIdx)
    {
        // temporary vector to store the outflow boundary fluxes
        PrimaryVariables flux(0.0);
        const BoundaryTypes &bcTypes = this->bcTypes_(scvIdx);

        // deal with outflow boundaries
        if (bcTypes.hasOutflow())
        {
            const BoundaryVariables boundaryVars(this->problem_(),
                                                 this->elem_(),
                                                 this->fvElemGeom_(),
                                                 boundaryFaceIdx,
                                                 this->curVolVars_());

            PrimaryVariables values(0.0);
            asImp_()->computeOutflowValues(values, boundaryVars, scvIdx, boundaryFaceIdx);
            Valgrind::CheckDefined(values);

            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            {
                if (!bcTypes.isOutflow(eqIdx) )
                    continue;
                this->residual_[scvIdx][eqIdx] += values[eqIdx];
            }
        }
    }

    /*!
     * \brief Evaluate the amount all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub-control volume divided by the volume)
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
        const ElementVolumeVariables &elemVolVars = usePrevSol ? this->prevVolVars_()
            : this->curVolVars_();
        const VolumeVariables &volVars = elemVolVars[scvIdx];

        // compute storage term of all components within all phases
        result = 0;

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (int compIdx = contiCompIdx1_(); compIdx <= contiCompIdx2_(); ++compIdx)
            {
                int eqIdx = (compIdx == lCompIdx) ? contiLEqIdx : contiGEqIdx;
                result[eqIdx] += volVars.density(phaseIdx)
                    * volVars.saturation(phaseIdx)
                    * volVars.fluidState().massFraction(phaseIdx, compIdx);
            }
            // this is only processed, if one component mass balance equation
            // is replaced by the total mass balance equation
            if (replaceCompEqIdx < numComponents)
                result[replaceCompEqIdx] +=
                    volVars.density(phaseIdx)
                    * volVars.saturation(phaseIdx);
        }
        result *= volVars.porosity();
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
        Valgrind::CheckDefined(flux);
        asImp_()->computeDiffusiveFlux(flux, vars);
        Valgrind::CheckDefined(flux);
    }

    /*!
     * \brief Evaluates the advective mass flux of all components over
     *        a face of a sub-control volume.
     *
     * \param flux The advective flux over the sub-control-volume face for each component
     * \param vars The flux variables at the current SCV face
     */
    void computeAdvectiveFlux(PrimaryVariables &flux, const FluxVariables &vars) const
    {
        ////////
        // advective fluxes of all components in all phases
        ////////
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // data attached to upstream and the downstream vertices
            // of the current phase
            const VolumeVariables &up =
                this->curVolVars_(vars.upstreamIdx(phaseIdx));
            const VolumeVariables &dn =
                this->curVolVars_(vars.downstreamIdx(phaseIdx));

            for (int compIdx = contiCompIdx1_(); compIdx <= contiCompIdx2_(); ++compIdx)
            {
                int eqIdx = (compIdx == lCompIdx) ? contiLEqIdx : contiGEqIdx;
                // add advective flux of current component in current
                // phase
                if (massUpwindWeight_ > 0.0)
                    // upstream vertex
                    flux[eqIdx] +=
                        vars.KmvpNormal(phaseIdx)
                        * massUpwindWeight_
                        * up.density(phaseIdx)
                        * up.mobility(phaseIdx)
                        * up.fluidState().massFraction(phaseIdx, compIdx);
                if (massUpwindWeight_ < 1.0)
                    // downstream vertex
                    flux[eqIdx] +=
                        vars.KmvpNormal(phaseIdx)
                        * (1 - massUpwindWeight_)
                        * dn.density(phaseIdx)
                        * dn.mobility(phaseIdx)
                        * dn.fluidState().massFraction(phaseIdx, compIdx);

                Valgrind::CheckDefined(vars.KmvpNormal(phaseIdx));
                Valgrind::CheckDefined(up.density(phaseIdx));
                Valgrind::CheckDefined(up.mobility(phaseIdx));
                Valgrind::CheckDefined(up.fluidState().massFraction(phaseIdx, compIdx));
                Valgrind::CheckDefined(dn.density(phaseIdx));
                Valgrind::CheckDefined(dn.mobility(phaseIdx));
                Valgrind::CheckDefined(dn.fluidState().massFraction(phaseIdx, compIdx));
            }
            // flux of the total mass balance;
            // this is only processed, if one component mass balance equation
            // is replaced by a total mass balance equation
            if (replaceCompEqIdx < numComponents)
            {
                // upstream vertex
                if (massUpwindWeight_ > 0.0)
                    flux[replaceCompEqIdx] +=
                        vars.KmvpNormal(phaseIdx)
                        * massUpwindWeight_
                        * up.density(phaseIdx)
                        * up.mobility(phaseIdx);
                // downstream vertex
                if (massUpwindWeight_ < 1.0)
                    flux[replaceCompEqIdx] +=
                        vars.KmvpNormal(phaseIdx)
                        * (1 - massUpwindWeight_)
                        * dn.density(phaseIdx)
                        * dn.mobility(phaseIdx);
                Valgrind::CheckDefined(vars.KmvpNormal(phaseIdx));
                Valgrind::CheckDefined(up.density(phaseIdx));
                Valgrind::CheckDefined(up.mobility(phaseIdx));
                Valgrind::CheckDefined(dn.density(phaseIdx));
                Valgrind::CheckDefined(dn.mobility(phaseIdx));

            }

        }
    }

    /*!
     * \brief Adds the diffusive mass flux of all components over
     *        a face of a sub-control volume.
     *
     * \param flux The diffusive flux over the sub-control-volume face for each component
     * \param vars The flux variables at the current sub control volume face
     */
    void computeDiffusiveFlux(PrimaryVariables &flux, const FluxVariables &vars) const
    {
        // add diffusive flux of gas component in liquid phase
        Scalar tmp = vars.molarConcGrad(lPhaseIdx)*vars.face().normal;
        tmp *= -1;
        tmp *=
            vars.porousDiffCoeff(lPhaseIdx) *
            vars.molarDensityAtIP(lPhaseIdx);
        // add the diffusive fluxes only to the component mass balance
        if (replaceCompEqIdx != contiGEqIdx)
            flux[contiGEqIdx] += tmp * FluidSystem::molarMass(gCompIdx);
        if (replaceCompEqIdx != contiLEqIdx)
            flux[contiLEqIdx] -= tmp * FluidSystem::molarMass(lCompIdx);

        // add diffusive flux of liquid component in gas phase
        tmp = vars.molarConcGrad(gPhaseIdx)*vars.face().normal;
        tmp *= -1;
        tmp *=
            vars.porousDiffCoeff(gPhaseIdx) *
            vars.molarDensityAtIP(gPhaseIdx);
        // add the diffusive fluxes only to the component mass balance
        if (replaceCompEqIdx != contiLEqIdx)
            flux[contiLEqIdx] += tmp * FluidSystem::molarMass(lCompIdx);
        if (replaceCompEqIdx != contiGEqIdx)
            flux[contiGEqIdx] -= tmp * FluidSystem::molarMass(gCompIdx);
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param q The source/sink in the sub-control volume for each component
     * \param localVertexIdx The index of the sub-control volume
     */
    void computeSource(PrimaryVariables &q, int localVertexIdx)
    {
        this->problem_().boxSDSource(q,
                                     this->elem_(),
                                     this->fvElemGeom_(),
                                     localVertexIdx,
                                     this->curVolVars_());
    }

    /*!
     * \brief Compute the fluxes at outflow boundaries. This does essentially the same
     *        as computeFluxes, but the fluxes are evaluated at the integration point
     *        of the boundary face. However, some variables are evaluated
     *        at the vertex (usually the ones which are upwinded).
     *
     * \param flux A temporary vector, where the outflow boundary fluxes are stored
     * \param boundaryVars The boundary variables object
     * \param scvIdx The index of the SCV containing the outflow boundary face
     * \param boundaryFaceIdx The index of the boundary face
     */
    void computeOutflowValues(PrimaryVariables &flux,
                              const BoundaryVariables &boundaryVars,
                              const int scvIdx,
                              const int boundaryFaceIdx)
    {
        const BoundaryTypes &bcTypes = this->bcTypes_(scvIdx);
        const VolumeVariables& vertVars = this->curVolVars_()[scvIdx];

        // component transport
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (int compIdx = contiCompIdx1_(); compIdx <= contiCompIdx2_(); ++compIdx)
            {
                int eqIdx = (compIdx == lCompIdx) ? contiLEqIdx : contiGEqIdx;
                if (bcTypes.isOutflow(eqIdx))
                {
                    // advective flux
                    flux[eqIdx] += boundaryVars.KmvpNormal(phaseIdx)
                        * vertVars.density(phaseIdx)
                        * vertVars.mobility(phaseIdx)
                        * vertVars.fluidState().massFraction(phaseIdx, compIdx);
                }
            }
        }

        // add diffusive flux of gas component in liquid phase
        Scalar tmp = boundaryVars.molarConcGrad(lPhaseIdx)*boundaryVars.boundaryFace().normal;
        tmp *= -1;
        tmp *= boundaryVars.porousDiffCoeff(lPhaseIdx) *
            boundaryVars.molarDensityAtIP(lPhaseIdx);
        // add the diffusive fluxes only to the component mass balance
        if (replaceCompEqIdx != contiGEqIdx && bcTypes.isOutflow(contiGEqIdx))
            flux[contiGEqIdx] += tmp * FluidSystem::molarMass(gCompIdx);
        if (replaceCompEqIdx != contiLEqIdx && bcTypes.isOutflow(contiLEqIdx))
            flux[contiLEqIdx] -= tmp * FluidSystem::molarMass(lCompIdx);

        // add diffusive flux of liquid component in gas phase
        tmp = boundaryVars.molarConcGrad(gPhaseIdx)*boundaryVars.boundaryFace().normal;
        tmp *= -1;
        tmp *= boundaryVars.porousDiffCoeff(gPhaseIdx) *
            boundaryVars.molarDensityAtIP(gPhaseIdx);
        // add the diffusive fluxes only to the component mass balance
        if (replaceCompEqIdx != contiLEqIdx && bcTypes.isOutflow(contiLEqIdx))
            flux[contiLEqIdx] += tmp * FluidSystem::molarMass(lCompIdx);
        if (replaceCompEqIdx != contiGEqIdx && bcTypes.isOutflow(contiGEqIdx))
            flux[contiGEqIdx] -= tmp * FluidSystem::molarMass(gCompIdx);

        // flux of the total mass balance;
        // this is only processed, if one component mass-balance equation
        // is replaced by a total mass-balance equation
        if (replaceCompEqIdx < numComponents)
        {
            if (bcTypes.isOutflow(replaceCompEqIdx))
                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                    flux[replaceCompEqIdx] += boundaryVars.KmvpNormal(phaseIdx)
                        *vertVars.density(phaseIdx)
                        *vertVars.mobility(phaseIdx);
        }
    }


 protected:
    void evalPhaseStorage_(int phaseIdx)
    {
        // evaluate the storage terms of a single phase
        for (int i=0; i < this->fvElemGeom_().numVertices; i++) {
            PrimaryVariables &result = this->storageTerm_[i];
            const ElementVolumeVariables &elemVolVars = this->curVolVars_();
            const VolumeVariables &volVars = elemVolVars[i];

            // compute storage term of all components within all phases
            result = 0;
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                int eqIdx = (compIdx == lCompIdx) ? contiLEqIdx : contiGEqIdx;
                result[eqIdx] += volVars.density(phaseIdx)
                    * volVars.saturation(phaseIdx)
                    * volVars.fluidState().massFraction(phaseIdx, compIdx);
            }

            result *= volVars.porosity();
            result *= this->fvElemGeom_().subContVol[i].volume;
        }
    }

    /*!
     * \brief Return the equation index of the first mass-balance equation
     *        of the component (used for loops); if one component mass balance
     *        is replaced by the total mass balance, this is the index
     *        of the remaining component mass-balance equation.
     */
    constexpr unsigned int contiCompIdx1_() const {
        switch (replaceCompEqIdx)
        {
        case contiLEqIdx: return contiGEqIdx;
        case contiGEqIdx: return contiLEqIdx;
        default:          return 0;
        }
    }

    /*!
     * \brief Return the equation index of the second mass balance
     *        of the component (used for loops);
     *        if one component mass balance is replaced by the total mass balance
     *        (replaceCompEqIdx < 2), this index is the same as contiCompIdx1().
     */
    constexpr unsigned int contiCompIdx2_() const {
        switch (replaceCompEqIdx)
        {
        case contiLEqIdx: return contiGEqIdx;
        case contiGEqIdx: return contiLEqIdx;
        default:          return numComponents-1;
        }
    }

    Implementation *asImp_()
    { return static_cast<Implementation *> (this); }
    const Implementation *asImp_() const
    { return static_cast<const Implementation *> (this); }

 private:
    Scalar massUpwindWeight_;
};

} // end namespace

#endif
