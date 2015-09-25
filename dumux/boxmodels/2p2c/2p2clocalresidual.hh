// $Id: 2p2clocalresidual.hh 3795 2010-06-25 16:08:04Z melanie $
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
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the two-phase two-component box model.
 *
 * This class is used to fill the gaps in BoxLocalResidual for the 2P-2C flow.
 */
template<class TypeTag>
class TwoPTwoCLocalResidual: public BoxLocalResidual<TypeTag>
{
protected:
    typedef TwoPTwoCLocalResidual<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalResidual)) Implementation;
    typedef BoxLocalResidual<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementSolutionVector)) ElementSolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;

    typedef TwoPTwoCFluidState<TypeTag> FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices)) Indices;

    enum
    {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),
        numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)),

        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,

        contiLEqIdx = Indices::contiLEqIdx,
        contiGEqIdx = Indices::contiGEqIdx,

        lPhaseIdx = Indices::lPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,

        lCompIdx = Indices::lCompIdx,
        gCompIdx = Indices::gCompIdx,

        lPhaseOnly = Indices::lPhaseOnly,
        gPhaseOnly = Indices::gPhaseOnly,
        bothPhases = Indices::bothPhases,

        plSg = TwoPTwoCFormulation::plSg,
        pgSl = TwoPTwoCFormulation::pgSl,
        formulation = GET_PROP_VALUE(TypeTag, PTAG(Formulation))
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementBoundaryTypes)) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVariables)) FluxVariables;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::ctype CoordScalar;

    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> Tensor;

    static const Scalar mobilityUpwindAlpha =
            GET_PROP_VALUE(TypeTag, PTAG(MobilityUpwindAlpha));

public:
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

        this->residual_.resize(fvGeom.numVertices);
        this->residual_ = 0;

        this->elemPtr_ = &element;
        this->fvElemGeomPtr_ = &fvGeom;
        this->bcTypesPtr_ = &bcTypes;
        this->prevVolVarsPtr_ = 0;
        this->curVolVarsPtr_ = &volVars;
        evalPhaseStorage_(phaseIdx);
    }

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
        const ElementVolumeVariables &elemVolVars = usePrevSol ? this->prevVolVars_()
                : this->curVolVars_();
        const VolumeVariables &volVars = elemVolVars[scvIdx];

        // compute storage term of all components within all phases
        result = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                int eqIdx = (compIdx == lCompIdx) ? contiLEqIdx : contiGEqIdx;
                result[eqIdx] += volVars.density(phaseIdx)
                        * volVars.saturation(phaseIdx)
                        * volVars.fluidState().massFrac(phaseIdx, compIdx);
            }
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

            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                int eqIdx = (compIdx == lCompIdx) ? contiLEqIdx : contiGEqIdx;
                // add advective flux of current component in current
                // phase
                if (mobilityUpwindAlpha > 0.0)
                    // upstream vertex
                    flux[eqIdx] += vars.KmvpNormal(phaseIdx)
                            * mobilityUpwindAlpha * (up.density(phaseIdx)
                            * up.mobility(phaseIdx) * up.fluidState().massFrac(
                            phaseIdx, compIdx));
                if (mobilityUpwindAlpha < 1.0)
                    // downstream vertex
                    flux[eqIdx] += vars.KmvpNormal(phaseIdx) * (1
                            - mobilityUpwindAlpha) * (dn.density(phaseIdx)
                            * dn.mobility(phaseIdx) * dn.fluidState().massFrac(
                            phaseIdx, compIdx));
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
        // add diffusive flux of gas component in liquid phase
        Scalar tmp =
            - vars.porousDiffCoeff(lPhaseIdx) *
            vars.molarDensityAtIP(lPhaseIdx) *
            (vars.molarConcGrad(lPhaseIdx) * vars.face().normal);
        flux[contiGEqIdx] += tmp * FluidSystem::molarMass(gCompIdx);
        flux[contiLEqIdx] -= tmp * FluidSystem::molarMass(lCompIdx);

        // add diffusive flux of liquid component in gas phase
        tmp =
            - vars.porousDiffCoeff(gPhaseIdx) *
            vars.molarDensityAtIP(gPhaseIdx) *
            (vars.molarConcGrad(gPhaseIdx) * vars.face().normal);
        flux[contiLEqIdx] += tmp * FluidSystem::molarMass(lCompIdx);
        flux[contiGEqIdx] -= tmp * FluidSystem::molarMass(gCompIdx);
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param q The source/sink in the SCV for each component
     * \param localVertexIdx The index of the SCV
     */
    void computeSource(PrimaryVariables &q, int localVertexIdx)
    {
        this->problem_().source(q,
                                this->elem_(),
                                this->fvElemGeom_(),
                                localVertexIdx);
    }

protected:
    void evalPhaseStorage_(int phaseIdx)
    {
        // evaluate the storage terms of a single phase
        for (int i=0; i < this->fvElemGeom_().numVertices; i++) {
            PrimaryVariables &result = this->residual_[i];
            const ElementVolumeVariables &elemVolVars = this->curVolVars_();
            const VolumeVariables &volVars = elemVolVars[i];

            // compute storage term of all components within all phases
            result = 0;
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                int eqIdx = (compIdx == lCompIdx) ? contiLEqIdx : contiGEqIdx;
                result[eqIdx] += volVars.density(phaseIdx)
                    * volVars.saturation(phaseIdx)
                    * volVars.fluidState().massFrac(phaseIdx, compIdx);
            }

            result *= volVars.porosity();
            result *= this->fvElemGeom_().subContVol[i].volume;
        }
    }


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
