// $Id: 2p2clocalresidual.hh 3795 2010-06-25 16:08:04Z melanie $
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Bernd Flemisch                               *
 *   Copyright (C) 2009-2010 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUMUX_NEW_2P2C_BOX_JACOBIAN_BASE_HH
#define DUMUX_NEW_2P2C_BOX_JACOBIAN_BASE_HH

#include <dumux/boxmodels/common/boxmodel.hh>
#include <dumux/common/math.hh>

#include "2p2cproperties.hh"

#include "2p2csecondaryvars.hh"

#include "2p2cfluxvars.hh"

#include "2p2cnewtoncontroller.hh"

#include <iostream>
#include <vector>

//#define VELOCITY_OUTPUT 1 // uncomment this line if an output of the velocity is needed

namespace Dumux
{
/*!
 * \ingroup TwoPTwoCModel
 * \brief 2P-2C specific details needed to approximately calculate
 *        the local jacobian in the BOX scheme.
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
    typedef typename GridView::Grid::ctype CoordScalar;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementSolutionVector)) ElementSolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVarVector)) PrimaryVarVector;
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

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SecondaryVars)) SecondaryVars;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementSecondaryVars)) ElementSecondaryVars;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVars)) FluxVars;

    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> Tensor;

    static const Scalar mobilityUpwindAlpha =
            GET_PROP_VALUE(TypeTag, PTAG(MobilityUpwindAlpha));

public:
    /*!
     * \brief Evaluate the amount all conservation quantites
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     */
    void computeStorage(PrimaryVarVector &result, int scvIdx, bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const ElementSecondaryVars &elemDat = usePrevSol ? this->prevSecVars_()
                : this->curSecVars_();
        const SecondaryVars &vertDat = elemDat[scvIdx];

        // compute storage term of all components within all phases
        result = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                int eqIdx = (compIdx == lCompIdx) ? contiLEqIdx : contiGEqIdx;
                result[eqIdx] += vertDat.density(phaseIdx)
                        * vertDat.saturation(phaseIdx)
                        * vertDat.fluidState().massFrac(phaseIdx, compIdx);
            }
        }
        result *= vertDat.porosity();
    }

    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a subcontrol volume.
     */
    void computeFlux(PrimaryVarVector &flux, int faceIdx) const
    {
        FluxVars vars(this->problem_(), 
                      this->elem_(),
                      this->fvElemGeom_(), 
                      faceIdx, 
                      this->curSecVars_());

        flux = 0;
        asImp_()->computeAdvectiveFlux(flux, vars);
        asImp_()->computeDiffusiveFlux(flux, vars);

        // the face normal points into the outward direction, so we
        // have to multiply all fluxes with -1
        flux *= -1;
    }

    /*!
     * \brief Evaluates the advective mass flux of all components over
     *        a face of a subcontrol volume.
     */
    void computeAdvectiveFlux(PrimaryVarVector &flux, const FluxVars &vars) const
    {
        ////////
        // advective fluxes of all components in all phases
        ////////
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // data attached to upstream and the downstream vertices
            // of the current phase
            const SecondaryVars &up =
                    this->curSecVars_(vars.upstreamIdx(phaseIdx));
            const SecondaryVars &dn = this->curSecVars_()[vars.downstreamIdx(
                    phaseIdx)];

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
     */
    void computeDiffusiveFlux(PrimaryVarVector &flux, const FluxVars &vars) const
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
     */
    void computeSource(PrimaryVarVector &q, int localVertexIdx)
    {
        this->problem_().source(q,
                                this->elem_(),
                                this->fvElemGeom_(),
                                localVertexIdx);
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
