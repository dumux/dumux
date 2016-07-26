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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Element-wise calculation of the residual for the two-phase fully implicit model.
 */
#ifndef DUMUX_TWOPMINC_LOCAL_RESIDUAL_BASE_HH
#define DUMUX_TWOPMINC_LOCAL_RESIDUAL_BASE_HH

#include <dumux/porousmediumflow/2p/implicit/localresidual.hh>
#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup TwoPBoxModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the two-phase fully implicit model.
 *
 * This class is also used for the non-isothermal model, which means
 * that it uses static polymorphism.
 */
template<class TypeTag>
class TwoPMincLocalResidual : public TwoPLocalResidual<TypeTag>
{
protected:
    typedef TwoPLocalResidual<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum
    {
        contiWEqIdx = Indices::contiWEqIdx,
        contiNEqIdx = Indices::contiNEqIdx,
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        //numEq = GET_PROP_VALUE(TypeTag, NumEq),
        numContinua = GET_PROP_VALUE(TypeTag, NumContinua),
        fractureIdx = 0, /*!< index of the fracture (continua are indexed nC = 0,1,..numContinua-1)*/
        // Options for choosing the nested volume elements
        constantVolFraction_ = 0,
        equidistantNestedVolElements_ = 1,
        DFM_volFraction_ = 2
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    /*!
     * \brief Constructor. Sets the upwind weight.
     */
    TwoPMincLocalResidual():
        interactingContinuaType_(GET_PARAM_FROM_GROUP(TypeTag, int, Problem, InteractingContinuaType))
    {
        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        massUpwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);
    }

    /*!
     * \brief Evaluate the amount all conservation quantities
     *        (e.g. phase mass) within a finite sub-control volume.
     *
     *  \param storage The phase mass within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    void computeStorage(PrimaryVariables &storage, int scvIdx, bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit Euler method.
        const ElementVolumeVariables &elemVolVars = usePrevSol ? this->prevVolVars_() : this->curVolVars_();
        const VolumeVariables &volVars = elemVolVars[scvIdx];
        for (int cIdx = 0; cIdx < numContinua; ++cIdx)
        {
            // wetting phase mass
            storage[Indices::contiWEqIdxc(cIdx)] = volVars.density(wPhaseIdx,cIdx)
                    * volVars.porosity(cIdx)
                    * volVars.volumeFraction(cIdx)
                    * volVars.saturation(wPhaseIdx,cIdx);

            // non-wetting phase mass
            storage[Indices::contiNEqIdxc(cIdx)] = volVars.density(nPhaseIdx,cIdx)
                    * volVars.porosity(cIdx)
                    * volVars.volumeFraction(cIdx)
                    * volVars.saturation(nPhaseIdx,cIdx);
        }
    }

    /*!
     * \brief Evaluates the advective mass flux of all components over
     *        a face of a sub-control volume.
     *
     * This method is called by compute flux and is mainly there for
     * derived models to ease adding equations selectively.
     * The upwind direction is determined by the (higher permeable) fracture
     * continuum (nC=0).
     *
     * \param flux The advective flux over the sub-control-volume face for each phase
     * \param fluxVars The flux variables at the current SCV
     */
    void computeAdvectiveFlux(PrimaryVariables &flux, const FluxVariables &fluxVars) const
    {
        // loop over all phases
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // data attached to upstream and the downstream vertices
            // of the current phase
            const VolumeVariables &up = this->curVolVars_(fluxVars.upstreamIdx(phaseIdx));
            const VolumeVariables &dn = this->curVolVars_(fluxVars.downstreamIdx(phaseIdx));

            // add advective flux of current phase
            int eqIdx = (phaseIdx == wPhaseIdx) ? contiWEqIdx : contiNEqIdx;
            flux[eqIdx] +=
                fluxVars.volumeFlux(phaseIdx)
                *
                ((    massUpwindWeight_)*up.density(phaseIdx, fractureIdx)
                 +
                 (1 - massUpwindWeight_)*dn.density(phaseIdx, fractureIdx));
        }
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param q The source/sink in the SCV for each phase
     * \param scvIdx The index of the SCV
     *
     */
    void computeSource(PrimaryVariables &q, int scvIdx)
    {
        ParentType::computeSource(q, scvIdx);
        q += transferFlux(scvIdx);
    }

    PrimaryVariables transferFlux(int scvIdx)
    {
        PrimaryVariables mincInterContinuaFlux;
        // retrieve the source term intrinsic to the problem
        const ElementVolumeVariables &elemDat = this->curVolVars_();
        const VolumeVariables &volVars = elemDat[scvIdx];
        pGrad_ = 0.0;
        for (int nC = 0; nC<numContinua-1; nC++)
        {
            /* Calculating the pressure gradients in the coarse block between the
             * successive nested elements
             */
            pGrad_[wPhaseIdx][nC] =
                    volVars.pressure(wPhaseIdx,nC)
                - volVars.pressure(wPhaseIdx,nC+1);

            pGrad_[nPhaseIdx][nC] =
                    volVars.pressure(nPhaseIdx,nC)
                - volVars.pressure(nPhaseIdx,nC+1);

            /*
             * Two approaches:
             * 1) the effective permeabilities are defined in the test/soil file
             * 2) the transmissibilities between nested volume elements are
             *    calculated in 2pFract>1phase - NoFlowProblem and written to
             *    transmissibility.csv
             *
             * if the transmissivites file given, than use those values
             *
             * Transmissibility is in this case T = K*A_ij / d_ij
             * K - permeability, A_ij - inteface between continuum i and
             * continuum j, and d_ij distance between continuum i and j.
             * Coming from grad P = (p_i - p_j)/d_ij
             */
            //wetting phase
            pGrad_[wPhaseIdx][nC] *= volVars.intrinsicPermeability(nC);
            pGrad_[wPhaseIdx][nC] /= this->problem_().model().getDistNestedContinua()[nC+1];
            pGrad_[wPhaseIdx][nC] *= this->problem_().model().getInterfaceArea()[nC];

            //non-wetting phase
            pGrad_[nPhaseIdx][nC] *= volVars.intrinsicPermeability(nC);
            pGrad_[nPhaseIdx][nC] /= this->problem_().model().getDistNestedContinua()[nC+1];
            pGrad_[nPhaseIdx][nC] *= this->problem_().model().getInterfaceArea()[nC];

            Valgrind::CheckDefined(pGrad_);
            Scalar mobilityW = 0;
            Scalar mobilityN = 0;

            if (pGrad_[wPhaseIdx][nC] > 0) // full upwinding
            {
                mobilityW = volVars.mobility(wPhaseIdx,nC);
            }
            else
            {
                mobilityW = volVars.mobility(wPhaseIdx,nC+1);
            }

            pGrad_[wPhaseIdx][nC] *= mobilityW;
            pGrad_[wPhaseIdx][nC] *= volVars.density(wPhaseIdx,nC);

            if (pGrad_[nPhaseIdx][nC] > 0) // full upwinding
            {
                mobilityN = volVars.mobility(nPhaseIdx,nC);
            }
            else
            {
                mobilityN = volVars.mobility(nPhaseIdx,nC+1);
            }

            pGrad_[nPhaseIdx][nC] *= mobilityN;
            pGrad_[nPhaseIdx][nC] *= volVars.density(nPhaseIdx,nC);
            Valgrind::CheckDefined(pGrad_);

            //wetting phase
            mincInterContinuaFlux[Indices::contiWEqIdxc(nC)] -= pGrad_[wPhaseIdx][nC];
            mincInterContinuaFlux[Indices::contiWEqIdxc(nC+1)] += pGrad_[wPhaseIdx][nC];

            //non-wetting phase
            mincInterContinuaFlux[Indices::contiNEqIdxc(nC)] -= pGrad_[nPhaseIdx][nC];
            mincInterContinuaFlux[Indices::contiNEqIdxc(nC+1)] += pGrad_[nPhaseIdx][nC];

        }

        return (mincInterContinuaFlux);
    }

private:
    Dune::FieldMatrix<Scalar, numPhases, numContinua> pGrad_;
    Scalar massUpwindWeight_;
    int interactingContinuaType_;

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

}

#endif
