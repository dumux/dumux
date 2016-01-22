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
 *        using the two-phase n-component box model.
 */

#ifndef DUMUX_2PNC_LOCAL_RESIDUAL_BASE_HH
#define DUMUX_2PNC_LOCAL_RESIDUAL_BASE_HH

#include "2pncproperties.hh"
#include "2pncvolumevariables.hh"
#include <dumux/nonlinear/newtoncontroller.hh>

#include <iostream>
#include <vector>

namespace Dumux
{
/*!
 * \ingroup TwoPNCModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the two-phase n-component fully implicit box model.
 *
 * This class is used to fill the gaps in ImplicitLocalResidual for the two-phase n-component flow.
 */
template<class TypeTag>
class TwoPNCLocalResidual: public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
protected:
    typedef TwoPNCLocalResidual<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;
    typedef BoxLocalResidual<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, ElementSolutionVector) ElementSolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    typedef CompositionalFluidState<Scalar, FluidSystem> FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum
    {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),

        replaceCompEqIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx),

        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,

        wPhaseIdx = FluidSystem::wPhaseIdx,
        nPhaseIdx = FluidSystem::nPhaseIdx,

        wCompIdx = FluidSystem::wCompIdx,
        nCompIdx = FluidSystem::nCompIdx,

        conti0EqIdx = Indices::conti0EqIdx,

        wPhaseOnly = Indices::wPhaseOnly,
        nPhaseOnly = Indices::nPhaseOnly,
        bothPhases = Indices::bothPhases,

        plSg = TwoPNCFormulation::plSg,
        pgSl = TwoPNCFormulation::pgSl,
        formulation = GET_PROP_VALUE(TypeTag, Formulation)
    };

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;


    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::ctype CoordScalar;


public:
    /*!
     * \brief Constructor. Sets the upwind weight.
     */
    TwoPNCLocalResidual()
    {
        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        massUpwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);
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
        FVElementGeometry fvGeometry;
        fvGeometry.update(this->gridView_(), element);
        ElementBoundaryTypes bcTypes;
        bcTypes.update(this->problem_(), element, fvGeometry);
        ElementVolumeVariables volVars;
        volVars.update(this->problem_(), element, fvGeometry, false);

        this->residual_.resize(fvGeometry.numScv);
        this->residual_ = 0;

        this->elemPtr_ = &element;
        this->fvElemGeomPtr_ = &fvGeometry;
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
     *  \param storage the mass of the component within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    void computeStorage(PrimaryVariables &storage, int scvIdx, bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const ElementVolumeVariables &elemVolVars = usePrevSol ? this->prevVolVars_()
                : this->curVolVars_();
        const VolumeVariables &volVars = elemVolVars[scvIdx];

        // Compute storage term of all fluid components in the fluid phases
        storage = 0;
        for (unsigned int phaseIdx = 0; phaseIdx < numPhases /*+ numSPhases*/; ++phaseIdx)
        {
            //if(phaseIdx< numPhases)
            //{
                for (unsigned int compIdx = 0; compIdx < numComponents; ++compIdx) //H2O, Air, Salt
                {
                    int eqIdx = conti0EqIdx + compIdx;
                    if (replaceCompEqIdx != eqIdx)
                    {
                        storage[eqIdx] += volVars.molarDensity(phaseIdx)
                                        * volVars.saturation(phaseIdx)
                                        * volVars.moleFraction(phaseIdx, compIdx)
                                        * volVars.porosity();
                    }
                    else
                    {
                        storage[replaceCompEqIdx] += volVars.molarDensity(phaseIdx)
                                                * volVars.saturation(phaseIdx)
                                                * volVars.porosity();
                    }
                }
        }
         Valgrind::CheckDefined(storage);
    }
    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a sub-control volume.
     *
     * \param flux The flux over the sub-control-volume face for each component
     * \param fIdx The index of the sub-control-volume face
     * \param onBoundary Evaluate flux at inner sub-control-volume face or on a boundary face
     */
    void computeFlux(PrimaryVariables &flux, const int fIdx, bool onBoundary = false) const
    {
        FluxVariables fluxVars(this->problem_(),
                      this->element_(),
                      this->fvGeometry_(),
                      fIdx,
                      this->curVolVars_(),
                      onBoundary);

        flux = 0;
        asImp_()->computeAdvectiveFlux(flux, fluxVars);
        asImp_()->computeDiffusiveFlux(flux, fluxVars);
        Valgrind::CheckDefined(flux);
    }
    /*!
     * \brief Evaluates the advective mass flux of all components over
     *        a face of a subcontrol volume.
     *
     * \param flux The advective flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current SCV
     */
    void computeAdvectiveFlux(PrimaryVariables &flux, const FluxVariables &fluxVars) const
    {
        ////////
        // advective fluxes of all components in all phases
        ////////
        for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
         // data attached to upstream and the downstream vertices
         // of the current phase
         const VolumeVariables &up = this->curVolVars_(fluxVars.upstreamIdx(phaseIdx));
         const VolumeVariables &dn = this->curVolVars_(fluxVars.downstreamIdx(phaseIdx));

         for (unsigned int compIdx = 0; compIdx < numComponents; ++compIdx)
         {
         // add advective flux of current component in current
         // phase
            unsigned int eqIdx = conti0EqIdx + compIdx;

         if (eqIdx != replaceCompEqIdx)
         {
            // upstream vertex
            flux[eqIdx] += fluxVars.KmvpNormal(phaseIdx)
                        * (massUpwindWeight_
                            * up.mobility(phaseIdx)
                            * up.molarDensity(phaseIdx)
                            * up.moleFraction(phaseIdx, compIdx)
                        +
                            (1.0 - massUpwindWeight_)
                            * dn.mobility(phaseIdx)
                            * dn.molarDensity(phaseIdx)
                            * dn.moleFraction(phaseIdx, compIdx));

            Valgrind::CheckDefined(fluxVars.KmvpNormal(phaseIdx));
            Valgrind::CheckDefined(up.molarDensity(phaseIdx));
            Valgrind::CheckDefined(dn.molarDensity(phaseIdx));
         }
         else
         {
         flux[replaceCompEqIdx] += fluxVars.KmvpNormal(phaseIdx)
                                * (massUpwindWeight_
                                    * up.molarDensity(phaseIdx)
                                    * up.mobility(phaseIdx)
                                    +
                                    (1.0 - massUpwindWeight_)
                                    * dn.molarDensity(phaseIdx)
                                    * dn.mobility(phaseIdx));


         Valgrind::CheckDefined(fluxVars.KmvpNormal(phaseIdx));
         Valgrind::CheckDefined(up.molarDensity(phaseIdx));
         Valgrind::CheckDefined(dn.molarDensity(phaseIdx));
         }
         }
       }
     }

    /*!
     * \brief Evaluates the diffusive mass flux of all components over
     *        a face of a sub-control volume.
     *
     * \param flux The flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current sub-control-volume face
     */
    void computeDiffusiveFlux(PrimaryVariables &flux, const FluxVariables &fluxVars) const
    {
        //Loop calculates the diffusive flux for every component in a phase. The amount of moles of a component
        //(eg Air in liquid) in a phase
        //which is not the main component (eg. H2O in the liquid phase) moved from i to j equals the amount of moles moved
        //from the main component in a phase (eg. H2O in the liquid phase) from j to i. So two fluxes in each component loop
        // are calculated in the same phase.

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                //add diffusive fluxes only to the component balances
                if (replaceCompEqIdx != (conti0EqIdx + compIdx))
                {
                    Scalar diffCont = - fluxVars.porousDiffCoeff(phaseIdx ,compIdx)
                                        * fluxVars.molarDensity(phaseIdx)
                                        * (fluxVars.moleFractionGrad(phaseIdx, compIdx)
                                            * fluxVars.face().normal);
                    flux[conti0EqIdx + compIdx] += diffCont;
                    flux[conti0EqIdx + phaseIdx] -= diffCont;
                }
            }
    }

protected:

    void evalPhaseStorage_(int phaseIdx)
    {
        // evaluate the storage terms of a single phase
        for (int i=0; i < this->fvGeometry_().numScv; i++)
        {
            PrimaryVariables &result = this->residual_[i];
            const ElementVolumeVariables &elemVolVars = this->curVolVars_();
            const VolumeVariables &volVars = elemVolVars[i];

            // compute storage term of all fluid components within all phases
            result = 0;
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                result[conti0EqIdx + compIdx] += volVars.density(phaseIdx)
                                * volVars.saturation(phaseIdx)
                                * volVars.massFraction(phaseIdx, compIdx)
                                * volVars.porosity();
            }
            result *= this->fvGeometry_().subContVol[i].volume;
        }
    }

    Implementation *asImp_()
    { return static_cast<Implementation *> (this); }
    const Implementation *asImp_() const
    { return static_cast<const Implementation *> (this); }

public:
   Scalar massUpwindWeight_;
};

} // end namespace

#endif
