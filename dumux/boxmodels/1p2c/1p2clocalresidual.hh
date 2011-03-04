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
#define VELOCITY_OUTPUT 1 //1 turns velocity output on, 0 turns it off

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
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalResidual)) Implementation;
    typedef BoxLocalResidual<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVariables)) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePTwoCIndices)) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

    enum
    {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),

        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,
        x1Idx = Indices::x1Idx,

        phaseIdx = Indices::phaseIdx,
        comp0Idx = Indices::comp0Idx,
        comp1Idx = Indices::comp1Idx,

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

        // storage term of the transport equation - molefractions
//        result[transEqIdx] =
//            volVars.concentration(comp1Idx) *
//            volVars.porosity();

        //storage term of the transport equation - massfractions
        result[transEqIdx] =
                   volVars.density() * volVars.massFrac(comp1Idx) * volVars.porosity();

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

        asImp_()->computeAdvectiveFlux(flux, fluxVars);
        asImp_()->computeDiffusiveFlux(flux, fluxVars);
    }

    /*!
         * \brief Evaluates the advective mass flux of all components over
         *        a face of a subcontrol volume.
         *
         * \param flux The advective flux over the sub-control-volume face for each component
         * \param vars The flux variables at the current SCV
         */
    void computeAdvectiveFlux(PrimaryVariables &flux, const FluxVariables &fluxVars) const
    {
        ////////
        // advective fluxes of all components in all phases
        ////////

        // data attached to upstream and the downstream vertices
       // of the current phase
       const VolumeVariables &up =
           this->curVolVars_(fluxVars.upstreamIdx());
       const VolumeVariables &dn =
           this->curVolVars_(fluxVars.downstreamIdx());

        // total mass flux
       //KmvpNormal is the Darcy velocity multiplied with the normal vector, calculated in 1p2cfluxvariables.hh
        flux[contiEqIdx] =
            fluxVars.KmvpNormal() *
            ((     upwindAlpha)*up.density()/up.viscosity()
             +
             ((1 - upwindAlpha)*dn.density()/dn.viscosity()));

        // advective flux of the second component -molefraction
//       flux[transEqIdx] +=
//           fluxVars.KmvpNormal() *
//           ((    upwindAlpha)*up.concentration(comp1Idx)/up.viscosity()
//            +
//            (1 - upwindAlpha)*dn.concentration(comp1Idx)/dn.viscosity());

        // advective flux of the second component - massfraction
        flux[transEqIdx] +=
            fluxVars.KmvpNormal() *
            ((    upwindAlpha)*up.density() * up.massFrac(comp1Idx)/up.viscosity()
             +
             (1 - upwindAlpha)*dn.density()*dn.massFrac(comp1Idx)/dn.viscosity());

    }

    /*!
         * \brief Adds the diffusive mass flux of all components over
         *        a face of a subcontrol volume.
         *
         * \param flux The diffusive flux over the sub-control-volume face for each component
         * \param vars The flux variables at the current SCV
         */
    void computeDiffusiveFlux(PrimaryVariables &flux, const FluxVariables &fluxVars) const
    {
        Scalar tmp(0);
        // diffusive flux of second component
      tmp -=
           fluxVars.porousDiffCoeff() *
           (fluxVars.concentrationGrad(comp1Idx) * fluxVars.face().normal);

       // dispersive flux of second component
       Vector normalDisp;
       fluxVars.dispersionTensor().mv(fluxVars.face().normal, normalDisp);
       tmp -=
           (normalDisp * fluxVars.concentrationGrad(comp1Idx));

       //molar
//       flux[transEqIdx] += tmp;
       //transform to mass fractions
       flux[transEqIdx] += tmp * FluidSystem::molarMass(comp1Idx);
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
    
    Implementation *asImp_()
    { return static_cast<Implementation *> (this); }
    const Implementation *asImp_() const
    { return static_cast<const Implementation *> (this); }
};

}

#endif
