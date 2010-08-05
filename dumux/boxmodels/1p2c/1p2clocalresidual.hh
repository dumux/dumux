// $Id: 1p2clocalresidual.hh 3784 2010-06-24 13:43:57Z bernd $
/*****************************************************************************
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
 * \brief Calculate the local Jacobian for the single-phase,
 *        two-component model in the BOX scheme.
 */
template<class TypeTag>
class OnePTwoCLocalResidual : public BoxLocalResidual<TypeTag>
{
protected:
    typedef OnePTwoCLocalResidual<TypeTag> ThisType;
    typedef BoxLocalResidual<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePTwoCIndices)) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    enum
    {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),
        numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)),

        konti = Indices::konti,
        transport = Indices::transport,
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;


    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementSolutionVector)) ElementSolutionVector;

    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVariables)) FluxVariables;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef Dune::FieldMatrix<Scalar, dim, dim> Tensor;

    static const Scalar upwindAlpha = GET_PROP_VALUE(TypeTag, PTAG(UpwindAlpha));

public:
    /*!
     * \brief Evaluate the amount all conservation quantites
     *        (e.g. phase mass) within a finite volume.
     */
    void computeStorage(PrimaryVariables &result, int scvIdx, bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const ElementVolumeVariables &elemDat = usePrevSol ? this->prevVolVars_() : this->curVolVars_();
        const VolumeVariables &vertDat = elemDat[scvIdx];

        // storage term of continuity equation
        result[konti] = 0;

        // storage term of the transport equation
        result[transport] = vertDat.molarDensity * vertDat.porosity * vertDat.molefraction;
    }

    /*!
     * \brief Evaluates the mass flux over a face of a subcontrol
     *        volume.
     */
    void computeFlux(PrimaryVariables &flux, int faceId) const
    {
        FluxVariables vars(this->problem_(),
                      this->elem_(),
                      this->fvElemGeom_(),
                      faceId,
                      this->curVolVars_());
        flux = 0;

        // data attached to upstream and the downstream vertices
        // of the current phase
        const VolumeVariables &up = this->curVolVars_(vars.upstreamIdx);
        const VolumeVariables &dn = this->curVolVars_(vars.downstreamIdx);

        flux[konti] = vars.vDarcyNormal / vars.viscosityAtIP;

        // advective flux
        flux[transport] +=
        vars.vDarcyNormal *
        ( upwindAlpha*
                ( up.molarDensity * up.molefraction/up.viscosity )
                +
                (1 - upwindAlpha)*
                ( dn.molarDensity * dn.molefraction/dn.viscosity ) );

        Dune::FieldVector<Scalar,dim> unitNormal(vars.face->normal);
        unitNormal/=vars.face->normal.two_norm();

        // diffusive flux
        flux[transport] +=
        vars.molarDensityAtIP * vars.diffCoeffPM *
        (vars.concentrationGrad * vars.face->normal);

        //multiply hydrodynamic dispersion tensor with the face normal
        Dune::FieldVector<Scalar,dim> normalDisp;
        vars.dispersionTensor.mv(vars.face->normal, normalDisp);

        //add dispersive flux
        flux[transport] +=
        vars.molarDensityAtIP * (normalDisp * vars.concentrationGrad);
    }

    /*!
     * \brief Calculate the source term of the equation
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
