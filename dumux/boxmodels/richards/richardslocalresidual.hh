// $Id: richardslocalresidual.hh 3738 2010-06-15 14:01:09Z lauser $
/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
 *   Copyright (C) 2009 by Onur Dogan                                        *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
#ifndef DUMUX_RICHARDS_BOX_JACOBIAN_HH
#define DUMUX_RICHARDS_BOX_JACOBIAN_HH

#include <dumux/boxmodels/common/boxlocalresidual.hh>

#include "richardssecondaryvars.hh"

#include "richardsfluxvars.hh"

namespace Dumux
{
/*!
 * \ingroup RichardsModel
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the Richards box model.
 */
template<class TypeTag>
class RichardsLocalResidual : public BoxLocalResidual<TypeTag>
{
    typedef RichardsLocalResidual<TypeTag> ThisType;
    typedef BoxLocalResidual<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;


    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementSolutionVector)) ElementSolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVarVector)) PrimaryVarVector;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(RichardsIndices)) Indices;

    enum {
        dim        = GridView::dimension,
        dimWorld   = GridView::dimensionworld,

        pW         = Indices::pW,
    };


    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SecondaryVars)) SecondaryVars;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVars)) FluxVars;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementSecondaryVars)) ElementSecondaryVars;

    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    static const Scalar mobilityUpwindAlpha = GET_PROP_VALUE(TypeTag, PTAG(MobilityUpwindAlpha));

public:
    /*!
     * \brief Evaluate the rate of change of all conservation
     *        quantites (e.g. phase mass) within a sub control
     *        volume of a finite volume element for the Richards
     *        model.
     *
     * This function should not include the source and sink terms.
     */
    void computeStorage(PrimaryVarVector &result, int scvIdx, bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const ElementSecondaryVars &elemDat = usePrevSol ? this->prevSecVars_()  : this->curSecVars_();
        const SecondaryVars  &vertDat = elemDat[scvIdx];

        // partial time derivative of the wetting phase mass
        result[pW] =
            vertDat.densityW
            * vertDat.porosity
            * this->prevSecVars_()[scvIdx].dSwdpC // TODO: use derivative for the current solution
            * (vertDat.pNreference - vertDat.pW);
    }


    /*!
     * \brief Evaluates the mass flux over a face of a subcontrol
     *        volume.
     */
    void computeFlux(PrimaryVarVector &flux, int faceId) const
    {
        FluxVars vars(this->problem_(),
                      this->elem_(),
                      this->fvElemGeom_(),
                      faceId,
                      this->curSecVars_());

        // data attached to upstream and the downstream vertices
        const SecondaryVars &up = this->curSecVars_(vars.upstreamIdx);
        const SecondaryVars &dn = this->curSecVars_(vars.downstreamIdx);

        flux[pW] =
            vars.vDarcyNormal*
            (  mobilityUpwindAlpha*
               (  up.densityW *
                  up.mobilityW)
               +
               (1 - mobilityUpwindAlpha)*
               (  dn.densityW*
                  dn.mobilityW));
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

    /*!
     * \brief Return the temperature given the solution vector of a
     *        finite volume.
     */
    template <class PrimaryVarVector>
    Scalar temperature(const PrimaryVarVector &sol)
    { return this->problem_.temperature(); /* constant temperature */ }

private:
    ThisType &asImp_()
    { return *static_cast<ThisType *>(this); }

    const ThisType &asImp_() const
    { return *static_cast<const ThisType *>(this); }
};

};

#endif
