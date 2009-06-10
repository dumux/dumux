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

#include <dumux/boxmodels/boxscheme/boxjacobian.hh>

#include "richardsvertexdata.hh"
#include "richardselementdata.hh"
#include "richardsfluxdata.hh"

namespace Dune
{
/*!
 * \ingroup RichardsBoxModel
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the Richards box model.
 */
template<class TypeTag>
class RichardsBoxJacobian : public BoxJacobian<TypeTag, RichardsBoxJacobian<TypeTag> >
{
    typedef RichardsBoxJacobian<TypeTag>         ThisType;
    typedef BoxJacobian<TypeTag, ThisType>   ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))   Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView))  GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))    Scalar;

    typedef typename GridView::template Codim<0>::Entity   Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::SolutionFunction        SolutionFunction;
    typedef typename SolutionTypes::SolutionOnElement       SolutionOnElement;
    typedef typename SolutionTypes::PrimaryVarVector        PrimaryVarVector;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(RichardsIndices)) Indices;

    enum {
        dim        = GridView::dimension,
        dimWorld   = GridView::dimensionworld,

        pWIdx      = Indices::pWIdx,
    };


    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexData))   VertexData;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxData))     FluxData;
    typedef std::vector<VertexData> VertexDataArray;

    typedef Dune::FieldVector<Scalar, dim>       LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld>  GlobalPosition;

    static const Scalar mobilityUpwindAlpha = GET_PROP_VALUE(TypeTag, PTAG(MobilityUpwindAlpha));

public:
    RichardsBoxJacobian(Problem &problem)
        : ParentType(problem)
    {};

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
        const VertexDataArray &elemDat = usePrevSol ? this->prevElemDat_  : this->curElemDat_;
        const VertexData  &vertDat = elemDat[scvIdx];

        // partial time derivative of the wetting phase mass
        result[pWIdx] =
            vertDat.densityW
            * vertDat.porosity
            * this->prevElemDat_[scvIdx].dSwdpC // TODO: use derivative for the current solution
            * (vertDat.pNreference - vertDat.pW);
    }


    /*!
     * \brief Evaluates the mass flux over a face of a subcontrol
     *        volume.
     */
    void computeFlux(PrimaryVarVector &flux, int faceId) const
    {
		FluxData vars(this->problem_,
                      this->curElement_(),
                      this->curElementGeom_,
                      faceId,
                      this->curElemDat_);

        // data attached to upstream and the downstream vertices
        const VertexData &up = this->curElemDat_[vars.upstreamIdx];
        const VertexData &dn = this->curElemDat_[vars.downstreamIdx];

        flux[pWIdx] =
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
        this->problem_.source(q,
                              this->curElement_(),
                              this->curElementGeom_,
                              localVertexIdx);
    }

    /*!
     * \brief Return the temperature given the solution vector of a
     *        finite volume.
     */
    template <class PrimaryVarVector>
    Scalar temperature(const PrimaryVarVector &sol)
    { return this->problem_.temperature(); /* constant temperature */ }

    /*!
     * \brief All relevant primary and secondary of a given
     *        solution to an ouput writer.
     */
    template <class MultiWriter>
    void addVtkFields(MultiWriter &writer, const SolutionFunction &globalSol)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

        // create the required scalar fields
        unsigned numVertices = this->problem_.numVertices();
        ScalarField *Sw =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *Sn =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *pC =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *pW =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *dSwdpC =       writer.template createField<Scalar, 1>(numVertices);
        ScalarField *rhoW =         writer.template createField<Scalar, 1>(numVertices);
        ScalarField *mobW =         writer.template createField<Scalar, 1>(numVertices);

        SolutionOnElement tmpSol;
        ElementIterator elementIt = this->problem_.elementBegin();
        ElementIterator endit = this->problem_.elementEnd();
        for (; elementIt != endit; ++elementIt)
        {
            int numLocalVerts = elementIt->template count<dim>();
            tmpSol.resize(numLocalVerts);

            setCurrentElement(*elementIt);
            this->restrictToElement(tmpSol, globalSol);
            this->setCurrentSolution(tmpSol);

            for (int i = 0; i < numLocalVerts; ++i)
            {
                int globalIdx = this->problem_.vertexIdx(*elementIt, i);

                (*Sw)[globalIdx] = this->curElemDat_[i].Sw;
                (*Sn)[globalIdx] = 1.0 - this->curElemDat_[i].Sw;
                (*pC)[globalIdx] = this->curElemDat_[i].pC;
                (*pW)[globalIdx] = this->curElemDat_[i].pW;
                (*dSwdpC)[globalIdx] = this->curElemDat_[i].dSwdpC;
                (*rhoW)[globalIdx] = this->curElemDat_[i].densityW;
                (*mobW)[globalIdx] = this->curElemDat_[i].mobilityW;
            };
        }


        writer.addVertexData(Sw, "Sw");
        writer.addVertexData(Sn, "Sn");
        writer.addVertexData(pC, "pC");
        writer.addVertexData(pW, "pW");
        writer.addVertexData(dSwdpC, "dSwdpC");
        writer.addVertexData(rhoW, "rhoW");
        writer.addVertexData(mobW, "mobW");
    }

private:
    ThisType &asImp_()
    { return *static_cast<ThisType *>(this); }

    const ThisType &asImp_() const
    { return *static_cast<const ThisType *>(this); }
};

};

#endif
