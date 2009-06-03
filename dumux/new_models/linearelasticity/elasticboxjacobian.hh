//$Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Melanie Darcis                                    *
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
#ifndef DUMUX_ELASTIC_BOX_JACOBIAN_HH
#define DUMUX_ELASTIC_BOX_JACOBIAN_HH

#include <dumux/new_models/boxscheme/boxjacobian.hh>

#include <dune/istl/scaledidmatrix.hh>
#include <tr1/array>

#include "elasticvertexdata.hh"


namespace Dune
{

/*!
 * \brief Calculate the local Jacobian for the linear elasticity
 *        model.
 */
template<class TypeTag>
class ElasticBoxJacobian : public BoxJacobian<TypeTag, ElasticBoxJacobian<TypeTag>>
{
protected:
    typedef ElasticBoxJacobian<TypeTag>    ThisType;
    typedef BoxJacobian<TypeTag, ThisType> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))        Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElasticIndices)) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))       Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView))     GridView;

    enum {
        dim            = GridView::dimension,
        dimWorld       = GridView::dimensionworld,

        numEq          = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
    };

    typedef typename GridView::template Codim<0>::Entity   Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef FieldVector<Scalar, dim>       LocalPosition;
    typedef FieldVector<Scalar, dimWorld>  GlobalPosition;

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVarVector        PrimaryVarVector;
    typedef typename SolutionTypes::SolutionFunction        SolutionFunction;
    typedef typename SolutionTypes::SolutionOnElement       SolutionOnElement;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexData)) VertexData;
    typedef std::vector<VertexData>                           VertexDataArray;

    typedef Dune::FieldVector<Scalar, dim>      DisplacementVector;
    typedef Dune::FieldMatrix<Scalar, dim, dim> DisplacementGradient;
    typedef Dune::FieldMatrix<Scalar, dim, dim> StressTensor;
    typedef Dune::FieldMatrix<Scalar, dim, dim> StrainTensor;
    typedef Dune::ScaledIdentityMatrix<Scalar, dim> IdentityMatrix;

public:
    ElasticBoxJacobian(Problem &problem)
        : ParentType(problem)
    {
    };

    /*!
     * \brief Evaluate the rate of change of all conservation
     *        quantites (e.g. phase mass) within a sub control
     *        volume of a finite volume element in the 2P-2C
     *        model.
     *
     * This function should not include the source and sink terms.
     */
    void computeStorage(PrimaryVarVector &result, int scvIdx, bool usePrevSol) const
    {
        result = Scalar(0);

        // if flag usePrevSol is set, the solution from the previous time step is used,
        // otherwise the current solution is used. The secondary variables are used accordingly.
        // This computes the derivative of the storage term.
        //const VertexDataArray &elemDat = usePrevSol ? this->prevElemDat_  : this->curElemDat_;
        //const VertexData  &vertDat = elemDat[scvIdx];

        // quasistationary conditions assumed
        result = 0.0;
    }

    /*!
     * \brief Evaluates the mass flux over a face of a subcontrol
     *        volume.
     */
    void computeFlux(PrimaryVarVector &flux, int faceId) const
    {
        // set flux vector to zero
        int i = this->curElementGeom_.subContVolFace[faceId].i;
        int j = this->curElementGeom_.subContVolFace[faceId].j;

        // normal vector, value of the area of the scvf
        const GlobalPosition &normal(this->curElementGeom_.subContVolFace[faceId].normal);

        const VertexData &vDat_i = this->curElemDat_[i];
        const VertexData &vDat_j = this->curElemDat_[j];

        Scalar epsilonTimesIdentity(0.0), divU(0.0);
        GlobalPosition tmp(0.0);
        DisplacementGradient gradU(0.0), gradUTransposed(0.0);
        StrainTensor epsilon(0.0);
        IdentityMatrix identity;
        StressTensor sigma(0.0), firstTerm(0.0), secondTerm(0.0);

        // calculate FE gradient (grad p for each phase)
        for (int k = 0; k < this->curElementGeom_.numVertices; k++) // loop over adjacent vertices

        {
            // FEGradient at vertex k
            const LocalPosition &feGrad = this->curElementGeom_.subContVolFace[faceId].grad[k];

            // compute sum of displacement gradients for each space direction
            for (int coordIdx = 0; coordIdx < dim; ++coordIdx) {
                tmp = feGrad;
                tmp *= this->curElemDat_[k].displacement[coordIdx];
                gradU[coordIdx] += tmp;
            }
        }

        divU = 0;
        for(int col=0; col < dim; col++)
        {
            divU += gradU[col][col];
            for(int row=0; row<dim; row++)
                gradUTransposed[row][col] = gradU[col][row];
        }

        epsilon += gradU;
        epsilon += gradUTransposed;
        epsilon *= 0.5;
        epsilonTimesIdentity = divU;

        // get Lame parameters
        Scalar lambda = 0.5*(vDat_i.lambda + vDat_j.lambda); //   CHECK which kind of mean is appropriate??
        Scalar mu = 0.5*(vDat_i.mu + vDat_j.mu); //   CHECK which kind of mean is appropriate??


        // displacement
        firstTerm += epsilon;
        firstTerm *= 2;
        firstTerm *= mu;

        for (int i = 0; i < dim; ++i)
            secondTerm[i][i] += 1.0;
        secondTerm *= lambda;
        secondTerm *= epsilonTimesIdentity;

        sigma += firstTerm;
        sigma += secondTerm;

        sigma.mv(normal, tmp);

        // TODO assign the tensions to the primary variables
        for (int i=0; i < dim; ++i)
            flux[Indices::disp2Primary(i)] = tmp[i];
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
     * \brief Add the mass fraction of air in water to VTK output of
     *        the current timestep.
     */
    template <class MultiWriter>
    void addVtkFields(MultiWriter &writer, const SolutionFunction &globalSol)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

        // create the required scalar fields
        unsigned numVertices = this->problem_.numVertices();
        ScalarField *ux = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *uy = (dim < 2)?NULL:writer.template createField<Scalar, 1>(numVertices);
        ScalarField *uz = (dim < 3)?NULL:writer.template createField<Scalar, 1>(numVertices);

        SolutionOnElement curSol;

        ElementIterator elementIt = this->problem_.elementBegin();
        ElementIterator endit = this->problem_.elementEnd();

        for (; elementIt != endit; ++elementIt)
        {
            int numLocalVerts = elementIt->template count<dim>();

            setCurrentElement(*elementIt);
            this->restrictToElement(curSol, globalSol);
            this->updateElementData_(this->curElemDat_, curSol, false);

            for (int i = 0; i < numLocalVerts; ++i)
            {
                int globalIdx = this->problem_.vertexIdx(*elementIt, i);

                (*ux)[globalIdx] = this->curElemDat_[i].displacement[0];
                if (dim >= 2)
                    (*uy)[globalIdx] = this->curElemDat_[i].displacement[1];
                if (dim >= 3)
                    (*uz)[globalIdx] = this->curElemDat_[i].displacement[2];
            };

        }

        writer.addVertexData(ux, "ux");
        if (dim >= 2)
            writer.addVertexData(uy, "uy");
        if (dim >= 3)
            writer.addVertexData(uz, "uz");
    }
};


}

#endif
