/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUMUX_RICHARDS_BOX_MODEL_HH
#define DUMUX_RICHARDS_BOX_MODEL_HH

#include <dumux/new_models/boxscheme/boxscheme.hh>
#include <dumux/new_models/boxscheme/p1boxtraits.hh>

#include <dumux/auxiliary/math.hh>

#include "richardstraits.hh"
#include "richardsvertexdata.hh"

#include <dumux/auxiliary/apis.hh>

namespace Dune
{

///////////////////////////////////////////////////////////////////////////
// RichardsBoxJacobian (evaluate the local jacobian for the newton method.)
///////////////////////////////////////////////////////////////////////////
/*!
 * \brief Local Jacobian for the single phase isothermal model
 */
template<class ProblemT, class BoxTraitsT, class RichardsTraitsT>
class RichardsBoxJacobian : public BoxJacobian<ProblemT,
                                               BoxTraitsT,
                                               RichardsBoxJacobian<ProblemT, BoxTraitsT, RichardsTraitsT>,
                                               RichardsVertexData<RichardsTraitsT, ProblemT> >
{
private:
    typedef RichardsBoxJacobian<ProblemT, BoxTraitsT, RichardsTraitsT>  ThisType;
    typedef RichardsVertexData<RichardsTraitsT, ProblemT>  VertexData;
    typedef BoxJacobian<ProblemT, BoxTraitsT, ThisType, VertexData >  ParentType;

    typedef ProblemT                                Problem;
    typedef typename Problem::DomainTraits          DomTraits;
    typedef BoxTraitsT                              BoxTraits;
    typedef RichardsTraitsT                         RichardsTraits;

    enum {
        dim            = DomTraits::dim,
        dimWorld       = DomTraits::dimWorld,

        numEq          = BoxTraits::numEq,

        numPhases      = RichardsTraits::numPhases,
        pWIdx          = RichardsTraits::pWIdx,
    };

    typedef typename DomTraits::Scalar              Scalar;
    typedef typename DomTraits::CoordScalar         CoordScalar;
    typedef typename DomTraits::Grid                Grid;
    typedef typename DomTraits::Element             Element;
    typedef typename DomTraits::ElementIterator     ElementIterator;
    typedef typename Element::EntityPointer         ElementPointer;
    typedef typename DomTraits::LocalPosition       LocalPosition;
    typedef typename DomTraits::GlobalPosition      GlobalPosition;

    typedef typename BoxTraits::FVElementGeometry   FVElementGeometry;
    typedef typename BoxTraits::SpatialFunction     SpatialFunction;
    typedef typename BoxTraits::LocalFunction       LocalFunction;
    typedef typename BoxTraits::SolutionVector      SolutionVector;

    typedef std::vector<VertexData>        VertexDataArray;
    typedef FieldMatrix<Scalar, dim, dim>  Tensor;

public:
    RichardsBoxJacobian(ProblemT &problem)
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
    void computeStorage(SolutionVector &result, int scvIdx, bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const VertexDataArray &elemDat = usePrevSol ? this->prevElemDat_  : this->curElemDat_;
        const VertexData  &vertDat = elemDat[scvIdx];

        // partial time derivative of the wetting phase mass
        result[pWIdx] = -vertDat.densityW
            * vertDat.porosity
            * this->prevElemDat_[scvIdx].dSwdpC // TODO: use derivative for the current solution
            * vertDat.pW;
    }


    /*!
     * \brief Evaluates the mass flux over a face of a subcontrol
     *        volume.
     */
    void computeFlux(SolutionVector &flux, int faceId) const
    {
        const typename FVElementGeometry::SubControlVolumeFace
            &face = ParentType::curElementGeom_.subContVolFace[faceId];

        const int i = face.i;
        const int j = face.j;

        // normal vector, value of the area of the scvf
        const GlobalPosition &normal(this->curElementGeom_.subContVolFace[faceId].normal);

        // get global coordinates of verts i,j
        const GlobalPosition &global_i = this->curElementGeom_.subContVol[i].global;
        const GlobalPosition &global_j = this->curElementGeom_.subContVol[j].global;

        // get local coordinates of verts i,j
        const LocalPosition &local_i = this->curElementGeom_.subContVol[i].local;
        const LocalPosition &local_j = this->curElementGeom_.subContVol[j].local;

        // calculate FE gradient
        Scalar densityIJ = 0;
        LocalPosition pGrad(0);
        for (int k = 0; k < this->curElementGeom_.numVertices; k++) {
            LocalPosition grad(face.grad[k]);
            grad *= this->curElemDat_[k].pW;
            pGrad += grad;

            densityIJ += this->curElemDat_[k].densityW*face.shapeValue[k];
        }

        // adjust pressure gradient by gravity force
        LocalPosition gravity = ParentType::problem_.gravity();
        gravity *= densityIJ;
        pGrad   -= gravity;

        // calculate darcy velocity
        Tensor K         = this->problem_.soil().K(global_i, ParentType::curElement_(), local_i);
        const Tensor &Kj = this->problem_.soil().K(global_j, ParentType::curElement_(), local_j);
        harmonicMeanK_(K, Kj);

        // magnitute of darcy velocity of each phase projected
        // on the normal of the sub-control volume's face
        Scalar vDarcyOut;
        // temporary vector for the Darcy velocity
        GlobalPosition vDarcy;

        K.mv(pGrad, vDarcy);  // vDarcy = K * grad p
        vDarcyOut = vDarcy*normal;

        // calculate the flux using upwind
        if (vDarcyOut < 0)
            flux[pWIdx] =
            	this->curElemDat_[i].densityW
                * this->curElemDat_[i].mobilityW
                * vDarcyOut;
        else
            flux[pWIdx] =
            	this->curElemDat_[j].densityW
                * this->curElemDat_[j].mobilityW
                * vDarcyOut;
    }

    /*!
     * \brief Calculate the source term of the equation
     */
    void computeSource(SolutionVector &q, int localVertexIdx)
    {
        this->problem_.source(q,
                              this->curElement_(),
                              this->curElementGeom_,
                              localVertexIdx);
    }

    /*!
     * \brief All relevant primary and secondary of a given
     *        solution to an ouput writer.
     */
    template <class MultiWriter>
    void addVtkFields(MultiWriter &writer, const SpatialFunction &globalSol)
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

        LocalFunction tmpSol;
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
                (*pW)[globalIdx] = this->problem_.pNreference() + this->curElemDat_[i].pW;
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
    // harmonic mean of the permeability computed directly.  the
    // first parameter is used to store the result.
    static void harmonicMeanK_(Tensor &Ki, const Tensor &Kj)
    {
        for (int kx=0; kx < Tensor::rows; kx++){
            for (int ky=0; ky< Tensor::cols; ky++){
                if (Ki[kx][ky] != Kj[kx][ky]) {
                    Ki[kx][ky] = harmonicMean_(Ki[kx][ky], Kj[kx][ky]);
                }
            }
        }
    }

    // returns the harmonic mean of two scalars
    static Scalar harmonicMean_(Scalar x, Scalar y)
    {
        if (x == 0 || y == 0)
            return 0;
        return (2*x*y)/(x + y);
    };

    ThisType &asImp_()
    { return *static_cast<ThisType *>(this); }

    const ThisType &asImp_() const
    { return *static_cast<const ThisType *>(this); }

 };


///////////////////////////////////////////////////////////////////////////
// RichardsBoxModel (The actual numerical model.)
///////////////////////////////////////////////////////////////////////////
/*!
 * \brief Adaption of the BOX scheme to the single phase isothermal flow model.
 */
template<class ProblemT>
class RichardsBoxModel : public BoxScheme< // The implementation of the model
    RichardsBoxModel<ProblemT>,

    // The Traits for the BOX method
    P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                typename ProblemT::DomainTraits::Grid,
                RichardsTraits::numEq>,

    // The actual problem we would like to solve
    ProblemT,

    // The local jacobian operator
    RichardsBoxJacobian<ProblemT,
                        P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                                    typename ProblemT::DomainTraits::Grid,
                                    RichardsTraits::numEq>,
                        RichardsTraits> >
{
    typedef typename ProblemT::DomainTraits::Grid   Grid;
    typedef typename ProblemT::DomainTraits::Scalar Scalar;
    typedef RichardsBoxModel<ProblemT>              ThisType;

public:
    typedef P1BoxTraits<Scalar, Grid, RichardsTraits::numEq> BoxTraits;
    typedef Dune::RichardsTraits                             RichardsTraits;

private:
    typedef RichardsBoxJacobian<ProblemT, BoxTraits, RichardsTraits>  RichardsLocalJacobian;
    typedef BoxScheme<ThisType,
                      BoxTraits,
                      ProblemT,
                      RichardsLocalJacobian>  ParentType;

public:
    typedef NewNewtonMethod<ThisType> NewtonMethod;

    RichardsBoxModel(ProblemT &prob)
        : ParentType(prob, richardsLocalJacobian_),
          richardsLocalJacobian_(prob)
    {
        Api::require<Api::BasicDomainTraits, typename ProblemT::DomainTraits>();
    }

    /*!
     * \brief All relevant primary and secondary of the current
     *        solution to an ouput writer.
     */
    template <class MultiWriter>
    void addVtkFields(MultiWriter &writer)
    {
        richardsLocalJacobian_.addVtkFields(writer, this->currentSolution());
    }

private:
    // calculates the jacobian matrix at a given position
    RichardsLocalJacobian  richardsLocalJacobian_;
};
}

#endif
