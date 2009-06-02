//$Id$

/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch, Melanie Darcis                   *
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
#ifndef DUMUX_NEW_ELASTIC_BOX_MODEL_HH
#define DUMUX_NEW_ELASTIC_BOX_MODEL_HH

#include <dumux/new_models/boxscheme/boxscheme.hh>
#include "elastictraits.hh"

#include <dumux/auxiliary/apis.hh>

#include <vector>

namespace Dune
{

///////////////////////////////////////////////////////////////////////////
// ElasticJacobian (evaluate the local jacobian for the newton method.)
///////////////////////////////////////////////////////////////////////////
/*!
 * \brief elastic deformation specific details needed to approximately calculate
 *        the local jacobian in the BOX scheme.
 *
 */
template<class ProblemT, class BoxTraitsT, class ElasticTraitsT,
        class Implementation>
class ElasticBoxJacobianBase: public BoxJacobian<ProblemT, BoxTraitsT,
        Implementation>
{
    protected:
        typedef ElasticBoxJacobianBase<ProblemT, BoxTraitsT, ElasticTraitsT,
                Implementation> ThisType;
        typedef BoxJacobian<ProblemT, BoxTraitsT, Implementation> ParentType;

        typedef ProblemT Problem;
typedef        typename Problem::DomainTraits DomTraits;
        typedef BoxTraitsT BoxTraits;
        typedef ElasticTraitsT ElasticTraits;

        enum
        {
            dim = DomTraits::dim,
            dimWorld = DomTraits::dimWorld,

            numEq = BoxTraits::numEq,

            uxIdx = ElasticTraits::uxIdx,
            uyIdx = ElasticTraits::uyIdx,
            uzIdx = ElasticTraits::uzIdx,
        };

        typedef typename DomTraits::Scalar Scalar;
        typedef typename DomTraits::CoordScalar CoordScalar;
        typedef typename DomTraits::Grid Grid;
        typedef typename DomTraits::Element Element;
        typedef typename DomTraits::ElementIterator ElementIterator;
        typedef typename Element::EntityPointer ElementPointer;
        typedef typename DomTraits::LocalPosition LocalPosition;
        typedef typename DomTraits::GlobalPosition GlobalPosition;
        typedef typename DomTraits::VertexIterator VertexIterator;

        typedef typename BoxTraits::SolutionVector SolutionVector;
        typedef typename BoxTraits::FVElementGeometry FVElementGeometry;
        typedef typename BoxTraits::SpatialFunction SpatialFunction;
        typedef typename BoxTraits::LocalFunction LocalFunction;

        typedef typename ElasticTraits::VariableVertexData VariableVertexData;
        typedef typename ElasticTraits::DisplacementVector SolidDisplacement;
        typedef typename ElasticTraits::DisplacementGradient DisplacementGradient;
        typedef typename ElasticTraits::StressTensor StressTensor;
        typedef typename ElasticTraits::StrainTensor StrainTensor;
        typedef typename ElasticTraits::IdentityMatrix IdentityMatrix;
        typedef FieldMatrix<Scalar, dim, dim> Tensor;

        /*!
         * \brief Cached data for the each vertex of the element.
         */
        struct ElementData
        {
            VariableVertexData vertex[BoxTraits::ShapeFunctionSetContainer::maxsize];
        };

        /*!
         * \brief Function to update variable data of the vertices of the
         *        the current element (essentially secondary variables)
         */
        void updateVarVertexData_(VariableVertexData &vertDat,
                const SolutionVector &vertSol,
                const Element &element,
                int localIdx,
                Problem &problem) const
        {
            const GlobalPosition &global = element.geometry().corner(localIdx);
            const LocalPosition &local =
            DomTraits::referenceElement(element.type()).position(localIdx,
                    dim);
            vertDat.ux = vertSol[uxIdx];
            vertDat.uy = vertSol[uyIdx];
            vertDat.uz = vertSol[uzIdx];
            vertDat.lambda = this->problem_.soil().lameParams(global, element, local)[0];
            vertDat.mu = this->problem_.soil().lameParams(global, element, local)[1];

        }

        public:
        ElasticBoxJacobianBase(ProblemT &problem)
        : ParentType(problem)
        {
        };

        /*!
         * \brief Set the current grid element.
         */
        void setCurrentElement(const Element &element)
        {
            ParentType::setCurrentElement_(element);
        };

        /*!
         * \brief Set the parameters for the calls to the remaining
         *        members.
         */
        void setParams(const Element &element, LocalFunction &curSol, LocalFunction &prevSol)
        {
            setCurrentElement(element);

            // TODO: scheme which allows not to copy curSol and
            // prevSol all the time
            curSol_ = curSol;
            updateElementData_(curElemDat_, curSol_, false);
            curSolDeflected_ = false;

            prevSol_ = prevSol;
            updateElementData_(prevElemDat_, prevSol_, true);
        };

        /*!
         * \brief Vary a single component of a single vertex of the
         *        local solution for the current element.
         *
         * This method is a optimization, since if varying a single
         * component at a degree of freedom not the whole element cache
         * needs to be recalculated. (Updating the element cache is very
         * expensive since material laws need to be evaluated.)
         */
        void deflectCurSolution(int vert, int component, Scalar value)
        {
            // make sure that the original state can be restored
            if (!curSolDeflected_)
            {
                curSolDeflected_ = true;

                curSolOrigValue_ = curSol_[vert][component];
                curSolOrigVarData_ = curElemDat_.vertex[vert];
            }

            int globalIdx = ParentType::problem_.vertexIdx(ParentType::curElement_(),
                    vert);

            curSol_[vert][component] = value;
            asImp_()->updateVarVertexData_(curElemDat_.vertex[vert],
                    curSol_[vert],
                    this->curElement_(),
                    vert,
                    this->problem_);
        }

        /*!
         * \brief Restore the local jacobian to the state before
         *        deflectCurSolution() was called.
         *
         * This only works if deflectSolution was only called with
         * (vertex, component) as arguments.
         */
        void restoreCurSolution(int vert, int component)
        {
            curSolDeflected_ = false;
            curSol_[vert][component] = curSolOrigValue_;
            curElemDat_.vertex[vert] = curSolOrigVarData_;
        };

        /*!
         * \brief Evaluate the rate of change of all conservation
         *        quantites (e.g. phase mass) within a sub control
         *        volume of a finite volume element in the 2P-2C
         *        model.
         *
         * This function should not include the source and sink terms.
         */
        void computeStorage(SolutionVector &result, int scvId, bool usePrevSol) const
        {
            result = Scalar(0);

            // if flag usePrevSol is set, the solution from the previous time step is used,
            // otherwise the current solution is used. The secondary variables are used accordingly.
            // This computes the derivative of the storage term.
            const LocalFunction &sol = usePrevSol ? this->prevSol_ : this->curSol_;
            const ElementData &elementCache = usePrevSol ? prevElemDat_ : curElemDat_;

            // quasistationary conditions assumed
            result = 0.0;
        }

        /*!
         * \brief Evaluates the mass flux over a face of a subcontrol
         *        volume.
         */
        void computeFlux(SolutionVector &flux, int faceId) const
        {
            // set flux vector to zero
            int i = this->curElementGeom_.subContVolFace[faceId].i;
            int j = this->curElementGeom_.subContVolFace[faceId].j;

            // normal vector, value of the area of the scvf
            const GlobalPosition &normal(this->curElementGeom_.subContVolFace[faceId].normal);

            const ElementData &elemDat = this->curElemDat_;
            const VariableVertexData &vDat_i = elemDat.vertex[i];
            const VariableVertexData &vDat_j = elemDat.vertex[j];

            Scalar epsilonTimesIdentity(0.0), divU(0.0);
            GlobalPosition tmp(0.0);
            SolidDisplacement u(0.0);
            DisplacementGradient gradU(0.0), gradUTransposed(0.0);
            StrainTensor epsilon(0.0);
            IdentityMatrix identity(0.0);
            StressTensor sigma(0.0), firstTerm(0.0), secondTerm(0.0);

            // calculate FE gradient (grad p for each phase)
            for (int k = 0; k < this->curElementGeom_.numVertices; k++) // loop over adjacent vertices

            {
                // FEGradient at vertex k
                const LocalPosition &feGrad = this->curElementGeom_.subContVolFace[faceId].grad[k];

                u[uxIdx] = elemDat.vertex[k].ux;
                u[uyIdx] = elemDat.vertex[k].uy;
                u[uzIdx] = elemDat.vertex[k].uz;
                // compute sum of displacement gradients for each space direction

                // the displacement gradients
                tmp = feGrad;
                tmp *= u[uxIdx];
                gradU[0] += tmp;

                tmp = feGrad;
                tmp *= u[uyIdx];
                gradU[1] += tmp;

                tmp = feGrad;
                tmp *= u[uzIdx];
                gradU[2] += tmp;

            }

            divU = gradU[0][0] + gradU[1][1] + gradU[2][2];

            for(int cols=0; cols<3; cols++)
            {
                for(int rows=0; rows<3; rows++)
                {
                    gradUTransposed[rows][cols] = gradU[cols][rows];
                    identity[rows][rows] = 1.0;
                }
            }

            epsilon += gradU;
            epsilon += gradUTransposed;
            epsilon *= 0.5;
            epsilonTimesIdentity = divU;

            // get Lame parameters
            Scalar lambda = 0.5*(vDat_i.lambda + vDat_j.lambda); //   CHECK which kind of mean is appropriate??
            Scalar mu = 0.5*(vDat_i.mu + vDat_j.mu); //   CHECK which kind of mean is appropriate??


            // flux in the wetting phase
            firstTerm += epsilon;
            firstTerm *= 2;
            firstTerm *= mu;

            secondTerm += identity;
            secondTerm *= lambda;
            secondTerm *= epsilonTimesIdentity;

            sigma += firstTerm;
            sigma += secondTerm;
            sigma.mv(normal, flux);

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
         * \brief Add the mass fraction of air in water to VTK output of
         *        the current timestep.
         */
        template <class MultiWriter>
        void addVtkFields(MultiWriter &writer, const SpatialFunction &globalSol)
        {
            typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

            // create the required scalar fields
            unsigned numVertices = this->problem_.numVertices();
            ScalarField *ux = writer.template createField<Scalar, 1>(numVertices);
            ScalarField *uy = writer.template createField<Scalar, 1>(numVertices);
            ScalarField *uz = writer.template createField<Scalar, 1>(numVertices);

            LocalFunction curSol(numVertices);
            ElementData elemDat;

            ElementIterator elementIt = this->problem_.elementBegin();
            ElementIterator endit = this->problem_.elementEnd();

            for (; elementIt != endit; ++elementIt)
            {
                int numLocalVerts = elementIt->template count<dim>();

                setCurrentElement(*elementIt);
                this->restrictToElement(curSol, globalSol);
                updateElementData_(elemDat, curSol, false);

                for (int i = 0; i < numLocalVerts; ++i)
                {
                    int globalIdx = this->problem_.vertexIdx(*elementIt, i);

                    (*ux)[globalIdx] = elemDat.vertex[i].ux;
                    (*uy)[globalIdx] = elemDat.vertex[i].uy;
                    (*uz)[globalIdx] = elemDat.vertex[i].uz;
                };

            }

            writer.addVertexData(ux, "ux");
            writer.addVertexData(uy, "uy");
            writer.addVertexData(uz, "uz");
        }

        protected:
        Implementation *asImp_()
        {   return static_cast<Implementation *>(this);}
        const Implementation *asImp_() const
        {   return static_cast<const Implementation *>(this);}

        void updateElementData_(ElementData &dest, const LocalFunction &sol, bool isOldSol)
        {
            int phaseState;
            int numVertices = this->curElement_().template count<dim>();
            for (int i = 0; i < numVertices; i++)
            {
                int iGlobal = ParentType::problem_.vertexIdx(ParentType::curElement_(), i);
                asImp_()->updateVarVertexData_(dest.vertex[i],
                        sol[i],
                        this->curElement_(),
                        i,
                        this->problem_);
            }
        }

        // returns the harmonic mean of two scalars
        static Scalar harmonicMean_(Scalar x, Scalar y)
        {
            if (x == 0 || y == 0)
            return 0;
            return (2*x*y)/(x + y);
        };

        // parameters given in constructor
        bool switchFlag_;
        int formulaion_;

        // current solution
        LocalFunction curSol_;
        ElementData curElemDat_;

        // needed for restoreCurSolution()
        bool curSolDeflected_;
        Scalar curSolOrigValue_;
        VariableVertexData curSolOrigVarData_;

        // previous solution
        LocalFunction prevSol_;
        ElementData prevElemDat_;
    };

    template<class ProblemT,
    class BoxTraitsT,
    class ElasticTraitsT>
    class ElasticBoxJacobian : public ElasticBoxJacobianBase<ProblemT,
    BoxTraitsT,
    ElasticTraitsT,
    ElasticBoxJacobian<ProblemT,
    BoxTraitsT,
    ElasticTraitsT> >
    {
        typedef ElasticBoxJacobian<ProblemT,
        BoxTraitsT,
        ElasticTraitsT> ThisType;
        typedef ElasticBoxJacobianBase<ProblemT,
        BoxTraitsT,
        ElasticTraitsT,
        ThisType> ParentType;

        typedef ProblemT Problem;

        typedef typename Problem::DomainTraits DomTraits;
        typedef ElasticTraitsT ElasticTraits;
        typedef BoxTraitsT BoxTraits;

        typedef typename DomTraits::Scalar Scalar;
        typedef typename DomTraits::Vertex Vertex;

        typedef typename DomTraits::GlobalPosition GlobalPosition;
        typedef typename DomTraits::LocalPosition LocalPosition;

        typedef typename BoxTraits::SolutionVector SolutionVector;
        typedef typename BoxTraits::FVElementGeometry FVElementGeometry;

        typedef typename ParentType::VariableVertexData VariableVertexData;
        typedef typename ParentType::LocalFunction LocalFunction;
        typedef typename ParentType::ElementData ElementData;

        public:
        ElasticBoxJacobian(ProblemT &problem)
        : ParentType(problem)
        {
        };

    };

    ///////////////////////////////////////////////////////////////////////////
    // ElasticBoxModel (The actual numerical model.)
    ///////////////////////////////////////////////////////////////////////////
    /**
     * \brief Isothermal two phase two component model with Pw and
     *        Sn/X as primary unknowns.
     *
     * This implements an isothermal two phase two component model
     * with Pw and Sn/X as primary unknowns.
     */

    template<class ProblemT, class ElasticTraitsT = ElasticTraits<typename ProblemT::DomainTraits::Scalar> >
    class ElasticBoxModel
    : public BoxScheme<ElasticBoxModel<ProblemT, ElasticTraitsT>, // Implementation of the box scheme

    // The Traits for the BOX method
    P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
    typename ProblemT::DomainTraits::Grid,
    ElasticTraitsT::numEq>,

    // The actual problem we would like to solve
    ProblemT,
    // The local jacobian operator
    ElasticBoxJacobian<ProblemT,
    P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
    typename ProblemT::DomainTraits::Grid,
    ElasticTraitsT::numEq>, ElasticTraitsT> >
    {
        typedef typename ProblemT::DomainTraits::Grid Grid;
        typedef typename ProblemT::DomainTraits::Scalar Scalar;
        typedef ElasticBoxModel<ProblemT,ElasticTraitsT> ThisType;

        public:
        typedef ElasticTraitsT ElasticTraits;
        typedef P1BoxTraits<Scalar, Grid, ElasticTraits::numEq> BoxTraits;

        private:
        typedef ElasticBoxJacobian<ProblemT, BoxTraits, ElasticTraits> ElasticLocalJacobian;
        typedef BoxScheme<ThisType,
        BoxTraits,
        ProblemT,
        ElasticLocalJacobian> ParentType;

        typedef typename ProblemT::DomainTraits DomTraits;
        typedef typename DomTraits::Element Element;
        typedef typename DomTraits::ElementIterator ElementIterator;
        typedef typename DomTraits::LocalPosition LocalPosition;
        typedef typename DomTraits::GlobalPosition GlobalPosition;

        enum
        {
            dim = DomTraits::dim,
            dimWorld = DomTraits::dimWorld
        };

        public:
        typedef NewNewtonMethod<ThisType> NewtonMethod;

        ElasticBoxModel(ProblemT &prob)
        : ParentType(prob, elasticLocalJacobian_),
        elasticLocalJacobian_(prob)
        {
            Api::require<Api::BasicDomainTraits, typename ProblemT::DomainTraits>();
        }

        /*!
         * \brief Called by the update() method if applying the newton
         *         method was unsuccessful.
         */
        void updateFailedTry()
        {
            ParentType::updateFailedTry();
        };

        /*!
         * \brief Called by the BoxScheme's update method.
         */
        void updateSuccessful()
        {
            ParentType::updateSuccessful();
        }

        /*!
         * \brief Add the mass fraction of air in water to VTK output of
         *        the current timestep.
         */
        template <class MultiWriter>
        void addVtkFields(MultiWriter &writer)
        {
            elasticLocalJacobian_.addVtkFields(writer, this->currentSolution());
        }

        private:
        // calculates the jacobian matrix at a given position
        ElasticLocalJacobian elasticLocalJacobian_;
    };
}

#endif
