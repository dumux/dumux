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
#ifndef DUMUX_1P_BOX_MODEL_HH
#define DUMUX_1P_BOX_MODEL_HH

#include <dumux/new_models/boxscheme/boxscheme.hh>
#include <dumux/new_models/boxscheme/p1boxtraits.hh>

#include <dumux/auxiliary/apis.hh>

namespace Dune
{
    ///////////////////////////////////////////////////////////////////////////
    // traits for the single phase isothermal model
    ///////////////////////////////////////////////////////////////////////////
    /*!
     * \brief The pw-Sn specific traits.
     */
    class OnePTraits
    {
    public:
        enum {
            numEq     = 1,        //!< Number of primary variables
            numPhases = 1,    //!< Number of fluid phases
        };
        enum {
            pIdx = 0  //!< Index for the fluid pressure in a field vector
        };
    };
    
    ///////////////////////////////////////////////////////////////////////////
    // OnePBoxJacobian (evaluate the local jacobian for the newton method.)
    ///////////////////////////////////////////////////////////////////////////
    /*!
     * \brief Local Jacobian for the single phase isothermal model
     */
    template<class ProblemT, class BoxTraitsT, class OnePTraitsT>
    class OnePBoxJacobian : public BoxJacobian<ProblemT, 
                                               BoxTraitsT, 
                                               OnePBoxJacobian<ProblemT,
                                                               BoxTraitsT,
                                                               OnePTraitsT> >
    {
    private:
        typedef OnePBoxJacobian<ProblemT, BoxTraitsT, OnePTraitsT>  ThisType;
        typedef BoxJacobian<ProblemT, BoxTraitsT, ThisType >        ParentType;

        typedef ProblemT                                Problem;
        typedef typename Problem::DomainTraits          DomTraits;
        typedef BoxTraitsT                              BoxTraits;
        typedef OnePTraitsT                             OnePTraits;

        enum {
            dim            = DomTraits::dim,
            dimWorld       = DomTraits::dimWorld,

            numEq          = BoxTraits::numEq,

            numPhases      = OnePTraits::numPhases,
            pIdx           = OnePTraits::pIdx,
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

        typedef FieldMatrix<Scalar, dim, dim>  Tensor;

        /*!
         * \brief Data which is attached to each vert of the and can
         *        be shared between multiple calculations and should
         *        thus be cached in order to increase efficency.
         */
        struct VariableVertexData
        {
            Scalar density;
            Scalar viscosity;
        };

        /*!
         * \brief Cached data for the each vert of the element.
         */
        struct ElementData
        {
            VariableVertexData  vertex[BoxTraits::ShapeFunctionSetContainer::maxsize];
        };

    public:
        OnePBoxJacobian(ProblemT &problem)
            : ParentType(problem)
        {};

        /*!
         * \brief Set the parameters for the calls to the remaining
         *        members.
         */
        void setParams(const Element &element, LocalFunction &curSol, LocalFunction &prevSol)
        {
            setCurrentElement(element);

            curSol_ = &curSol;
            updateElementData_(curElemDat_, *curSol_);
            curSolDeflected_ = false;

            prevSol_ = &prevSol;
            updateElementData_(prevElemDat_, *prevSol_);
        };

        /*!
         * \brief Vary a single component of a single vert of the
         *        local solution for the current element.
         *
         * This method is a optimization, since if varying a single
         * component at a degree of freedom not the whole element cache
         * needs to be recalculated. (Updating the element cache is very
         * expensive since material laws need to be evaluated.)
         */
        void deflectCurSolution(int vert, int component, Scalar value)
        {
            // make sure that the orignal state can be restored
            if (!curSolDeflected_) {
                curSolDeflected_ = true;

                curSolOrigValue_ = (*curSol_)[vert][component];
                curSolOrigVarData_ = curElemDat_.vertex[vert];
            }

            (*curSol_)[vert][component] = value;
            updateVarVertexData_(curElemDat_.vertex[vert],
                                 (*curSol_)[vert],
                                 vert,
                                 ParentType::problem_.vertexIdx(ParentType::curElement_(),
                                                                vert));

        }

        /*!
         * \brief Restore the local jacobian to the state before
         *        deflectCurSolution() was called.
         *
         * This only works if deflectSolution was only called with
         * (vert, component) as arguments.
         */
        void restoreCurSolution(int vert, int component)
        {
            curSolDeflected_ = false;
            (*curSol_)[vert][component] = curSolOrigValue_;
            curElemDat_.vertex[vert] = curSolOrigVarData_;
        };

        /*!
         * \brief Evaluate the rate of change of all conservation
         *        quantites (e.g. phase mass) within a sub control
         *        volume of a finite volume element for the OneP 
         *        model.
         *
         * This function should not include the source and sink terms.
         */
        void computeStorage(SolutionVector &result, int scvId, bool usePrevSol) const
        {
//            LocalFunction *sol = usePrevSol?prevSol_:curSol_;
            const VariableVertexData &vertDat = usePrevSol?prevElemDat_.vertex[scvId]:curElemDat_.vertex[scvId];

            // partial time derivative of the wetting phase mass
            result[pIdx] = - vertDat.density
                * ParentType::problem_.porosity(*this->curElementPtr_, scvId);
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
            Scalar viscosityIJ = 0;
            LocalPosition pGrad(0);
            for (int k = 0; k < ParentType::curElementGeom_.numVertices; k++) {
                LocalPosition grad(face.grad[k]);
                grad *= (*curSol_)[k][pIdx];
                pGrad += grad;
                    
                densityIJ += curElemDat_.vertex[k].density*face.shapeValue[k];
                viscosityIJ += curElemDat_.vertex[k].viscosity*face.shapeValue[k];
            }

            // adjust pressure gradient by gravity force
            LocalPosition gravity = ParentType::problem_.gravity();
            gravity *= densityIJ;
            pGrad   -= gravity;

            // calculate darcy velocity
            Tensor K         = this->problem_.soil().K(global_i, ParentType::curElement_(), local_i);
            const Tensor &Kj = this->problem_.soil().K(global_j, ParentType::curElement_(), local_j);
            harmonicMeanK_(K, Kj);

            // temporary vector for the Darcy velocity
            GlobalPosition vDarcy;

            K.mv(pGrad, vDarcy);  // vDarcy = K * grad p
            vDarcy /= viscosityIJ;
               
            flux[pIdx] = vDarcy*normal;
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
        void addVtkFields(MultiWriter &writer, const SpatialFunction &globalSol) const
            {
                typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

                // create the required scalar fields
                unsigned numVertices = this->problem_.numVertices();
                ScalarField *p = writer.template createField<Scalar, 1>(numVertices);

                VariableVertexData tmp;
                ElementIterator it = this->problem_.elementBegin();
                ElementIterator endit = this->problem_.elementEnd();
                for (; it != endit; ++it) {
                    for (int i = 0; i < it->template count<dim>(); ++i) {
                        int globalI = this->problem_.vertexIdx(*it, i);

                        asImp_().updateVarVertexData_(tmp,
                                                      (*globalSol)[globalI],
                                                      i,
                                                      globalI);

                        (*p)[globalI] = (*globalSol)[globalI][pIdx];
                    };
                }

                writer.addVertexData(p, "p");
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

        /*!
         * \brief Pre-compute the element cache data.
         *
         * This method is called by BoxJacobian (which in turn is
         * called by the operator assembler) every time the current
         * element changes.
         */
        void updateElementData_(ElementData &dest, const LocalFunction &sol)
        {
            for (int i = 0; i < ParentType::curElementGeom_.numVertices; i++) {
                int iGlobal = ParentType::problem_.vertexIdx(ParentType::curElement_(), i);
                updateVarVertexData_(dest.vertex[i],
                                     sol[i],
                                     i,        // index of sub volume to update,
                                     iGlobal); // global vert index of the sub volume's grid vert
            }
        }


        void updateVarVertexData_(VariableVertexData   &vertexData,
                                  const SolutionVector &sol,
                                  int                  i, // local index of the subvolume/grid vertex
                                  int                  iGlobal) const // global index of the sub-volume's vertex
        {
//            const GlobalPosition &global = this->curElementGeom_.subContVol[i].global;
//            const LocalPosition &local = this->curElementGeom_.subContVol[i].local;

            vertexData.density = this->problem_.fluid().density(asImp_().temperature_(), 
                                                                sol[pIdx]);
            vertexData.viscosity = this->problem_.fluid().density(asImp_().temperature_(), 
                                                                  sol[pIdx]);
        }
        
        Scalar temperature_() const
        { return this->problem_.temperature(); }

        ThisType &asImp_() 
        { return *static_cast<ThisType *>(this); }

        const ThisType &asImp_() const
        { return *static_cast<const ThisType *>(this); }

        LocalFunction   *curSol_;
        ElementData      curElemDat_;

        bool               curSolDeflected_;
        Scalar             curSolOrigValue_;
        VariableVertexData curSolOrigVarData_;

        LocalFunction    *prevSol_;
        ElementData       prevElemDat_;
    };


    ///////////////////////////////////////////////////////////////////////////
    // OnePBoxModel (The actual numerical model.)
    ///////////////////////////////////////////////////////////////////////////
    /*!
     * \brief Adaption of the BOX scheme to the single phase isothermal flow model.
     */
    template<class ProblemT>
    class OnePBoxModel : public BoxScheme< // The implementation of the model
        OnePBoxModel<ProblemT>,

        // The Traits for the BOX method
        P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                    typename ProblemT::DomainTraits::Grid,
                    OnePTraits::numEq>,

        // The actual problem we would like to solve
        ProblemT,

        // The local jacobian operator
        OnePBoxJacobian<ProblemT,
                            P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                                        typename ProblemT::DomainTraits::Grid,
                                        OnePTraits::numEq>,
                            OnePTraits> >
    {
        typedef typename ProblemT::DomainTraits::Grid   Grid;
        typedef typename ProblemT::DomainTraits::Scalar Scalar;
        typedef OnePBoxModel<ProblemT>              ThisType;

    public:
        typedef P1BoxTraits<Scalar, Grid, OnePTraits::numEq> BoxTraits;
        typedef Dune::OnePTraits                             OnePTraits;

    private:
        typedef OnePBoxJacobian<ProblemT, BoxTraits, OnePTraits>  OnePLocalJacobian;
        typedef BoxScheme<ThisType,
                          BoxTraits,
                          ProblemT,
                          OnePLocalJacobian>  ParentType;

    public:
        typedef NewNewtonMethod<ThisType> NewtonMethod;

        OnePBoxModel(ProblemT &prob)
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
        void addVtkFields(MultiWriter &writer) const
            {
                richardsLocalJacobian_.addVtkFields(writer, this->currentSolution());
            }

    private:
        // calculates the jacobian matrix at a given position
        OnePLocalJacobian  richardsLocalJacobian_;
    };
}

#endif
