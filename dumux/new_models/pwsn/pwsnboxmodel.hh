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
#ifndef DUMUX_PWSN_BOX_MODEL_HH
#define DUMUX_PWSN_BOX_MODEL_HH

#include <dumux/new_models/boxscheme/boxscheme.hh>
#include <dumux/new_models/boxscheme/p1boxtraits.hh>

#include <dumux/auxiliary/apis.hh>

namespace Dune
{
    ///////////////////////////////////////////////////////////////////////////
    // PwSn traits (central place for names and indices required by the
    // PwSnBoxJacobian and PwSnBoxModel)
    ///////////////////////////////////////////////////////////////////////////
    /*!
     * \brief The pw-Sn specific traits.
     */
    class PwSnTraits
    {
    public:
        enum {
            numEq = 2 //!< Number of primary variables
        };
        enum {
            pWIdx = 0,  //!< Idx for the wetting phase pressure in a field vector
            snIdx = 1   //!< Idx for the non-wetting phase saturation in a field vector
        };
    };

    ///////////////////////////////////////////////////////////////////////////
    // PwSnBoxJacobian (evaluate the local jacobian for the newton method.)
    ///////////////////////////////////////////////////////////////////////////
    /*!
     * \brief pw-Sn formulation specific details needed to
     *        approximately calculate the local jacobian in the BOX
     *        scheme.
     *
     * This class is used to fill the gaps in BoxJacobian for the Pw-Sn twophase flow.
     */
    template<class ProblemT, class BoxTraitsT, class PwSnTraitsT>
    class PwSnBoxJacobian : public BoxJacobian<ProblemT, BoxTraitsT, PwSnBoxJacobian<ProblemT, BoxTraitsT, PwSnTraitsT> >
    {
    private:
        typedef PwSnBoxJacobian<ProblemT, BoxTraitsT, PwSnTraitsT>  ThisType;
        typedef BoxJacobian<ProblemT, BoxTraitsT, ThisType >        ParentType;

        typedef ProblemT                                Problem;
        typedef typename Problem::DomainTraits          DomTraits;
        typedef BoxTraitsT                              BoxTraits;
        typedef PwSnTraitsT                             PwSnTraits;

        enum {
            dim            = DomTraits::dim,
            dimWorld       = DomTraits::dimWorld,
            numEq          = BoxTraits::numEq,
            pWIdx          = PwSnTraits::pWIdx,
            snIdx          = PwSnTraits::snIdx
        };

        typedef typename DomTraits::Scalar              Scalar;
        typedef typename DomTraits::CoordScalar         CoordScalar;
        typedef typename DomTraits::Grid                Grid;
        typedef typename DomTraits::Element                Element;
        typedef typename Element::EntityPointer            ElementPointer;
        typedef typename DomTraits::LocalPosition          LocalPosition;

        typedef typename BoxTraits::SolutionVector      SolutionVector;
        typedef typename BoxTraits::FVElementGeometry   FVElementGeometry;
        typedef typename BoxTraits::LocalFunction       LocalFunction;

        typedef FieldMatrix<Scalar, dim, dim>  Tensor;

        /*!
         * \brief Data which is attached to each vert of the and can
         *        be shared between multiple calculations and should
         *        thus be cached in order to increase efficency.
         */
        struct VariableVertexData
        {
            Scalar Sw;
            Scalar pC;
            Scalar pN;

            SolutionVector mobility;  //FieldVector with the number of phases
        };

        /*!
         * \brief Cached data for the each vert of the element.
         */
        struct ElementData
        {
            VariableVertexData  vertex[BoxTraits::ShapeFunctionSetContainer::maxsize];
        };

    public:
        PwSnBoxJacobian(ProblemT &problem)
            : ParentType(problem)
            {};

        /*!
         * \brief Set the current grid element.
         */
        void setCurrentElement(const Element &element)
            {
                if (ParentType::setCurrentElement_(element)) {
                    curElementPorosity_ = ParentType::problem_.porosity(ParentType::curElement_());
                }
            };

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
                partialElementDataUpdate_(curElemDat_,
                                        ParentType::problem_.elementIdx(ParentType::curElement_()),
                                        *curSol_,
                                        vert,
                                        ParentType::problem_.vertIdx(ParentType::curElement_(),
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
         *        volume of a finite volume element in the pw-Sn
         *        formulation.
         *
         * This function should not include the source and sink terms.
         */
        void computeStorage(SolutionVector &result, int scvId, bool usePrevSol) const
            {
                LocalFunction *sol = usePrevSol?prevSol_:curSol_;

                // partial time derivative of the wetting phase mass
                result[pWIdx] = -ParentType::problem_.densityW()
                                   * curElementPorosity_
                                   * (*sol)[scvId][snIdx];
                // partial time derivative of the non-wetting phase mass
                result[snIdx] = ParentType::problem_.densityN()
                                  * curElementPorosity_
                                  * (*sol)[scvId][snIdx];

            }


        /*!
         * \brief Evaluates the mass flux over a face of a subcontrol
         *        volume.
         */
        void computeFlux(SolutionVector &flux, int faceId) const
            {
                Api::require<Api::BasicDomainTraits, typename ProblemT::DomainTraits>();
                assert(numEq == 2);

                const typename FVElementGeometry::SubControlVolumeFace
                    &face = ParentType::curElementGeom_.subContVolFace[faceId];
                const int i = face.i;
                const int j = face.j;

                LocalPosition Kij(0);

                // Kij = K*normal
                ParentType::problem_.applyPermeabilityTensor(Kij,
                                                             ParentType::curElement_(),
                                                             face.normal);

                for (int phase = 0; phase < numEq; phase++) {
                    // calculate FE gradient
                    LocalPosition pGrad(0);
                    for (int k = 0; k < ParentType::curElementGeom_.numVertices; k++) {
                        LocalPosition grad(face.grad[k]);
                        if (phase == snIdx)
                            grad *= curElemDat_.vertex[k].pN;
                        else
                            grad *= (*curSol_)[k][pWIdx];

                        pGrad += grad;
                    }

                    // adjust pressure gradient by gravity force
                    Scalar phaseDensity = ParentType::problem_.density(phase);
                    LocalPosition gravity = ParentType::problem_.gravity();
                    gravity *= phaseDensity;
                    pGrad   -= gravity;

                    // calculate the flux using upwind
                    Scalar outward = pGrad*Kij;
                    if (outward < 0)
                        flux[phase] = phaseDensity*curElemDat_.vertex[i].mobility[phase]*outward;
                    else
                        flux[phase] = phaseDensity*curElemDat_.vertex[j].mobility[phase]*outward;
                }
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



    private:
        /*!
         * \brief Pre-compute the element cache data.
         *
         * This method is called by BoxJacobian (which in turn is
         * called by the operator assembler) every time the current
         * element changes.
         */
        void updateElementData_(ElementData &dest, const LocalFunction &sol)
            {
                assert(numEq == 2);

                int elementIdx   = ParentType::problem_.elementIdx(ParentType::curElement_());
                for (int i = 0; i < ParentType::curElementGeom_.numVertices; i++) {
                    int iGlobal = ParentType::problem_.vertIdx(ParentType::curElement_(), i);
                    partialElementDataUpdate_(dest,
                                            elementIdx,
                                            sol,
                                            i,        // index of sub volume to update,
                                            iGlobal); // global vert index of the sub volume's grid vert
                }
            }


        void partialElementDataUpdate_(ElementData           &dest,
                                     int                  elementIdx,
                                     const LocalFunction &sol,
                                     int                  i, // index of the subvolume/grid vert
                                     int                  iGlobal) // global index of the sub-volume's vert
            {
                // Current cache at sub-controlvolume
                VariableVertexData &scvCache = dest.vertex[i];
                // Current solution for sub-controlvolume
                const SolutionVector &scvSol = sol[i];

                scvCache.Sw = 1.0 - scvSol[snIdx];
                scvCache.pC = ParentType::problem_.pC(ParentType::curElement_(),
                                                      elementIdx,
                                                      i,
                                                      iGlobal,
                                                      scvCache.Sw);
                scvCache.pN = scvSol[pWIdx] + scvCache.pC;
                scvCache.mobility[pWIdx] = ParentType::problem_.mobilityW(ParentType::curElement_(),
                                                                            elementIdx,
                                                                            i,
                                                                            iGlobal,
                                                                            scvCache.Sw);
                scvCache.mobility[snIdx] = ParentType::problem_.mobilityN(ParentType::curElement_(),
                                                                            elementIdx,
                                                                            i,
                                                                            iGlobal,
                                                                            scvSol[snIdx]);
            }

//        ElementPointer     curElement_;
        Scalar          curElementPorosity_;
        Tensor         *curElementPermeability_;

        LocalFunction   *curSol_;
        ElementData        curElemDat_;

        bool             curSolDeflected_;
        Scalar           curSolOrigValue_;
        VariableVertexData curSolOrigVarData_;

        LocalFunction  *prevSol_;
        ElementData       prevElemDat_;
    };


    ///////////////////////////////////////////////////////////////////////////
    // PwSnBoxModel (The actual numerical model.)
    ///////////////////////////////////////////////////////////////////////////
    /*!
     * \brief Adaption of the BOX scheme to the pW-Sn twophase flow model.
     */
    template<class ProblemT>
    class PwSnBoxModel : public BoxScheme< // The implementation of the model
                                           PwSnBoxModel<ProblemT>,

                                           // The Traits for the BOX method
                                           P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                                                       typename ProblemT::DomainTraits::Grid,
                                                       PwSnTraits::numEq>,

                                           // The actual problem we would like to solve
                                           ProblemT,

                                           // The local jacobian operator
                                           PwSnBoxJacobian<ProblemT,
                                                           P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                                                                       typename ProblemT::DomainTraits::Grid,
                                                                       PwSnTraits::numEq>,
                                                           PwSnTraits> >
    {
        typedef typename ProblemT::DomainTraits::Grid   Grid;
        typedef typename ProblemT::DomainTraits::Scalar Scalar;
        typedef PwSnBoxModel<ProblemT>                  ThisType;

    public:
        typedef P1BoxTraits<Scalar, Grid, PwSnTraits::numEq> BoxTraits;
        typedef Dune::PwSnTraits                                   PwSnTraits;

    private:
        typedef PwSnBoxJacobian<ProblemT, BoxTraits, PwSnTraits>  PwSnLocalJacobian;
        typedef BoxScheme<ThisType,
                          BoxTraits,
                          ProblemT,
                          PwSnLocalJacobian>  ParentType;

    public:
        typedef NewNewtonMethod<ThisType> NewtonMethod;

        PwSnBoxModel(ProblemT &prob)
            : ParentType(prob, pwSnLocalJacobian_),
              pwSnLocalJacobian_(prob)
            {
                Api::require<Api::BasicDomainTraits, typename ProblemT::DomainTraits>();
            }

    private:
        // calculates the jacobian matrix at a given position
        PwSnLocalJacobian  pwSnLocalJacobian_;
    };
}

#endif
