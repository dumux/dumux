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
#ifndef DUMUX_TWOP_BOX_MODEL_HH
#define DUMUX_TWOP_BOX_MODEL_HH

#include <dumux/new_models/boxscheme/boxscheme.hh>
#include <dumux/new_models/boxscheme/p1boxtraits.hh>

#include "2ptraits.hh"

#include <dumux/auxiliary/apis.hh>

namespace Dune
{
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
template<class ProblemT, class BoxTraitsT, class TwoPTraitsT>
class TwoPBoxJacobian : public BoxJacobian<ProblemT, 
                                           BoxTraitsT,
                                           TwoPBoxJacobian<ProblemT, BoxTraitsT, TwoPTraitsT>,
                                           typename TwoPTraitsT::ElementData,
                                           typename TwoPTraitsT::VertexData >
{
private:
    typedef TwoPBoxJacobian<ProblemT,
                            BoxTraitsT,
                            TwoPTraitsT>                   ThisType;
    typedef BoxJacobian<ProblemT, 
                        BoxTraitsT, 
                        ThisType,
                        typename TwoPTraitsT::ElementData,
                        typename TwoPTraitsT::VertexData>   ParentType;

    typedef ProblemT                       Problem;
    typedef typename Problem::DomainTraits DomTraits;
    typedef BoxTraitsT                     BoxTraits;
    typedef TwoPTraitsT                    TwoPTraits;

    enum {
        dim            = DomTraits::dim,
        dimWorld       = DomTraits::dimWorld,

        numEq          = BoxTraits::numEq,
        numPhases      = TwoPTraits::numPhases,

        pressureIdx    = TwoPTraits::pressureIdx,
        saturationIdx  = TwoPTraits::saturationIdx,

        wMassIdx       = TwoPTraits::wMassIdx,
        nMassIdx       = TwoPTraits::nMassIdx,
        
        formulation    = TwoPTraits::formulation,

        wPhase         = TwoPTraits::wPhase,
        nPhase         = TwoPTraits::nPhase,
    };

    typedef typename DomTraits::Scalar              Scalar;
    typedef typename DomTraits::CoordScalar         CoordScalar;
    typedef typename DomTraits::Grid                Grid;
    typedef typename DomTraits::Element             Element;
    typedef typename DomTraits::ElementIterator     ElementIterator;
    typedef typename Element::EntityPointer         ElementPointer;

    typedef typename DomTraits::LocalPosition       LocalPosition;
    typedef typename DomTraits::GlobalPosition      GlobalPosition;

    typedef typename BoxTraits::SolutionVector      SolutionVector;
    typedef typename BoxTraits::FVElementGeometry   FVElementGeometry;
    typedef typename BoxTraits::SpatialFunction     SpatialFunction;
    typedef typename BoxTraits::LocalFunction       LocalFunction;

    typedef FieldMatrix<Scalar, dim, dim>  Tensor;

    typedef typename TwoPTraits::VertexData   VertexData;
    typedef typename TwoPTraits::ElementData  ElementData;

public:
    TwoPBoxJacobian(ProblemT &problem)
        : ParentType(problem)
    {};

    /*!
     * \brief Evaluate the amount all conservation quantites
     *        (e.g. phase mass) within a finite volume.
     */
    void computeStorage(SolutionVector &result, int scvId, bool usePrevSol) const
    {
        // wetting phase mass
        result[wMassIdx] =
            this->curElemDat_[scvId].density[wPhase]
            * this->curElemDat_[scvId].porosity
            * this->curElemDat_[scvId].satW;

        // non-wetting phase mass
        result[nMassIdx] = 
            this->curElemDat_[scvId].density[nPhase]
            * this->curElemDat_[scvId].porosity
            * this->curElemDat_[scvId].satN;
    }

    /*!
     * \brief Evaluates the mass flux over a face of a subcontrol
     *        volume.
     */
    void computeFlux(SolutionVector &flux, int faceId) const
    {
        assert(numEq == 2);
        
        const typename FVElementGeometry::SubControlVolumeFace
            &face = ParentType::curElementGeom_.subContVolFace[faceId];
        const int i = face.i;
        const int j = face.j;

        // get global coordinates of verts i,j
        const GlobalPosition &global_i = this->curElementGeom_.subContVol[i].global;
        const GlobalPosition &global_j = this->curElementGeom_.subContVol[j].global;

        // get local coordinates of verts i,j
        const LocalPosition &local_i = this->curElementGeom_.subContVol[i].local;
        const LocalPosition &local_j = this->curElementGeom_.subContVol[j].local;

        // calculate the permeability tensor
        Tensor K         = this->problem_.soil().K(global_i, ParentType::curElement_(), local_i);
        const Tensor &Kj = this->problem_.soil().K(global_j, ParentType::curElement_(), local_j);
        harmonicMeanK_(K, Kj);
        
        GlobalPosition Kij;
        K.mv(face.normal, Kij); // Kij = K*normal
        
        for (int phase = 0; phase < numPhases; phase++) {
            // calculate FE gradient
            LocalPosition pGrad(0);
            Scalar densityIJ = 0.0;
            for (int k = 0; k < ParentType::curElementGeom_.numVertices; k++) {
                LocalPosition grad(face.grad[k]);
                grad *= this->curElemDat_[k].pressure[phase];
                pGrad += grad;

                densityIJ += 
                    this->curElemDat_[k].density[phase] *
                    this->curElementGeom_.subContVolFace[faceId].shapeValue[k];
            }

            // adjust pressure gradient by gravity force
            GlobalPosition tmp = ParentType::problem_.gravity();
            tmp   *= densityIJ;
            pGrad -= tmp;

            // calculate the flux using upwind
            Scalar vDarcyOut = pGrad*Kij;
            int upwindIdx;
            if (vDarcyOut < 0)
                upwindIdx = i;
            else 
                upwindIdx = j;
            
            Scalar fluxVal =
                this->curElemDat_[upwindIdx].density[phase]*
                this->curElemDat_[upwindIdx].mobility[phase]*
                vDarcyOut;

            if (phase == wPhase)
                flux[wMassIdx] = fluxVal;
            else if (phase == nPhase)
                flux[nMassIdx] = fluxVal;
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
        ScalarField *pW =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *pN =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *pC =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *Sw =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *Sn =           writer.template createField<Scalar, 1>(numVertices);
        
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

                (*pW)[globalIdx] = this->curElemDat_[i].pressure[wPhase];
                (*pN)[globalIdx] = this->curElemDat_[i].pressure[nPhase];
                (*pC)[globalIdx] = this->curElemDat_[i].pC;
                (*Sw)[globalIdx] = this->curElemDat_[i].satW;
                (*Sn)[globalIdx] = this->curElemDat_[i].satN;
            };
        }

        writer.addVertexData(pW, "pW");
        writer.addVertexData(pN, "pN");
        writer.addVertexData(pC, "pC");
        writer.addVertexData(Sw, "SW");
        writer.addVertexData(Sn, "SN");
    }

    void updateVertexData_(ElementData         &dest,
                           const LocalFunction &sol,
                           int                  vertIdx,
                           bool                 isOldSol) // index of the subvolume/grid vertex
    {
        // Current cache at sub-controlvolume
        VertexData &vertDat = dest[vertIdx];
        // Current solution for sub-controlvolume
        const SolutionVector &vertSol = sol[vertIdx];
        
        // get local coordinates of the vertex
        const GlobalPosition &globalPos = this->curElementGeom_.subContVol[vertIdx].global;
        const LocalPosition &localPos = this->curElementGeom_.subContVol[vertIdx].local;

        if (formulation == TwoPTraits::pWsN) {
            vertDat.satN = vertSol[saturationIdx];
            vertDat.satW = 1.0 - vertDat.satN;
            vertDat.pC = this->problem_.materialLaw().pC(vertDat.satW,
                                                         globalPos,
                                                         this->curElement_(),
                                                         localPos);
            vertDat.pressure[wPhase] = vertSol[pressureIdx];
            vertDat.pressure[nPhase] = vertDat.pressure[wPhase] + vertDat.pC;
        }
        else if (formulation == TwoPTraits::pNsW) {
            vertDat.satW = vertSol[saturationIdx];
            vertDat.satN = 1.0 - vertDat.satW;
            vertDat.pC = this->problem_.materialLaw().pC(vertDat.satW,
                                                         globalPos,
                                                         this->curElement_(),
                                                         localPos);
            vertDat.pressure[nPhase] = vertSol[pressureIdx];
            vertDat.pressure[wPhase] = vertDat.pressure[nPhase] - vertDat.pC;
        }

        vertDat.porosity = this->problem_.porosity(this->curElement_(), 
                                                   vertIdx);

        vertDat.density[wPhase] = this->problem_.wettingPhase().density(283.15, // temperature,
                                                                        vertDat.pressure[wPhase]);
        vertDat.density[nPhase] = this->problem_.nonwettingPhase().density(283.15, //temperature,
                                                                           vertDat.pressure[nPhase]);

        vertDat.mobility[wPhase] = this->problem_.materialLaw().mobW(vertDat.satW,
                                                                     globalPos,
                                                                     this->curElement_(),
                                                                     localPos,
                                                                     283.15, // temperature,
                                                                     vertDat.pressure[wPhase]);
        vertDat.mobility[nPhase] = this->problem_.materialLaw().mobN(vertDat.satN,
                                                                     globalPos,
                                                                     this->curElement_(),
                                                                     localPos,
                                                                     283.15, // temperature,
                                                                     vertDat.pressure[nPhase]);
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
};


///////////////////////////////////////////////////////////////////////////
// TwoPBoxModel (The actual numerical model.)
///////////////////////////////////////////////////////////////////////////
/*!
 * \brief Adaption of the BOX scheme to the pW-Sn twophase flow model.
 */
template<class ProblemT, 
         class TwoPTraitsT = PwSnTwoPTraits<typename ProblemT::DomainTraits::Scalar> >
class TwoPBoxModel : public BoxScheme<TwoPBoxModel<ProblemT>, // Implementation
                                      
                                      // The Traits for the BOX method
                                      P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                                                  typename ProblemT::DomainTraits::Grid,
                                                  TwoPTraitsT::numEq>,
                                      
                                      // The actual problem we would like to solve
                                      ProblemT,
                                      
                                      // The local jacobian operator
                                      TwoPBoxJacobian<ProblemT,
                                                      P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                                                                  typename ProblemT::DomainTraits::Grid,
                                                                  TwoPTraitsT::numEq>,
                                                      TwoPTraitsT> >
{
    typedef typename ProblemT::DomainTraits::Grid   Grid;
    typedef typename ProblemT::DomainTraits::Scalar Scalar;
    typedef TwoPBoxModel<ProblemT>                  ThisType;

public:
    typedef TwoPTraitsT                                  TwoPTraits;
    typedef P1BoxTraits<Scalar, Grid, TwoPTraits::numEq> BoxTraits;

private:
    typedef TwoPBoxJacobian<ProblemT, BoxTraits, TwoPTraits>  TwoPLocalJacobian;
    typedef BoxScheme<ThisType,
                      BoxTraits,
                      ProblemT,
                      TwoPLocalJacobian>  ParentType;

public:
    typedef NewNewtonMethod<ThisType> NewtonMethod;

    TwoPBoxModel(ProblemT &prob)
        : ParentType(prob, twoPLocalJacobian_),
          twoPLocalJacobian_(prob)
    {
        Api::require<Api::BasicDomainTraits, typename ProblemT::DomainTraits>();
    }

    /*!
     * \brief Add the mass fraction of air in water to VTK output of
     *        the current timestep.
     */
    template <class MultiWriter>
    void addVtkFields(MultiWriter &writer)
    {
        twoPLocalJacobian_.addVtkFields(writer, this->currentSolution());
    }


private:
    // calculates the jacobian matrix at a given position
    TwoPLocalJacobian  twoPLocalJacobian_;
};
}

#endif
