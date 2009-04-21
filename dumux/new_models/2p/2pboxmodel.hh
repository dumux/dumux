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

#include <dumux/auxiliary/math.hh>

#include "2ptraits.hh"
#include "2pvertexdata.hh"

#include <dumux/auxiliary/apis.hh>

namespace Dune
{
/*!
 * \brief Calculate the local Jacobian for the two phase model in the
 *        BOX scheme.
 */
template<class ProblemT, class BoxTraitsT, class TwoPTraitsT>
class TwoPBoxJacobian : public BoxJacobian<ProblemT, 
                                           BoxTraitsT,
                                           TwoPBoxJacobian<ProblemT, BoxTraitsT, TwoPTraitsT>,
                                           TwoPVertexData<TwoPTraitsT, ProblemT> >
{
private:
    typedef TwoPBoxJacobian<ProblemT,
                            BoxTraitsT,
                            TwoPTraitsT>            ThisType;
    typedef TwoPVertexData<TwoPTraitsT, ProblemT>   VertexData;
    typedef BoxJacobian<ProblemT, 
                        BoxTraitsT, 
                        ThisType,
                        VertexData>   ParentType;

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

    typedef std::vector<VertexData>        VertexDataArray;
    typedef FieldMatrix<Scalar, dim, dim>  Tensor;

public:
    TwoPBoxJacobian(ProblemT &problem)
        : ParentType(problem)
    {};

    /*!
     * \brief Evaluate the amount all conservation quantites
     *        (e.g. phase mass) within a finite volume.
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
   
        // wetting phase mass
        result[wMassIdx] =
            vertDat.density[wPhase]
            * vertDat.porosity
            * vertDat.satW;

        // non-wetting phase mass
        result[nMassIdx] = 
            vertDat.density[nPhase]
            * vertDat.porosity
            * vertDat.satN;
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

        GlobalPosition pressureGrad[numPhases];
        Scalar densityAtIP[numPhases];
        for (int i = 0; i < numPhases; ++i) {
            pressureGrad[i] = 0;
            densityAtIP[i] = 0;
        }

        // calculate pressure gradients
        GlobalPosition tmp(0.0);
        for (int idx = 0; 
             idx < this->curElementGeom_.numVertices;
             idx++) // loop over adjacent vertices
        {
            // FE gradient at vertex idx
            const LocalPosition &feGrad = face.grad[idx];

            // compute sum of pressure gradients for each phase
            for (int phase = 0; phase < numPhases; phase++)
            {
                // the pressure gradient
                tmp = feGrad;
                tmp *= this->curElemDat_[idx].pressure[phase];
                pressureGrad[phase] += tmp;

                // phase density
                densityAtIP[phase] 
                    += 
                    this->curElemDat_[idx].density[phase] *
                    face.shapeValue[idx];
            }
        }
        
        // correct the pressure gradients by the hydrostatic
        // pressure due to gravity
        for (int phase=0; phase < numPhases; phase++)
        {
            tmp = this->problem_.gravity();
            tmp *= densityAtIP[phase];
            
            pressureGrad[phase] -= tmp;
        }


        // calculate the permeability tensor
        const Tensor &Ki = this->problem_.soil().K(global_i, ParentType::curElement_(), local_i);
        const Tensor &Kj = this->problem_.soil().K(global_j, ParentType::curElement_(), local_j);
        Tensor K;
        harmonicMeanMatrix(K, Ki, Kj);
        
        for (int phase = 0; phase < numPhases; phase++) 
        {
            GlobalPosition vDarcy;
            K.mv(pressureGrad[phase], vDarcy);  // vDarcy = K * grad p

            // calculate the flux using upwind
            Scalar vDarcyOut = vDarcy*face.normal;
            int upwindIdx;
            if (vDarcyOut < 0)
                upwindIdx = i;
            else 
                upwindIdx = j;
            
            Scalar fluxVal =
                1000.0* //this->curElemDat_[upwindIdx].density[phase]*
                this->curElemDat_[upwindIdx].mobility[phase]*
                vDarcyOut;
//            std::cout << "density: " << this->curElemDat_[upwindIdx].density[phase] << "\n";
//            std::cout << "mobility: " << this->curElemDat_[upwindIdx].mobility[phase] << "\n";
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
     * \brief Return the temperature given the solution vector of a
     *        finite volume.
     */
    template <class SolutionVector>
    Scalar temperature(const SolutionVector &sol)
    {
        return this->problem_.temperature(); /* constant temperature */ 
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
