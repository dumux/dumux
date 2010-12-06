/*****************************************************************************
 *   Copyright (C) 2009-2010 by Andreas Lauser                               *
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
/*!
 * \file
 * \brief A newton solver specific to the Richards problem.
 */
#ifndef DUMUX_RICHARDS_NEWTON_CONTROLLER_HH
#define DUMUX_RICHARDS_NEWTON_CONTROLLER_HH

#include "richardsproperties.hh"

#include <dumux/nonlinear/newtoncontroller.hh>

namespace Dumux {
/*!
 * \brief A Richards model specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * and can thus do update smarter than the plain Newton controller.
 */
template <class TypeTag>
class RichardsNewtonController : public NewtonController<TypeTag>
{
    typedef RichardsNewtonController<TypeTag> ThisType;
    typedef NewtonController<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw)) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLawParams)) MaterialLawParams;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(RichardsIndices)) Indices;
    enum {
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
        dim = GridView::dimension,
        
        pwIdx = Indices::pwIdx,
    };

    enum { enablePartialReassemble = GET_PROP_VALUE(TypeTag, PTAG(EnablePartialReassemble)) };

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef Dune::FieldVector<Scalar, dim> GlobalPosition;

public:
    /*!
     * \brief Constructor
     */
    RichardsNewtonController()
    { };

    /*!
     * \brief Update the current solution of the newton method
     *
     * This is basically the step
     * \f[ u^{k+1} = u^k - \Delta u^k \f]
     *
     * \param deltaU When the method is called, this contains the
     *               vector of differences between the current
     *               iterative solution and the next \f$\Delta
     *               u\f$. After the method is finished it should
     *               contain the next iterative solution \f$ u^{k+1} \f$.
     * \param uOld The current iterative solution \f$ u^k \f$
     */
    void newtonUpdate(SolutionVector &deltaU, const SolutionVector &uOld)
    {
        this->writeConvergence_(uOld, deltaU);
        this->newtonUpdateRelError(uOld, deltaU);

        // compute the vertex and element colors for partial
        // reassembly
        if (enablePartialReassemble) {
            Scalar reassembleTol = 0.3*Dumux::geometricMean(this->error_, 
                                                            this->tolerance_);
            reassembleTol = std::max(reassembleTol, this->tolerance_);
            this->model_().jacobianAssembler().updateDiscrepancy(uOld, deltaU);
            this->model_().jacobianAssembler().computeColors(reassembleTol);
        }

        if (GET_PROP_VALUE(TypeTag, PTAG(NewtonUseLineSearch)))
            lineSearchUpdate_(deltaU, uOld);
        else {
            // update the solution vector
            deltaU *= -1;
            deltaU += uOld;
            
            // clamp saturation change to at most 20% per iteration
            FVElementGeometry fvElemGeom;
            const GridView &gv = this->problem_().gridView();
            ElementIterator eIt = gv.template begin<0>();
            const ElementIterator eEndIt = gv.template end<0>();
            for (; eIt != eEndIt; ++eIt) {
                fvElemGeom.update(gv, *eIt);
                for (int i = 0; i < fvElemGeom.numVertices; ++i) {
                    int globI = this->problem_().vertexMapper().map(*eIt, i, dim);

                    // calculate the old wetting phase saturation
                    const SpatialParameters &sp = this->problem_().spatialParameters();
                    const MaterialLawParams &mp = sp.materialLawParams(*eIt, fvElemGeom, i);
                    Scalar pcMin = MaterialLaw::pC(mp, 1.0);
                    Scalar pW = uOld[globI][pwIdx];
                    Scalar pN = std::max(this->problem_().referencePressure(*eIt, fvElemGeom, i), 
                                         pW + pcMin);
                    Scalar pcOld = pN - pW;
                    Scalar SwOld = std::max(0.0, MaterialLaw::Sw(mp, pcOld));

                    // convert into minimum and maximum wetting phase
                    // pressures
                    Scalar pwMin = pN - MaterialLaw::pC(mp, SwOld - 0.2);
                    Scalar pwMax = pN - MaterialLaw::pC(mp, SwOld + 0.2);
                    
                    // clamp the result
                    deltaU[globI][pwIdx] = std::max(pwMin, std::min(deltaU[globI][pwIdx], pwMax));
                }
            }
        }
    }

private:
    void lineSearchUpdate_(SolutionVector &u, const SolutionVector &uOld)
    {
       Scalar lambda = 1.0;
       Scalar globDef;
       SolutionVector tmp(this->model_(), 0.0);
       Scalar oldGlobDef = this->model_().globalResidual(tmp, uOld);

       int n = 0;
       while (true) {
           u *= -lambda;
           u += uOld;
           globDef = this->model_().globalResidual(tmp);

           if (globDef < oldGlobDef || lambda <= 1.0/64) {
               this->endIterMsg() << ", defect " << oldGlobDef << "->"  << globDef << "@lambda=2^-" << n;
               return;
           }

           // undo the last iteration
           u -= uOld;
           u /= - lambda;

           // try with a smaller update
           lambda /= 2;
           ++n;
       }
    };
};
}

#endif
