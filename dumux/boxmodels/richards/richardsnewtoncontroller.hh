/*****************************************************************************
 *   Copyright (C) 2009-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
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
     * \param uCurrentIter The solution after the current Newton iteration \f$ u^{k+1} \f$
     * \param uLastIter The solution after the last Newton iteration \f$ u^k \f$
     * \param deltaU The vector of differences between the last
     *               iterative solution and the next one \f$ \Delta u^k \f$
     */
    void newtonUpdate(SolutionVector &uCurrentIter,
                      const SolutionVector &uLastIter,
                      const SolutionVector &deltaU)
    {
        this->writeConvergence_(uLastIter, deltaU);
        this->newtonUpdateRelError(uLastIter, deltaU);

        // compute the vertex and element colors for partial
        // reassembly
        if (enablePartialReassemble) {
            Scalar reassembleTol = Dumux::geometricMean(this->error_, 0.1*this->tolerance_);
            reassembleTol = std::max(reassembleTol, 0.1*this->tolerance_);
            this->model_().jacobianAssembler().updateDiscrepancy(uLastIter, deltaU);
            this->model_().jacobianAssembler().computeColors(reassembleTol);
        }

        if (GET_PROP_VALUE(TypeTag, PTAG(NewtonUseLineSearch)))
            lineSearchUpdate_(uCurrentIter, uLastIter, deltaU);
        else {
            // update the solution vector
            uCurrentIter = uLastIter;
            uCurrentIter -= deltaU;

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
                    Scalar pW = uLastIter[globI][pwIdx];
                    Scalar pN = std::max(this->problem_().referencePressure(*eIt, fvElemGeom, i),
                                         pW + pcMin);
                    Scalar pcOld = pN - pW;
                    Scalar SwOld = std::max(0.0, MaterialLaw::Sw(mp, pcOld));

                    // convert into minimum and maximum wetting phase
                    // pressures
                    Scalar pwMin = pN - MaterialLaw::pC(mp, SwOld - 0.2);
                    Scalar pwMax = pN - MaterialLaw::pC(mp, SwOld + 0.2);

                    // clamp the result
                    pW = uCurrentIter[globI][pwIdx];
                    pW = std::max(pwMin, std::min(pW, pwMax));
                    uCurrentIter[globI][pwIdx] = pW;

                }
            }
        }
    }

private:
    void lineSearchUpdate_(SolutionVector &uCurrentIter,
                           const SolutionVector &uLastIter,
                           const SolutionVector &deltaU)
    {
       Scalar lambda = 1.0;
       Scalar globDef;
       SolutionVector tmp(uLastIter);
       Scalar oldGlobDef = this->model_().globalResidual(tmp, uLastIter);

       while (true) {
           uCurrentIter = deltaU;
           uCurrentIter *= -lambda;
           uCurrentIter += uLastIter;

           // calculate the residual of the current solution
           globDef = this->model_().globalResidual(tmp, uCurrentIter);

           if (globDef < oldGlobDef || lambda <= 1.0/64) {
               this->endIterMsg() << ", defect " << oldGlobDef << "->"  << globDef << "@lambda=" << lambda;
               return;
           }

           // try with a smaller update
           lambda /= 2;
       }
    };
};
}

#endif
