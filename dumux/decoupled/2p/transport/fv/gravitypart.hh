// $Id$
/*****************************************************************************
 *   Copyright (C) 2008-2010 by Markus Wolff                                 *
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
#ifndef DUMUX_GRAVITYPART_HH
#define DUMUX_GRAVITYPART_HH

#include "dumux/decoupled/2p/transport/fv/convectivepart.hh"


/**
 * @file
 * @brief  Class for defining the gravity term of a saturation equation
 * @author Markus Wolff
 */

namespace Dumux
{
/*!\ingroup Saturation2p
 * @brief  Class for defining the gravity term of a saturation equation
 *
 * Defines the gravity term of the form
 *
 * \f[\bar \lambda \boldsymbol{K} \, (\rho_n - \rho_w) \, g \, \text{grad} \, z,\f]
 *
 * where \f$\bar \lambda = \lambda_w f_n = \lambda_n f_w\f$ and \f$\lambda\f$ is a phase mobility and \f$f\f$ a phase fractional flow function,
 * \f$ \boldsymbol{K} \f$ is the intrinsic permeability, \f$\rho\f$ is a phase density and  \f$g\f$ is the gravity constant.
 *
 * @tparam TypeTag The Type Tag
 */
template<class TypeTag>
class GravityPart: public ConvectivePart<TypeTag>
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
      typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
      typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
      typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

      typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;
      typedef typename SpatialParameters::MaterialLaw MaterialLaw;

      typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
      typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
    };

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef Dune::FieldVector<Scalar, dim> FieldVector;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

public:
    //! Returns the gravity term
    /*! Returns convective term for current element face
     *  @param[in] element        entity of codim 0
     *  @param[in] indexInInside  face index in reference element
     *  @param[in] satI           saturation of current element
     *  @param[in] satJ           saturation of neighbor element
     *  \return     gravity term of a saturation equation
     */
    FieldVector operator() (const Element& element, const int indexInInside, const Scalar satI, const Scalar satJ) const
    {
        // cell geometry type
        Dune::GeometryType gt = element.geometry().type();

        // cell center in reference element
        typedef Dune::GenericReferenceElements<Scalar, dim> ReferenceElements;
        const LocalPosition& localPos = ReferenceElements::general(gt).position(0,0);

        // get global coordinate of cell center
        const GlobalPosition& globalPos = element.geometry().global(localPos);

        // get absolute permeability of cell
        FieldMatrix permeability(problem_.spatialParameters().intrinsicPermeability(globalPos,element));

        Scalar temperature = problem_.temperature(globalPos, element);
        Scalar referencePressure = problem_.referencePressure(globalPos, element);
        FluidState fluidState;
        fluidState.update(satI, referencePressure, referencePressure, temperature);//not for compressible flow -> thus constant

        IntersectionIterator isItEnd = problem_.gridView().iend(element);
        IntersectionIterator isIt = problem_.gridView().ibegin(element);
        for (; isIt != isItEnd; ++isIt)
        {
            if(isIt->indexInInside() == indexInInside)
            break;
        }
        int globalIdxI = problem_.variables().index(element);

        // get geometry type of face
        Dune::GeometryType faceGT = isIt->geometryInInside().type();

        Scalar potentialW = problem_.variables().potentialWetting(globalIdxI, indexInInside);
        Scalar potentialNW = problem_.variables().potentialNonwetting(globalIdxI, indexInInside);

        //get lambda_bar = lambda_n*f_w
        Scalar lambdaWI = 0;
        Scalar lambdaNWI = 0;
        Scalar lambdaWJ = 0;
        Scalar lambdaNWJ = 0;
        Scalar densityWI = 0;
        Scalar densityNWI = 0;
        Scalar densityWJ = 0;
        Scalar densityNWJ = 0;

        if (preComput_)
        {
            lambdaWI=problem_.variables().mobilityWetting(globalIdxI);
            lambdaNWI=problem_.variables().mobilityNonwetting(globalIdxI);
            densityWI = problem_.variables().densityWetting(globalIdxI);
            densityNWI = problem_.variables().densityNonwetting(globalIdxI);
        }
        else
        {
            lambdaWI = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(globalPos, element), satI);
            lambdaWI /= FluidSystem::phaseViscosity(wPhaseIdx, temperature, referencePressure, fluidState);
            lambdaNWI = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(globalPos, element), satI);
            lambdaNWI /= FluidSystem::phaseViscosity(nPhaseIdx, temperature, referencePressure, fluidState);
            densityWI = FluidSystem::phaseDensity(wPhaseIdx, temperature, referencePressure, fluidState);
            densityNWI = FluidSystem::phaseDensity(nPhaseIdx, temperature, referencePressure, fluidState);
        }

        if (isIt->neighbor())
        {
            // access neighbor
            ElementPointer neighborPointer = isIt->outside();

            int globalIdxJ = problem_.variables().index(*neighborPointer);

            // compute factor in neighbor
            Dune::GeometryType neighborGT = neighborPointer->geometry().type();
            typedef Dune::GenericReferenceElements<Scalar, dim> ReferenceElements;
            const LocalPosition& localPosNeighbor = ReferenceElements::general(neighborGT).position(0,0);

            // neighbor cell center in global coordinates
            const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().global(localPosNeighbor);

            // take arithmetic average of absolute permeability
            permeability += problem_.spatialParameters().intrinsicPermeability(globalPosNeighbor, *neighborPointer);
            permeability *= 0.5;

            //get lambda_bar = lambda_n*f_w
            if (preComput_)
            {
                lambdaWJ=problem_.variables().mobilityWetting(globalIdxJ);
                lambdaNWJ=problem_.variables().mobilityNonwetting(globalIdxJ);
                densityWJ = problem_.variables().densityWetting(globalIdxJ);
                densityNWJ = problem_.variables().densityNonwetting(globalIdxJ);
            }
            else
            {
                lambdaWJ = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(globalPosNeighbor, *neighborPointer), satJ);
                lambdaWJ /= FluidSystem::phaseViscosity(wPhaseIdx, temperature, referencePressure, fluidState);
                lambdaNWJ = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(globalPosNeighbor, *neighborPointer), satJ);
                lambdaNWJ /= FluidSystem::phaseViscosity(nPhaseIdx, temperature, referencePressure, fluidState);
                densityWJ = FluidSystem::phaseDensity(wPhaseIdx, temperature, referencePressure, fluidState);
                densityNWJ = FluidSystem::phaseDensity(nPhaseIdx, temperature, referencePressure, fluidState);
            }
        }
        else
        {
            //calculate lambda_n*f_w at the boundary
            lambdaWJ = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(globalPos, element), satJ);
            lambdaWJ /= FluidSystem::phaseViscosity(wPhaseIdx, temperature, referencePressure, fluidState);
            lambdaNWJ = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(globalPos, element), satJ);
            lambdaNWJ /= FluidSystem::phaseViscosity(nPhaseIdx, temperature, referencePressure, fluidState);
            densityWJ = FluidSystem::phaseDensity(wPhaseIdx, temperature, referencePressure, fluidState);
            densityNWJ = FluidSystem::phaseDensity(nPhaseIdx, temperature, referencePressure, fluidState);
        }

        // set result to K*grad(pc)
        FieldVector result(0);
        permeability.umv(problem_.gravity(), result);

        Scalar lambdaW = (potentialW >= 0) ? lambdaWI : lambdaWJ;
        lambdaW = (potentialW == 0) ? 0.5 * (lambdaWI + lambdaWJ) : lambdaW;
        Scalar lambdaNW = (potentialNW >= 0) ? lambdaNWI : lambdaNWJ;
        lambdaNW = (potentialNW == 0) ? 0.5 * (lambdaNWI + lambdaNWJ) : lambdaNW;
        Scalar densityW = (potentialW >= 0) ? densityWI : densityWJ;
        densityW = (potentialW == 0.) ? 0.5 * (densityWI + densityWJ) : densityW;
        Scalar densityNW = (potentialNW >= 0) ? densityNWI : densityNWJ;
        densityNW = (potentialNW == 0.) ? 0.5 * (densityNWI + densityNWJ) : densityNW;

        // set result to f_w*lambda_n*K*grad(pc)
        result *= lambdaW*lambdaNW/(lambdaW+lambdaNW);
        result *= (densityW - densityNW);

        return result;
    }
    /*! @brief Constructs a GravityPart object
     *  @param problem an object of class Dumux::TransportProblem or derived
     *  @param preComput if preCompute = true previous calculated mobilities are taken, if preCompute = false new mobilities will be computed (for implicit Scheme)
     */
    GravityPart (Problem& problem, const bool preComput = true)
    : ConvectivePart<TypeTag>(problem), problem_(problem), preComput_(preComput)
    {}

private:
    Problem& problem_;//problem data
    const bool preComput_;//if preCompute = true the mobilities are taken from the variable object, if preCompute = false new mobilities will be taken (for implicit Scheme)
};
}

#endif
