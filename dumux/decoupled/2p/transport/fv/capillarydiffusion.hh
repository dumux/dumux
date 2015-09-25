// $Id$
/*****************************************************************************
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
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
#ifndef DUMUX_CAPILLARYDIFFUSION_HH
#define DUMUX_CAPILLARYDIFFUSION_HH

#include "dumux/decoupled/2p/transport/fv/diffusivepart.hh"

/**
 * @file
 * @brief  Class for defining the diffusive capillary pressure term of a 2p saturation equation
 * @author Bernd Flemisch, Markus Wolff
 */
namespace Dumux
{
/*!\ingroup Saturation2p
 * @brief  Class for defining the diffusive capillary pressure term of a saturation equation
 *
 * Defines the diffusive capillary pressure term of the form
 *
 * \f[\bar \lambda \boldsymbol{K} \text{grad} \, p_c,\f]
 *
 * where \f$\bar \lambda = \lambda_w f_n = \lambda_n f_w\f$ and \f$\lambda\f$ is a phase mobility and \f$f\f$ a phase fractional flow function,
 * \f$\boldsymbol{K}\f$ is the intrinsic permeability and \f$p_c = p_c(S_w) \f$ the capillary pressure.
 *
 * @tparam TypeTag The Type Tag
 */
template<class TypeTag>
class CapillaryDiffusion: public DiffusivePart<TypeTag>
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
    //! Returns capillary diffusion term
    /*! Returns capillary diffusion term for current element face
     *  @param[in] element        entity of codim 0
     *  @param[in] indexInInside  face index in reference element
     *  @param[in] satI           saturation of current element
     *  @param[in] satJ           saturation of neighbor element
     *  @param[in] pcGradient     gradient of capillary pressure between element I and J
     *  \return     capillary pressure term of the saturation equation
     */
    FieldVector operator() (const Element& element, const int indexInInside, Scalar satI, Scalar satJ, const FieldVector& pcGradient) const
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

        Scalar temperature = problem_.temperature(globalPos, element);
        Scalar referencePressure = problem_.referencePressure(globalPos, element);

        //get lambda_bar = lambda_n*f_w
        Scalar mobBar = 0;
        Scalar mobilityWI = 0;
        Scalar mobilityNWI = 0;

        if (preComput_)
        {
            mobilityWI = problem_.variables().mobilityWetting(globalIdxI);
            mobilityNWI = problem_.variables().mobilityNonwetting(globalIdxI);
        }
        else
        {
            FluidState fluidState;
            fluidState.update(satI, referencePressure, referencePressure, temperature);
            mobilityWI = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(globalPos, element), satI);
            mobilityWI /= FluidSystem::phaseViscosity(wPhaseIdx, temperature, referencePressure, fluidState);
            mobilityNWI = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(globalPos, element), satI);
            mobilityNWI /= FluidSystem::phaseViscosity(nPhaseIdx, temperature, referencePressure, fluidState);
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

            // distance vector between barycenters
            FieldVector distVec = globalPosNeighbor - globalPos;

            // compute distance between cell centers
            Scalar dist = distVec.two_norm();

            FieldVector unitDistVec(distVec);
            unitDistVec /= dist;

            // get absolute permeability
            FieldMatrix permeabilityJ(problem_.spatialParameters().intrinsicPermeability(globalPosNeighbor, *neighborPointer));

            // harmonic mean of permeability
            for (int x = 0;x<dim;x++)
            {
                for (int y = 0; y < dim;y++)
                {
                    if (permeability[x][y] && permeabilityJ[x][y])
                    {
                        permeability[x][y]= 2*permeability[x][y]*permeabilityJ[x][y]/(permeability[x][y]+permeabilityJ[x][y]);
                    }
                }
            }
            Scalar mobilityWJ = 0;
            Scalar mobilityNWJ = 0;
            //get lambda_bar = lambda_n*f_w
            if(preComput_)
            {
                mobilityWJ = problem_.variables().mobilityWetting(globalIdxJ);
                mobilityNWJ = problem_.variables().mobilityNonwetting(globalIdxJ);
            }
            else
            {
                FluidState fluidState;
                fluidState.update(satJ, referencePressure, referencePressure, temperature);
                mobilityWJ = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(globalPosNeighbor, *neighborPointer), satJ);
                mobilityWJ /= FluidSystem::phaseViscosity(wPhaseIdx, temperature, referencePressure, fluidState);
                mobilityNWJ = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(globalPosNeighbor, *neighborPointer), satJ);
                mobilityNWJ /= FluidSystem::phaseViscosity(nPhaseIdx, temperature, referencePressure, fluidState);
            }
            Scalar mobilityWMean = 0.5*(mobilityWI + mobilityWJ);
            Scalar mobilityNWMean = 0.5*(mobilityNWI + mobilityNWJ);
            mobBar = mobilityWMean*mobilityNWMean/(mobilityWMean+mobilityNWMean);
         }//end intersection with neighbor
        else
        {
            Scalar mobilityWJ = 0;
            Scalar mobilityNWJ = 0;

            //calculate lambda_n*f_w at the boundary
            FluidState fluidState;
            fluidState.update(satJ, referencePressure, referencePressure, temperature);
            mobilityWJ = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(globalPos, element), satJ);
            mobilityWJ /= FluidSystem::phaseViscosity(wPhaseIdx, temperature, referencePressure, fluidState);
            mobilityNWJ = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(globalPos, element), satJ);
            mobilityNWJ /= FluidSystem::phaseViscosity(nPhaseIdx, temperature, referencePressure, fluidState);

            Scalar mobWMean = 0.5 * (mobilityWI + mobilityWJ);
            Scalar mobNWMean = 0.5 * (mobilityNWI + mobilityNWJ);

            mobBar = mobWMean * mobNWMean / (mobWMean + mobNWMean);
        }

        // set result to K*grad(pc)
        FieldVector result(0);
        permeability.umv(pcGradient, result);

        // set result to f_w*lambda_n*K*grad(pc)
        result *= mobBar;

        return result;
    }

    /*! @brief Constructs a CapillaryDiffusion object
     *  @param problem an object of class Dumux::TransportProblem or derived
     *  @param preComput if preCompute = true previous calculated mobilities are taken, if preCompute = false new mobilities will be computed (for implicit Scheme)
     */
    CapillaryDiffusion (Problem& problem, const bool preComput = true)
    : DiffusivePart<TypeTag>(problem), problem_(problem), preComput_(preComput)
    {}

private:
    Problem& problem_;//problem data
    const bool preComput_;//if preCompute = true the mobilities are taken from the variable object, if preCompute = false new mobilities will be taken (for implicit Scheme)
};
}

#endif
