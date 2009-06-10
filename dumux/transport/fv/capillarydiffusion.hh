// $Id$
/*****************************************************************************
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
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
#ifndef DUNE_CAPILLARYDIFFUSION_HH
#define DUNE_CAPILLARYDIFFUSION_HH

#include "dumux/transport/fv/diffusivepart.hh"
#include "dumux/transport/transportproblem.hh"

/**
 * @file
 * @brief  Class for defining the diffusive capillary pressure term of a saturation equation
 * @author Bernd Flemisch, Markus Wolff
 */
namespace Dune
{
/*!\ingroup diffPart
 * @brief  Class for defining the diffusive capillary pressure term of a saturation equation
 *
 * Defines the diffusive capillary pressure term of the form
 * \f[
 *  \bar \lambda \boldsymbol{K} \text{grad} \, p_c,
 * \f]
 * where \f$\bar \lambda = \lambda_w f_n = \lambda_n f_w\f$ and \f$\lambda\f$ is a phase mobility and \f$f\f$ a phase fractional flow function,
 * \f$\boldsymbol{K}\f$ is the intrinsic permeability and \f$p_c = p_c(S_w) \f$ the capillary pressure.
 *
 * Template parameters are:

 - GridView      a DUNE gridview type
 - Scalar        type used for scalar quantities
 - VC            type of a class containing different variables of the model
 - Problem       class defining the physical problem
 */
template<class GridView, class Scalar, class VC,
        class Problem = TransportProblem<GridView, Scalar, VC> >
class CapillaryDiffusion: public DiffusivePart<GridView, Scalar>
{
    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
typedef    typename GridView::Grid Grid;
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
    virtual FieldVector operator() (const Element& element, const int indexInInside, Scalar satI, Scalar satJ, const FieldVector& pcGradient) const
    {
        // cell geometry type
        GeometryType gt = element.geometry().type();

        // cell center in reference element
        const LocalPosition& localPos = ReferenceElements<Scalar,dim>::general(gt).position(0,0);

        // get global coordinate of cell center
        const GlobalPosition& globalPos = element.geometry().global(localPos);

        // get absolute permeability of cell
        FieldMatrix permeability(soil_.K(globalPos,element,localPos));

        IntersectionIterator isItEnd = element.ilevelend();
        IntersectionIterator isIt = element.ilevelbegin();
        for (; isIt != isItEnd; ++isIt)
        {
            if(isIt->indexInInside() == indexInInside)
            break;
        }
        int globalIdxI = problem_.variables().indexTransport(element);

        // get geometry type of face
        GeometryType faceGT = isIt->geometryInInside().type();

        //get lambda_bar = lambda_n*f_w
        Scalar mobBarI, mobBarJ;

	if (preComput_)
            mobBarI=problem_.variables().mobilityWetting()[globalIdxI]*problem_.variables().fracFlowFuncNonWetting()[globalIdxI];
	else
	    mobBarI=problem_.materialLaw().mobN(1-satI,globalPos,element,localPos)*problem_.materialLaw().fractionalW(satI,globalPos,element,localPos);

        if (isIt->neighbor())
        {
            // access neighbor
            ElementPointer neighborPointer = isIt->outside();

            int globalIdxJ = problem_.variables().indexTransport(*neighborPointer);

            // compute factor in neighbor
            GeometryType neighborGT = neighborPointer->geometry().type();
            const LocalPosition& localPosNeighbor = ReferenceElements<Scalar,dim>::general(neighborGT).position(0,0);

            // neighbor cell center in global coordinates
            const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().global(localPosNeighbor);

            // take arithmetic average of absolute permeability
            permeability += soil_.K(globalPosNeighbor, *neighborPointer, localPosNeighbor);
            permeability *= 0.5;

            //get lambda_bar = lambda_n*f_w
	    if(preComput_)
        	mobBarJ=problem_.variables().mobilityWetting()[globalIdxJ]*problem_.variables().fracFlowFuncNonWetting()[globalIdxJ];
	    else
        	mobBarJ=problem_.materialLaw().mobN(1-satJ, globalPosNeighbor, *neighborPointer, localPosNeighbor)*problem_.materialLaw().fractionalW(satJ, globalPosNeighbor, *neighborPointer, localPosNeighbor);
        }
        else
        {
            //calculate lambda_n*f_w at the boundary: use a regularization with regularizationparamter eps_
            std::vector<Scalar> mobI = problem_.materialLaw().mob((eps_*satJ + (1-eps_)*satI), globalPos, element, localPos);
            std::vector<Scalar> mobJ = problem_.materialLaw().mob((eps_*satI + (1-eps_)*satJ), globalPos, element, localPos);
            mobBarI = mobI[0]*mobI[1]/(mobI[0]+mobI[1]);
            mobBarJ = mobJ[0]*mobJ[1]/(mobJ[0]+mobJ[1]);
        }

        // set result to K*grad(pc)
        FieldVector result(0);
        permeability.umv(pcGradient, result);

        // set result to f_w*lambda_n*K*grad(pc)
        result *= (mobBarI+mobBarJ)*0.5;

        return result;
    }

    /*! @brief Constructs a CapillaryDiffusion object
     *  @param problem an object of class Dune::TransportProblem or derived
     *  @param soil implementation of the solid matrix
     *  @param preComput if preCompute = true previous calculated mobilities are taken, if preCompute = false new mobilities will be computed (for implicit Scheme)
     */
    CapillaryDiffusion (Problem& problem, Matrix2p<Grid, Scalar>& soil, const bool preComput = true)
    : problem_(problem), soil_(soil), preComput_(preComput), eps_(0.5)
    {}

private:
    Problem& problem_;//problem data
    Matrix2p<Grid, Scalar>& soil_;//object derived from Dune::Matrix2p
    const bool preComput_;//if preCompute = true the mobilities are taken from the variable object, if preCompute = false new mobilities will be taken (for implicit Scheme)
    Scalar eps_;//weighting factor for averaging of the mobilities -> default is central weight
};
}

#endif
