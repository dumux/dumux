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
#ifndef DUNE_GRAVITYPART_HH
#define DUNE_GRAVITYPART_HH

#include "dumux/transport/fv/convectivepart.hh"
#include "dumux/transport/transportproblem.hh"

/**
 * @file
 * @brief  Class for defining the gravity term of a saturation equation
 * @author Markus Wolff
 */

namespace Dune
{
/*!\ingroup convPart
 * @brief  Class for defining the gravity term of a saturation equation
 *
 * Defines the gravity term of the form
 * \f[
 *  \bar \lambda \boldsymbol{K} \, (\rho_n - \rho_w) \, g \, \text{grad} \, z,
 * \f]
 *
 * where \f$\bar \lambda = \lambda_w f_n = \lambda_n f_w\f$ and \f$\lambda\f$ is a phase mobility and \f$f\f$ a phase fractional flow function,
 * \f$ \boldsymbol{K} \f$ is the intrinsic permeability, \f$\rho\f$ is a phase density and  \f$g\f$ is the gravity constant.

 * Template parameters are:

 - GridView      a DUNE gridview type
 - Scalar        type used for scalar quantities
 - VC            type of a class containing different variables of the model
 - Problem       class defining the physical problem
 */
template<class GridView, class Scalar, class VC,
        class Problem = TransportProblem<GridView, Scalar, VC> >
class GravityPart: public ConvectivePart<GridView, Scalar>
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
    //! Returns the gravity term
    /*! Returns convective term for current element face
     *  @param[in] element        entity of codim 0
     *  @param[in] satI           saturation of current element
     *  @param[in] satJ           saturation of neighbor element
     *  @param[in] indexInInside  face index in reference element
     *  \return     gravity term of a saturation equation
     */
    virtual FieldVector operator() (const Element& element, const Scalar satI, const Scalar satJ, const int indexInInside) const
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

        Scalar potentialW = problem_.variables().potentialWetting()[globalIdxI][indexInInside];
        Scalar potentialNW = problem_.variables().potentialNonWetting()[globalIdxI][indexInInside];

        //get lambda_bar = lambda_n*f_w
        Scalar lambdaWI = 0;
        Scalar lambdaNWI = 0;
        Scalar lambdaWJ = 0;
        Scalar lambdaNWJ = 0;

        if (preComput_)
        {
            lambdaWI=problem_.variables().mobilityWetting()[globalIdxI];
            lambdaNWI=problem_.variables().mobilityNonWetting()[globalIdxI];
        }
        else
        {
            std::vector<Scalar> mobilities = problem_.materialLaw().mob(satI,globalPos,element,localPos);
            lambdaWI = mobilities[0];
            lambdaNWI = mobilities[1];
        }

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
            if (preComput_)
            {
                lambdaWJ=problem_.variables().mobilityWetting()[globalIdxJ];
                lambdaNWJ=problem_.variables().mobilityNonWetting()[globalIdxJ];
            }
            else
            {
                std::vector<Scalar> mobilities = problem_.materialLaw().mob(satJ,globalPos,element,localPos);
                lambdaWJ = mobilities[0];
                lambdaNWJ = mobilities[1];
            }
        }
        else
        {
            std::vector<Scalar> mobilities = problem_.materialLaw().mob(satJ,globalPos,element,localPos);
            lambdaWJ = mobilities[0];
            lambdaNWJ = mobilities[1];
        }

        // set result to K*grad(pc)
        FieldVector result(0);
        permeability.umv(gravity_, result);

        Scalar lambdaW = (potentialW >= 0) ? lambdaWI : lambdaWJ;
        Scalar lambdaNW = (potentialNW >= 0) ? lambdaNWI : lambdaNWJ;

        // set result to f_w*lambda_n*K*grad(pc)
        result *= lambdaW*lambdaNW/(lambdaW+lambdaNW);

        return result;
    }
    /*! @brief Constructs a GravityPart object
     *  @param problem an object of class Dune::TransportProblem or derived
     *  @param soil implementation of the solid matrix
     *  @param preComput if preCompute = true previous calculated mobilities are taken, if preCompute = false new mobilities will be computed (for implicit Scheme)
     */
    GravityPart (Problem& problem, Matrix2p<Grid, Scalar>& soil, const bool preComput = true)
    : problem_(problem), soil_(soil), preComput_(preComput)
    {
        double rhoDiff = problem.wettingPhase().density() - problem.nonWettingPhase().density();
        gravity_ = problem.gravity();
        gravity_ *= rhoDiff;
    }

private:
    Problem& problem_;//problem data
    Matrix2p<Grid, Scalar>& soil_;//object derived from Dune::Matrix2p
    const bool preComput_;//if preCompute = true the mobilities are taken from the variable object, if preCompute = false new mobilities will be taken (for implicit Scheme)
    FieldVector gravity_;//gravity vector
};
}

#endif
