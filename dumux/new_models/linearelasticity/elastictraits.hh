//$Id$

/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch, Melanie Darcis    *
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

#ifndef DUMUX_ELASTICTRAITS_HH
#define DUMUX_ELASTICTRAITS_HH

namespace Dune
{

///////////////////////////////////////////////////////////////////////////
// elastic deformation traits (central place for names and
// indices required by the ElasticBoxJacobian and ElasticBoxModel)
///////////////////////////////////////////////////////////////////////////
/*!
 * \brief The elastic deformation specific traits.
 */
template<class Scalar>
class ElasticTraits
{
    public:
        enum
        {
            numEq = 3,
        //!< Number of primary variables / equations
        };
        enum
        { // Primary variable indices
            uxIdx = 0, //!< Idx for local displacement in x direction
            uyIdx = 1, //!< Idx for local displacement in y direction
            uzIdx = 2
        //!< Idx for local displacement in z direction
        };

        typedef FieldVector<Scalar, 3> DisplacementVector;
        typedef FieldMatrix<Scalar, 3, 3> DisplacementGradient;
        typedef FieldMatrix<Scalar, 3, 3> StressTensor;
        typedef FieldMatrix<Scalar, 3, 3> StrainTensor;
        typedef FieldMatrix<Scalar, 3, 3> IdentityMatrix;

        /*!
         * \brief Data which is attached to each vert of the and can
         *        be shared between multiple calculations and should
         *        thus be cached in order to increase efficency.
         */
        struct VariableVertexData
        {
                Scalar ux;
                Scalar uy;
                Scalar uz;
                Scalar lambda;
                Scalar mu;
        };

};

}

#endif
