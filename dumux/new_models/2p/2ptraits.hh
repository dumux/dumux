/*****************************************************************************
 *   Copyright (C) 2008 by Klaus Mosthaf, Andreas Lauser, Bernd Flemisch                    *
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
#ifndef DUMUX_2PTRAITS_HH
#define DUMUX_2PTRAITS_HH

#include <dune/common/fvector.hh>

namespace Dune
{
///////////////////////////////////////////////////////////////////////////
// The traits for the twophase model
///////////////////////////////////////////////////////////////////////////
/*!
 * \brief The traits for the twophase model.
 */
template <class Scalar, int numPhases_ = 2>
class TwoPBaseTraits
{
public:
    // basic constants
    static const int numEq = 2;              //!< Number of primary variables
    static const int numPhases = numPhases_; //!< Number of fluid phases
    
    // Indices of the primary variables
    static const int pressureIdx = 0;   //!< Index for either wetting or non-wetting phase pressure (depending on formulation) in a solution vector
    static const int saturationIdx = 1; //!< Index for the non-wetting phase saturation in a field vector
    
    // formulations
    static const int pWsN = 0; //!< Pw and Sn as primary variables
    static const int pNsW = 1; //!< Pn and Sw as primary variables

    // phase indices
    static const int wPhase = 0; //!< Phase index of the wetting phase
    static const int nPhase = 1;  //!< Phase index of the non-wetting phase

    typedef FieldVector<Scalar, numPhases>  PhasesVector;
};

/*!
 * \brief The traits for the Pw-Sn formulation of the twophase model.
 */
template <class Scalar>
class PwSnTwoPTraits : public TwoPBaseTraits<Scalar>
{
    typedef TwoPBaseTraits<Scalar> ParentType;

public:
    // set the actual formulation to pW-Sn
    static const int formulation = ParentType::pWsN;

    // the index of the wetting phase mass in solution vectors
    static const int wMassIdx = ParentType::pressureIdx;
    // the index of the non-wetting phase mass in solution vectors
    static const int nMassIdx = ParentType::saturationIdx;
};

/*!
 * \brief The traits for the Pn-Sw formulation of the twophase model.
 */
template <class Scalar>
class PnSwTwoPTraits : public TwoPBaseTraits<Scalar>
{
    typedef TwoPBaseTraits<Scalar> ParentType;

public:
    // set the actual formulation to pW-Sn
    static const int formulation = ParentType::pNsW;

    // the index of the wetting phase mass in solution vectors
    static const int wMassIdx = ParentType::saturationIdx;
    // the index of the non-wetting phase mass in solution vectors
    static const int nMassIdx = ParentType::pressureIdx;
};
}

#endif

