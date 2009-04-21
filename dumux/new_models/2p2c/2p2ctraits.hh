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

#ifndef DUMUX_2P2CTRAITS_HH
#define DUMUX_2P2CTRAITS_HH

namespace Dune
{
///////////////////////////////////////////////////////////////////////////
// two-phase two-component traits (central place for names and
// indices required by the TwoPTwoCBoxJacobian and TwoPTwoCBoxModel)
///////////////////////////////////////////////////////////////////////////
/*!
 * \brief The 2P-2C specific traits.
 */
template <class Scalar>
class TwoPTwoCBaseTraits
{
public:
    enum {
        numEq         = 2,  //!< Number of primary variables / equations
        numPhases     = 2,  //!< Number of fluid phases
        numComponents = 2   //!< Number of fluid components within a phase
    };
    enum { // Primary variable indices
        pressureIdx = 0,     //!< Index for wetting/non-wetting phase pressure (depending on formulation) in a solution vector
        switchIdx   = 1,     //!< Index of the either the saturation or the mass fraction of the non-wetting/wetting phase
    };
    enum { // present phases
        nPhaseOnly = 0, //!< Only the non-wetting phase is present
        wPhaseOnly = 1, //!< Only the wetting phase is present
        bothPhases = 2  //!< Both phases are present
    };
    enum { // formulation
        pWsN = 0, //!< Pw and Sn as primary variables
        pNsW = 1,  //!< Pn and Sw as primary variables
    };

    typedef FieldVector<Scalar, numPhases>  PhasesVector;

    /*!
     * \brief Data which is attached to each vertex of an
     *        element. These quantities are coincidental with the
     *        averaged quantities inside a FV box.
     */
    struct VertexData
    {
        PhasesVector saturation;
        PhasesVector pressure;
        Scalar pC;
        Scalar porosity;
        PhasesVector mobility;
        PhasesVector density;
        PhasesVector diffCoeff; // diffusion coefficients for the phases
        Dune::FieldMatrix<Scalar, numComponents, numPhases> massfrac;
    };

    /*!
     * \brief Data at the each vertex of an element.
     */
    typedef std::vector<VertexData> ElementData;
};


template <class Scalar>
class TwoPTwoCPwSnTraits : public TwoPTwoCBaseTraits<Scalar>
{
    typedef TwoPTwoCBaseTraits<Scalar>     ParentT;

public:
    typedef typename ParentT::PhasesVector PhasesVector;

    enum {
        numEq         = ParentT::numEq,  //!< Number of primary variables / equations
        numPhases     = ParentT::numPhases,  //!< Number of fluid phases
        numComponents = ParentT::numComponents,   //!< Number of fluid components within a phase

        // Primary variable indices
        pressureIdx = ParentT::pressureIdx,   //!< Idx for either wetting or non-wetting phase pressure (depending on formulation) in a solution vector
        switchIdx = ParentT::switchIdx,       //!< Idx for the wetting phase quantity
        // present phases
        nPhaseOnly = ParentT::nPhaseOnly, //!< Only the non-wetting phase is present
        wPhaseOnly = ParentT::wPhaseOnly, //!< Only the wetting phase is present
        bothPhases = ParentT::bothPhases,  //!< Both phases are present
        // formulation
        pWsN = ParentT::pWsN, //!< Pw and Sn as primary variables
        pNsW = ParentT::pNsW,  //!< Pn and Sw as primary variables

        formulation      = pWsN,
        wPhase           = 0,
        nPhase           = 1,

        wComp            = 0,
        nComp            = 1
    };
};

template <class Scalar>
class TwoPTwoCPnSwTraits : public TwoPTwoCBaseTraits<Scalar>
{
    typedef TwoPTwoCBaseTraits<Scalar>     ParentT;

public:
    typedef typename ParentT::PhasesVector PhasesVector;

    enum {
        numEq         = ParentT::numEq,  //!< Number of primary variables / equations
        numPhases     = ParentT::numPhases,  //!< Number of fluid phases
        numComponents = ParentT::numComponents,   //!< Number of fluid components within a phase

        // Primary variable indices
        pressureIdx = ParentT::pressureIdx,   //!< Idx for either wetting or non-wetting phase pressure (depending on formulation) in a solution vector
        switchIdx = ParentT::switchIdx,       //!< Idx for the wetting phase quantity
        // present phases
        nPhaseOnly = ParentT::nPhaseOnly, //!< Only the non-wetting phase is present
        wPhaseOnly = ParentT::wPhaseOnly, //!< Only the wetting phase is present
        bothPhases = ParentT::bothPhases,  //!< Both phases are present
        // formulation
        pWsN = ParentT::pWsN, //!< Pw and Sn as primary variables
        pNsW = ParentT::pNsW, //!< Pn and Sw as primary variables

        formulation      = pNsW,
        wPhase           = 1,
        nPhase           = 0,

        wComp            = 1,
        nComp            = 0

    };

};

///////////////////////////////////////////////////////////////////////////
// two-phase two-component traits (central place for names and
// indices required by the TwoPTwoCNIBoxJacobian and TwoPTwoCNIBoxModel)
///////////////////////////////////////////////////////////////////////////
/*!
 * \brief The 2P-2C specific traits.
 */
template <class Scalar,
          class BaseTraits = TwoPTwoCPwSnTraits<Scalar> >
class TwoPTwoCNITraits : public BaseTraits
{
public:
    enum {
        numEq         = 3,  //!< Number of primary variables
        numPhases     = 2,  //!< Number of fluid phases
        numComponents = 2   //!< Number of fluid components within a phase
    };
    enum { // Primary variable indices
        pressureIdx = BaseTraits::pressureIdx,
        switchIdx = BaseTraits::switchIdx,
        temperatureIdx = 2
    };
    enum { // Phase Indices
        wPhase = BaseTraits::wPhase,
        nPhase = BaseTraits::nPhase
    };
    enum { // Component indices
        wComp = BaseTraits::wComp,
        nComp = BaseTraits::nComp
    };
    enum { // present phases
        nPhaseOnly = BaseTraits::nPhaseOnly,
        wPhaseOnly = BaseTraits::wPhaseOnly,
        bothPhases = BaseTraits::bothPhases
    };

    typedef typename BaseTraits::PhasesVector PhasesVector;

    struct VariableVertexData : public BaseTraits::VariableVertexData
    {
        PhasesVector intenergy;
        PhasesVector enthalpy;
        Scalar       lambda;
    };
};



}

#endif
