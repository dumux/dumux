// $Id$

#ifndef DUMUX_1P2CTRAITS_HH
#define DUMUX_1P2CTRAITS_HH

namespace Dune
{

///////////////////////////////////////////////////////////////////////////
// one-phase two-component traits (central place for names and
// indices required by the OnePTwoCBoxJacobian and OnePTwoCBoxModel)
///////////////////////////////////////////////////////////////////////////
/*!
 * \brief The 1P-2C specific traits.
 */
template <class Scalar>
class OnePTwoCTraits
{
public:
    enum {
        numEq         = 2,  //!< Number of primary variables / equations
        numPhases     = 1,  //!< Number of fluid phases
        numComponents = 2   //!< Number of fluid components within a phase
    };
    enum { // Primary variable indices
        konti = 0,       // pressure in the solution vector
        transport = 1,       //mole fraction in the solution vector
    };

    typedef FieldVector<Scalar, numPhases>         PhasesVector;
};



}

#endif
