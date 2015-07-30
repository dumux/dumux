#ifndef NACL_HH
#define NACL_HH
/*!
 * \file
 *
 * \brief A class for the NaCl properties
 */
#ifndef DUMUX_NACL_HH
#define DUMUX_NACL_HH

#include <dumux/common/exceptions.hh>
#include <dumux/material/components/component.hh>

#include <cmath>
#include <iostream>


namespace Dumux
{
/*!
 * \brief A class for the NaCl properties
 */
template <class Scalar>
class NaCl : public Component<Scalar, NaCl<Scalar> >
{
public:
    /*!
     * \brief A human readable name for the NaCl.
     */
    static const char *name()
    { 
        return "NaCl"; 
    }

    /*!
     * \brief The mass in [kg] for one mole of NaCl.
     */
    static Scalar molarMass()
    { 
        return 58.4428e-3 ; 
    } // kg/mol

    /*!
     * \brief The diffusion Coefficient of NaCl in water.
     */
    static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure)
    { 
        return 2e-9; 
    }

    /*!
     * \brief The mass density of NaCl.
     */
    static Scalar Density()
    {
        return (2165); /* 2165 kg/m³*/
    }
};

} // end namespace

#endif



#endif // NACL_HH
