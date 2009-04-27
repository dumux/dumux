// $Id$

#ifndef DUNE_COMPUTESTORAGE_HH
#define DUNE_COMPUTESTORAGE_HH

//! \ingroup transport
//! \defgroup storage transport
/**
 * @file
 * @brief  Base class for defining the storage term of an advection-diffusion equation
 * @author Yufei Cao
 */

namespace Dune
{
//! \ingroup transport
//! Compute the storage term of the transport equation
template<class Scalar>
class ComputeStorage
{
public:
    /*! \brief Realizes the storage term.
     *
     *  \param sat first argument, usually the saturation value from one cell
     */

    virtual Scalar operator() (const Scalar sat) const
    {
       return sat;
    }

    virtual ~ComputeStorage()
    { }   
};
}
#endif
