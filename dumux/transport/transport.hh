// $Id$

#ifndef DUNE_TRANSPORT_HH
#define DUNE_TRANSPORT_HH

#include "dumux/transport/transportproblem.hh"

/**
 * @file
 * @brief  Base class for defining an instance of a numerical diffusion model
 * @author Bernd Flemisch
 * \defgroup transport Transport
 */

namespace Dune
{
//! \ingroup transport
//! Base class for defining an instance of a numerical transport model.
/*! An interface for defining a numerical transport model for the
 *  solution of equations of the form
 *  \f$S_t - \text{div}\, (f_\text{w}(S) \boldsymbol{v}_\text{total}) = 0\f$,
 * \f$S = g\f$ on \f$\Gamma_1\f$, and \f$S(t = 0) = S_0\f$. Here,
 * \f$S\f$ denotes the wetting phase saturation,
 * \f$\boldsymbol{v}_\text{total}\f$ the total velocity,
 * and \f$f_\text{w}\f$ the wetting phase fractional flow function.

 - Grid      a DUNE grid type
 - RT        type used for return values
 - RepresentationType   type of the vector holding the saturation values
 - VelType   type of the vector holding the velocity values

*/
template<class GridView, class Scalar, class VC, class Problem =  TransportProblem<GridView, Scalar, VC> >
class Transport
{
public:

    typedef    typename VC::ScalarVectorType RepresentationType;

    //! \brief Calculate the update vector.
    /*!
     *  \param[in]  t         time
     *  \param[out] dt        time step size
     *  \param[out] updateVec vector for hte update values
     *
     *  Calculate the update vector, i.e., the discretization
     *  of \f$\text{div}\, (f_\text{w}(S) \boldsymbol{v}_t)\f$.
     */

    virtual int update(const Scalar t, Scalar& dt, RepresentationType& updateVec, Scalar& CLFFac, bool impes) = 0;

    void initial()
    {
        initialTransport();
        return;
    }

    //! \brief Sets the initial solution \f$S_0\f$.
    virtual void initialTransport() = 0;

    //! return const reference to saturation vector
    virtual const RepresentationType& operator* () const
    {
        return transProblem.variables().saturation();
    }

    //! return reference to saturation vector
    virtual RepresentationType& operator* ()
    {
        return transProblem.variables().saturation();
    }

    virtual void postProcessUpdate(Scalar t, Scalar dt)
    {
        return;
    }

    virtual void updateMaterialLaws()
    {
        return;
    }

    virtual Problem& problem()
    {
        return transProblem;
    }

    virtual void vtkout(const char* name, int k)const = 0;

    //! always define virtual destructor in abstract base class
    virtual ~Transport ()
    {}

    /*! @brief constructor
     *  @param gridView a DUNE gridview object
     *  @param prob an object of class TransportProblem or derived
     */
    Transport(const GridView& gridView, Problem& problem)
        : gridView(gridView), transProblem(problem)
    {}

protected:
    const GridView& gridView;
    Problem& transProblem; //!< problem data
};

}
#endif
