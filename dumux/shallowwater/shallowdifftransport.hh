// $Id$

#ifndef DUNE_SHALLOWTRANSPORT_HH
#define DUNE_SHALLOWTRANSPORT_HH

#include "dumux/shallowwater/shallowproblembase.hh" //Hier wird Basisklasse des Problems eingebunden

/**
 * @file
 * @brief  Base class for defining an instance of a numerical diffusion model
 * @author Bernd Flemisch
 * \defgroup transport Transport
 */

namespace Dune
{
/*- Grid      a DUNE grid type
  - RT        type used for return values
  - RepresentationType   type of the vector holding the saturation values
  - VelType   type of the vector holding the velocity values */

template<class Grid, class Scalar, class VC>
class ShallowTransport {
public:

    typedef BlockVector< Dune::FieldVector<Scalar,1> > RepresentationType;
    typedef BlockVector< Dune::FieldVector<Dune::FieldVector<Scalar,1>, size > > UpdatevectorType;
    ShallowProblemBase<Grid, Scalar, VC>& transproblem; //!< problem data

    //! \brief Calculate the update vector.
    /*!
     *  \param[in]  t         time
     *  \param[out] dt        time step size
     *  \param[out] updateVec vector for hte update values
     */
    virtual int update(const Scalar t, Scalar& dt, UpdatevectorType& updateVec, Scalar& CLFFac) = 0;

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
        return transproblem.variables.waterdepth;
    }

    //! return reference to saturation vector
    virtual RepresentationType& operator* ()
    {
        return transproblem.variables.wDepth;
    }

    virtual void vtkout(const char* name, int k) const {
        transproblem.variables.vtkout(name, k);
        return;
    }

    //! always define virtual destructor in abstract base class
    virtual ~ShallowTransport () {}

    /*! @brief constructor
     *  @param g a DUNE grid object
     *  @param prob an object of class TransportProblem or derived
     */
    ShallowTransport(const Grid& g, ShallowProblemBase<Grid, Scalar, VC>& prob)
        : grid(g), transproblem(prob), level_(g.maxLevel())
    { }

    /*! @brief constructor
     *
     *  This constructor gives the additional possibility to specify a grid level on which
     *  the transport equation shall be solved. This especially important for multiscale methods.
     *  @param g a DUNE grid object
     *  @param prob an object of class TransportProblem or derived
     *  @param lev the grid level on which the Transport equation is to be solved.
     */
    ShallowTransport(const Grid& g, ShallowProblemBase<Grid, Scalar, VC>& prob, int lev)
        : transproblem(prob), grid(g), level_(lev)
    { }

    //! returns the level on which the transport eqution is solved.
    int level() const
    {
        return level_;
    }

    const Grid& g;

protected:
    const int level_;
};

}
#endif
