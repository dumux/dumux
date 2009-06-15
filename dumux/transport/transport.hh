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
#ifndef DUNE_TRANSPORT_HH
#define DUNE_TRANSPORT_HH

#include "dumux/transport/transportproblem.hh"

/**
 * @file
 * @brief  Base class for defining an instance of a numerical transport model
 * @author Bernd Flemisch, Markus Wolff
 */

/*!
 * \ingroup fracflow
 * \defgroup transport Transport
 */

namespace Dune
{
//! \ingroup transport
//! Base class for defining an instance of a numerical transport model.
/*! An interface for defining a numerical transport model for the
 *  solution of equations of the form
 *  \f[
 *    \frac{\partial S}{\partial t} + \text{div}\, \boldsymbol{v_{\alpha}} = 0,
 *  \f]
 *  where \f$S\f$ denotes a phase saturation and \f$\boldsymbol{v_{\alpha}}\f$ is a phase velocity.

 Template parameters are:

 - GridView      a DUNE gridview type
 - Scalar        type used for scalar quantities
 - VC            type of a class containing different variables of the model
 - Problem       class defining the physical problem

*/
template<class GridView, class Scalar, class VC, class Problem =  TransportProblem<GridView, Scalar, VC> >
class Transport
{
public:
    typedef    typename VC::ScalarVectorType RepresentationType;//!< Data type for a Vector of Scalars

    //! Calculate the update vector.
    /*!
     *  \param[in]  t         time
     *  \param[in] dt         time step size
     *  \param[in] updateVec  vector for the update values
     *  \param[in] CLFFac     security factor for the time step criterion (0 < CLFFac <= 1)
     *  \param[in] impes      variable is true if an impes algorithm is used and false if the transport part is solved independently
     *
     *  Calculate the update vector
     */
    virtual int update(const Scalar t, Scalar& dt, RepresentationType& updateVec, Scalar& CLFFac, bool impes) = 0;

    //! Call the initialization function.
    /*!
     *  Call the initialization function -> called from Dune::TimeLoop if the transport part is solved independently
     */
    void initial()
    {
        initialTransport();
        return;
    }

    //! Sets the initial solution \f$S_0\f$.
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

    //! Start a post-processing procedure at the end of a timestep.
    /*!
     *  \param t time
     *  \param dt time step
     *
     *  If an explicit Euler time discretization is used this function will be called
     *  at the end of each time step.
     */
    virtual void postProcessUpdate(Scalar t, Scalar dt)
    {
        return;
    }

    //! Update the values of the material laws and constitutive relations.
    /*!
     *  Constitutive relations like capillary pressure-saturation relationships, mobility-saturation relationships... are updated and stored in the variable class
     *  of type Dune::VariableClass2P. The update has to be done when new saturation are available.
     */
    virtual void updateMaterialLaws(RepresentationType& saturation, bool iterate)
    {
        return;
    }

    //! Returns a reference to the problem
    virtual Problem& problem()
    {
        return transProblem;
    }

    //! \brief Write data files
    /*!
     *  \param name file name
     *  \param k format parameter
     */
    virtual void vtkout(const char* name, int k)const = 0;

    //! always define virtual destructor in abstract base class
    virtual ~Transport ()
    {}

    /*! @brief Constructs a Transport object
     *  @param gridView a DUNE gridview object
     *  @param prob an object of class TransportProblem or derived
     */
    Transport(const GridView& gridView, Problem& problem)
        : gridView(gridView), transProblem(problem)
    {}

protected:
    const GridView& gridView; //!< object of type Dune::GridView
    Problem& transProblem; //!< problem data
};

}
#endif
