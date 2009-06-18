// $Id:$
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
/*!
 * \file
 *
 * \brief Base class for all problems which use the non-isothermal
 *        two-phase box model
 */
#ifndef DUMUX_2PNI_BOX_PROBLEM_HH
#define DUMUX_2PNI_BOX_PROBLEM_HH

#include <dumux/boxmodels/tags.hh>
#include <dumux/boxmodels/boxscheme/boxproblem.hh>
#include <dumux/material/twophaserelations.hh>

#include <dumux/boxmodels/2p/2pboxproblem.hh>


namespace Dune
{
/*!
 * \ingroup TwoPNIProblems
 * \brief  Base class for all problems which use the non-isothermal
 *         two-phase box model.
 *
 * \todo Please doc me more!
 */
template<class TypeTag, class Implementation>
class TwoPNIBoxProblem : public TwoPBoxProblem<TypeTag, Implementation>
{
    typedef TwoPBoxProblem<TypeTag, Implementation> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

public:
    TwoPNIBoxProblem(const GridView &gridView)
        : ParentType(gridView)
    {
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This is a non-isothermal model, so this method just throws an
     * exception. This method MUST NOT be overwritten by the actual
     * problem.
     */
    Scalar temperature() const
    { DUNE_THROW(Dune::Exception, "temperature() method called for a 2p2cni problem"); };
  
    // \}
};

}

#endif
