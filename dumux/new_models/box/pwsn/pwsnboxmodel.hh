/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
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
#ifndef DUMUX_PWSN_BOX_MODEL_HH
#define DUMUX_PWSN_BOX_MODEL_HH

#include <dumux/new_models/box/boxmodel.hh>
#include <dumux/new_models/box/pwsn/pwsnboxjacobian.hh>
#include <dumux/new_models/box/pwsn/pwsnboxtraits.hh>

#include <dumux/auxiliary/apis.hh>

#include <dune/istl/operators.hh>
#include <dune/disc/operators/p1operator.hh>

#include <boost/format.hpp>

namespace Dune
{

    /*!
     * \brief The base class for the BOX hybrid finite element/finite volume discretization model
     */
    template<class ProblemT>
    class PwSnBoxModel : public BoxModel< PwSnBoxTraits<typename ProblemT::DomainTraits::Scalar,
                                                        typename ProblemT::DomainTraits::Grid>,
                                          ProblemT, 
                                          PwSnBoxJacobian<ProblemT, 
                                                           PwSnBoxTraits<typename ProblemT::DomainTraits::Scalar,
                                                                         typename ProblemT::DomainTraits::Grid> > >
    {
        typedef typename ProblemT::DomainTraits::Grid   Grid;
        typedef typename ProblemT::DomainTraits::Scalar Scalar;

    public:
        typedef PwSnBoxTraits<Scalar, Grid>             BoxTraits;
        
    private:
        typedef PwSnBoxJacobian<ProblemT, BoxTraits>    PwSnLocalJacobian;
        typedef BoxModel<BoxTraits,
                         ProblemT, 
                         PwSnLocalJacobian>  ParentType;

    public:
        typedef NewNewtonMethod<ParentType> NewtonMethod;

        PwSnBoxModel(ProblemT &prob)
            : ParentType(prob, _pwSnLocalJacobian),
              _pwSnLocalJacobian(prob, false)
            {
                Api::require<Api::BasicDomainTraits, typename ProblemT::DomainTraits>();
            }

    private:
        // calculates the jacobian matrix at a given position
        PwSnLocalJacobian  _pwSnLocalJacobian;
    };
}

#endif
