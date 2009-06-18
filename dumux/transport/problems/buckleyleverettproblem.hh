// $Id:$
/*****************************************************************************
 *   Copyright (C) <YEARS> by <ADD_AUTHOR_HERE>                              *
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

#ifndef DUNE_BUCKLEYLEVERETTPROBLEM_HH
#define DUNE_BUCKLEYLEVERETTPROBLEM_HH

#include "dumux/transport/transportproblem_deprecated.hh"

namespace Dune
{
//! \ingroup transportProblems
//! @brief example class for a transport problem
template<class G, class RT, class VC>
class BuckleyLeverettProblem
    : public DeprecatedTransportProblem<G, RT,VC> {

    typedef typename G::ctype DT;
    enum {n=G::dimension, m=1, blocksize=2*G::dimension};
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef Dune::FieldVector<double, n> R1;

private:
    DT left;
    DT right;

public:
    BoundaryConditions::Flags bctype (const FieldVector<DT,n>& x, const Entity& e,
                                      const FieldVector<DT,n>& xi) const
    {
        if (x[0] > right-1E-8 || x[0] < left+1e-8)
            return Dune::BoundaryConditions::dirichlet;
        else
            return Dune::BoundaryConditions::neumann;
    }

    RT dirichlet (const FieldVector<DT,n>& x, const Entity& e,
                  const FieldVector<DT,n>& xi) const
    {
        if (x[0] < left+1e-8)
            return 0.8;
        else
            return 0.2;
    }

    RT initSat (const FieldVector<DT,n>& x, const Entity& e,
                const FieldVector<DT,n>& xi) const
    {
        return 0.2;
    }


    BuckleyLeverettProblem(VC& variableobj, DeprecatedTwoPhaseRelations& law = *(new DeprecatedLinearLaw),
                           const int level = 0, const bool cap = false)
        : DeprecatedTransportProblem<G, RT, VC>(variableobj,law, cap), left((variableobj.grid.lowerLeft())[0]), right((variableobj.grid.upperRight())[0])
    {}
};

}
#endif
