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

#ifndef LEVELHETPROBLEM_HH
#define LEVELHETPROBLEM_HH

#include "dumux/diffusion/diffusionproblem.hh"
#include "dumux/material/randompermeability.hh"
#include <dune/istl/bvector.hh>

namespace Dune {
//! \ingroup diffusionProblems
//! example class for diffusion problems
template<class G, class RT, class VC> class LevelHetProblem :
        public DeprecatedDiffusionProblem<G,RT,VC> {
    typedef typename G::ctype DT;
    enum {n=G::dimension};
    typedef typename G::Traits::template Codim<0>::Entity Entity;

public:
    LevelHetProblem(VC& variableobj, const int level,
                    const char* name = "permeab.dat", const bool create = true,
                    DeprecatedTwoPhaseRelations& law = *(new DeprecatedLinearLaw), const bool cap = false) :
        DeprecatedDiffusionProblem<G, RT, VC>(variableobj, law, cap),
        permeability(variableobj.grid, level, name, create) {
    }

    Dune::FieldMatrix<DT,n,n>& K(const Dune::FieldVector<DT,n>& x,
                                 const Entity& e, const Dune::FieldVector<DT,n>& xi) {

        return  permeability.K(e);
    }

    RT q(const Dune::FieldVector<DT,n>& x, const Entity& e,
         const Dune::FieldVector<DT,n>& xi) {
        return 0;
    }

    typename Dune::BoundaryConditions::Flags bctype(
                                                    const Dune::FieldVector<DT,n>& x, const Entity& e,
                                                    const Dune::FieldVector<DT,n>& xi) const {
        if (x[0] > 300-1e-6|| x[0] < 1e-6)
            return Dune::BoundaryConditions::dirichlet;
        // all other boundaries
        return Dune::BoundaryConditions::neumann;
    }

    RT dirichlet(const Dune::FieldVector<DT,n>& x, const Entity& e,
         const Dune::FieldVector<DT,n>& xi) const {
        return (x[0] < 1e-6) ? 2e6 : 1e6;
    }

    RT gSat(const FieldVector<DT,n>& x, const Entity& e,
            const FieldVector<DT,n>& xi) const {
        if (x[0] < 1e-6)
            return 0.8;
        else
            return 0.2;
    }

    RT neumann(const Dune::FieldVector<DT,n>& x, const Entity& e,
         const Dune::FieldVector<DT,n>& xi) const {
        if (x[0]<1e-6)
            return 1e-4;
        return 0;
    }

    LevelRandomPermeability<G> permeability;
};
}

#endif
