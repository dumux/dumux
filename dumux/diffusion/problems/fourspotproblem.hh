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

#ifndef FOURSPOTPROBLEM_HH
#define FOURSPOTPROBLEM_HH

#include "dumux/diffusion/diffusionproblem.hh"
#include "dumux/material/randompermeability.hh"

namespace Dune
{
//! \ingroup diffusionProblems
//! example class for diffusion problems
template<class G, class RT, class VC>
class FourSpotProblem : public DiffusionProblem<G,RT,VC>
{
    typedef typename G::ctype DT;
    enum {n=G::dimension};
    typedef typename G::Traits::template Codim<0>::Entity Entity;

    Dune::FieldMatrix<DT,n,n> K_;

public:
    FourSpotProblem(G& g, const int level, const char* name = "permeab.dat", const bool create = true,
                    TwoPhaseRelations& law = *(new DeprecatedLinearLaw), const bool cap = false)
        : grid(g), DiffusionProblem<G,RT,VC>(law, cap), permeability(g, level, name, create)
    { }

    const Dune::FieldMatrix<DT,n,n>& K (const Dune::FieldVector<DT,n>& x, const Entity& e,
                                        const Dune::FieldVector<DT,n>& xi)
    {
        K_[0][0] = K_[1][1] = 1e-10;
        //          Dune::FieldVector<DT,n> center(155);
        //          double radius = 50;
        //          if ((x-center).two_norm()<radius) K_[0][0] = K_[1][1] = 1e-6;
        //          if (x[0]>150 && x[0]<200 && x[1]<145 || x[1]>150 && x[1]<200 && x[0]<145) K_[0][0] = K_[1][1] = 1e-14;
        return K_;
        return permeability.K(e);
    }

    RT sat (const Dune::FieldVector<DT,n>& x, const Entity& e,
            const Dune::FieldVector<DT,n>& xi)
    {
        (*saturation)[grid.levelIndexSet(e.level()).index(e)];
    }

    RT q   (const Dune::FieldVector<DT,n>& x, const Entity& e,
            const Dune::FieldVector<DT,n>& xi)
    {
        return 0;
    }

    typename Dune::BoundaryConditions::Flags bctype (const Dune::FieldVector<DT,n>& x, const Entity& e,
                                                     const Dune::FieldVector<DT,n>& xi) const
    {
        if ( x[1]>285 && x[0]<15 || x[1]<15 && x[0]>285) return Dune::BoundaryConditions::dirichlet;
        return Dune::BoundaryConditions::neumann;
    }

    RT dirichlet(const Dune::FieldVector<DT,n>& x, const Entity& e,
          const Dune::FieldVector<DT,n>& xi) const
    {
        //          if (x[0]<15 && x[1]<15 || x[0]>285 && x[1]>285) return 1e7;
        if ( x[1]>285 && x[0]<15 ) return 2e5;
        return  1e5;
    }

    RT neumann(const Dune::FieldVector<DT,n>& x, const Entity& e,
          const Dune::FieldVector<DT,n>& xi) const
    {
        return 0;
    }

    LevelRandomPermeability<G> permeability;
private:
    G& grid;
public:
    Dune::BlockVector<Dune::FieldVector<RT,1> >* saturation;
};
}

#endif
