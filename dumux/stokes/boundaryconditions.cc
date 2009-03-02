// $Id$

#include "config.h"
#include"boundaryconditions.hh"

template<class Grid>
double
DirichletBoundary<Grid>::dirichletValue(int variable,  const Point& global, Point& local) const
{
    return exact.velocity(variable,global);
}





