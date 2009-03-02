// $Id$

#ifndef BOUNDARYCONDITIONS_HH
#define BOUNDARYCONDITIONS_HH

#include"testfunctions.hh"

/** \todo Please doc me! */

template<class Grid>
class DirichletBoundary
{

    enum{dim=Grid::dimension};
    typedef typename Grid::ctype ct;
public:
    DirichletBoundary(ExactSolution<ct,dim>& ex): exact(ex){}
    typedef Dune::FieldVector<ct,dim> Point;
    double dirichletValue(int comp,const Point& global, Point& local) const;
protected:
    ExactSolution<ct,dim>& exact;
};

#include"boundaryconditions.cc"
#endif

