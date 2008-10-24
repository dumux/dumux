// $Id$

#ifndef RHS_HH
#define RHS_HH
#include "testfunctions.hh"

template<class Grid>
class RightHandSide
{
  enum{dim=Grid::dimension};
  typedef typename Grid::ctype ct;
public:
  
   RightHandSide(ExactSolution<ct,dim>& ex): exact(ex){}
  
   typedef Dune::FieldVector<ct,dim> Point;
  double rhsValue(int variable,  Point& global, const Point& local) const;
 
protected:
   ExactSolution<ct,dim>& exact;
};



#include"rhs.cc"

#endif
