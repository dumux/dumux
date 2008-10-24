#ifndef TESTFUNCTIONS_HH
#define TESTFUNCTIONS_HH
#include<dune/common/misc.hh>
#include<dune/common/fvector.hh>

template<class ct,int dim>
class ExactSolution
{
public:
  typedef Dune::FieldVector<ct,dim> Point;
  typedef Dune::FieldVector< ct, dim> Gradient;
  ExactSolution(){}
 
  virtual ct velocity(int variable,const Point & global) const = 0;
  virtual Gradient velocityGradient(int variable,const Point & global) const = 0;
  virtual ct pressure(const Point & global) const = 0;
  virtual ct rhsvalue(int variable, const Point& global) const =0;
 
  virtual ~ExactSolution() {}
};

#include"testfunctions.cc"
#endif


