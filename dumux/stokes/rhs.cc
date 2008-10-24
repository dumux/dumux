#include"rhs.hh"

template<class Grid>
double
RightHandSide<Grid>::rhsValue(int variable,Point& global,const Point& local) const
{
  return exact.rhsvalue(variable,global);
}
