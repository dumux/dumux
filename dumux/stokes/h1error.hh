// $Id$

#include<dune/disc/stokes/dgstokes.hh>


template <class G,int v_order, int p_order>
double
Dune::DGFiniteElementMethod<G,v_order,p_order>::evaluateH1error(int variable, const Entity& element, const LocalVectorBlock& xe)const

{
  // stokes system has dim+1 variables (dim velocity comps and 1 pressure)
  double error[dim+1],error1[dim+1],error2[dim+1];
  error[variable]=0.0;
  Dune::FieldVector<ctype, dim> qp_loc(0.0);
  Dune::FieldVector<ctype, dim> qp_glob(0.0);
  Dune::GeometryType gt = element.type();
  // #warning fixed quadrature order
  int qord=6;

  for (int qp=0;qp<Dune::QuadratureRules<ctype,dim>::rule(gt,qord).size();++qp)
	{
	  qp_loc = Dune::QuadratureRules<ctype,dim>::rule(gt,qord)[qp].position();

	  qp_glob =element.geometry().global(qp_loc);
	  //std::cout<<"qp_loc: "<<qp_loc<<" qp_glob: "<<qp_glob<<std::endl;
	  double weight = Dune::QuadratureRules<ctype,dim>::rule(gt,qord)[qp].weight();
	  double detjac = element.geometry().integrationElement(qp_loc);
	  if (variable<dim)
		{
		  error1[variable]=((problem_.velocity(qp_glob))[variable]-evaluateSolution(variable,element,qp_loc,xe))
		    *((problem_.velocity(qp_glob))[variable]-evaluateSolution(variable,element,qp_loc,xe));

		  error2[variable]=((problem_.velocityGradient(qp_glob))[variable]-evaluateGradient(variable,element,qp_loc,xe)).two_norm2();
		  error[variable]+=weight*detjac*(error1[variable]+error2[variable]);


		}
	}
  return error[variable];

}


template<class G,int v_order, int p_order>
double Dune::DGStokes<G,v_order,p_order>::h1errorStokesSystem(int variable) const
{
  // stokes system has dim+1 variables (dim velocity comps and 1 pressure)
  double error[dim+1];
  error[variable]= 0.0;
    ElementLevelIterator it = grid.template lbegin<0>(level);
  ElementLevelIterator itend = grid.template lend<0>(level);
  //std::cout<<"level:::"<<level;
  for (; it != itend; ++it)
	{
	  int eid = grid.levelIndexSet(level).index(*it);
	  //std::cout<<" eid: "<<eid<<std::endl;
	  error[variable]+=dgfem.evaluateH1error(variable,*it,solution[eid]);
	}
  return sqrt(error[variable]);
}


