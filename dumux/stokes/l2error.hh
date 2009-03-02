// $Id$

#include<dumux/stokes/dgstokes.hh>


template <class G,int v_order, int p_order>
double
Dune::DGFiniteElementMethod<G,v_order,p_order>::evaluateL2error(int variable, const Entity& element, const LocalVectorBlock& xe)const

{
    // stokes system has dim+1 variables (dim velocity comps and 1 pressure)
    double error[dim+1];
    error[variable]=0.0;
    Dune::FieldVector<ctype, dim> qp_loc(0.0);
    Dune::FieldVector<ctype, dim> qp_glob(0.0);
    Dune::GeometryType gt = element.type();
    // #warning fixed quadrature order
    int qord=12;
    //int eid = grid.levelIndexSet(grid.maxLevel()).index(element);

    for (unsigned int qp=0;qp<Dune::QuadratureRules<ctype,dim>::rule(gt,qord).size();++qp)
    {
        qp_loc = Dune::QuadratureRules<ctype,dim>::rule(gt,qord)[qp].position();
        qp_glob =element.geometry().global(qp_loc);
        double weight = Dune::QuadratureRules<ctype,dim>::rule(gt,qord)[qp].weight();
        double detjac = element.geometry().integrationElement(qp_loc);
        if (variable<dim)
        {
            error[variable]+=weight*detjac
                *((problem_.velocity(qp_glob))[variable]-evaluateSolution(variable,element,qp_loc,xe))
                *((problem_.velocity(qp_glob))[variable]-evaluateSolution(variable,element,qp_loc,xe));
        }
        if(variable==dim)
        {
            error[variable]+=weight*detjac
                *(problem_.pressure(qp_glob)-evaluateSolution(variable,element,qp_loc,xe))
                *(problem_.pressure(qp_glob)-evaluateSolution(variable,element,qp_loc,xe));
        }

    }
    return error[variable];

}


template<class G,int v_order, int p_order>
double Dune::DGStokes<G,v_order,p_order>::l2errorStokesSystem(int variable) const
{
    // stokes system has dim+1 variables (dim velocity comps and 1 pressure)
    double error[dim+1];
    error[variable]= 0.0;
    ElementLevelIterator it = grid.template lbegin<0>(level);
    ElementLevelIterator itend = grid.template lend<0>(level);
    for (; it != itend; ++it)
    {
        int eid = grid.levelIndexSet(level).index(*it);
        error[variable]+=dgfem.evaluateL2error(variable,*it,solution[eid]);
    }
    return sqrt(error[variable]);
}


