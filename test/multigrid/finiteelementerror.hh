#ifndef FINITEELEMENTERROR_HH
#define FINITEELEMENTERROR_HH

#include <dune/disc/functions/functions.hh>

// compute L2 norm of finite element interpolation error
template<class G, class RT, int m, class Functor>
double L2Error (G& grid, const Dune::GridFunction<G, RT, m> &gf, Functor f, int p)
{
    // first we extract the dimensions of the grid
    const int dim = G::dimension;

    // type used for coordinates in the grid
    typedef typename G::ctype ct;

    // the grid has an iterator providing the access to
    // all elements (better codim 0 entities) on a grid level
    // Note the use of the typename and template keywords.
    typedef typename G::template Codim<0>::LeafIterator LeafIterator;

    // iterate through all entities of codim 0 on the given level
    double sum = 0.0;                                                     // the final result to be summed
    int order = p;
    LeafIterator eendit = grid.template leafend<0>();                     // preallocation for efficiency
    for (LeafIterator it = grid.template leafbegin<0>(); it!=eendit; ++it)// traverse leaf elements of the grid
        {
            Dune::GeometryType gt = it->geometry().type();                       // extract type of element
            p = order;
            if (gt.isPyramid()) p = 1; // a hack ...
            for (size_t i=0; i<Dune::QuadratureRules<ct,dim>::rule(gt,p).size(); ++i) // run through all quadrature points
                {
                    const Dune::FieldVector<ct,dim>& ippos = Dune::QuadratureRules<ct,dim>::rule(gt,p)[i].position();
                    double exact = f(it->geometry().global(ippos));
                    double approx = gf.evallocal(0,*it,ippos);
                    double weight = Dune::QuadratureRules<ct,dim>::rule(gt,p)[i].weight();
                    double refvolume = Dune::ReferenceElements<ct,dim>::general(gt).volume();
                    double detjac = it->geometry().integrationElement(ippos);
                    sum += (exact-approx)*(exact-approx)*weight*refvolume*detjac;
                }
        }

    return sqrt(sum);
}

// compute H1 norm of finite element interpolation error
template<class G, class RT, int m, class Functor>
double H1Error (G& grid, const Dune::DifferentiableGridFunction<G, RT, m> &gf, Functor f, int p)
{
    // first we extract the dimensions of the grid
    const int dim = G::dimension;
    const int dimworld = G::dimensionworld;

    // type used for coordinates in the grid
    typedef typename G::ctype ct;

    // the grid has an iterator providing the access to
    // all elements (better codim 0 entities) on a grid level
    // Note the use of the typename and template keywords.
    typedef typename G::template Codim<0>::LeafIterator LeafIterator;

    // iterate through all entities of codim 0 on the given level
    double sum = 0.0;                                                      // the final result to be summed
    LeafIterator eendit = grid.template leafend<0>();                     // preallocation for efficiency
    for (LeafIterator it = grid.template leafbegin<0>(); it!=eendit; ++it)// traverse leaf elements of the grid
        {
            Dune::GeometryType gt = it->geometry().type();                       // extract type of element

            for (int i=0; i<Dune::QuadratureRules<ct,dim>::rule(gt,p).size(); ++i) // run through all quadrature points
                {
                    const Dune::FieldVector<ct,dim>& ippos = Dune::QuadratureRules<ct,dim>::rule(gt,p)[i].position();
                    double weight = Dune::QuadratureRules<ct,dim>::rule(gt,p)[i].weight();
                    double refvolume = Dune::ReferenceElements<ct,dim>::general(gt).volume();
                    double detjac = it->geometry().integrationElement(ippos);

                    // function value
                    double exact = f(it->geometry().global(ippos));
                    double approx = gf.evallocal(0,*it,ippos);
                    sum += (exact-approx)*(exact-approx)*weight*refvolume*detjac;

                    // derivatives
                    Dune::FieldVector<RT,dim> gexact = f.grad(it->geometry().global(ippos));
                    for (int dir=0; dir<dim; dir++)
                        {
                            Dune::FieldVector<int,dim> d(0);
                            d[dir] = 1;
                            double gapprox = gf.derivativelocal(0,d,*it,ippos);
                            sum += (gexact[dir]-gapprox)*(gexact[dir]-gapprox)*weight*refvolume*detjac;
                        }
                }
        }

    return sqrt(sum);
}

#endif
