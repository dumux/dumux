#ifndef DUNE_TRANSPORTPROBLEM_HH
#define DUNE_TRANSPORTPROBLEM_HH

#include"dumux/stokes/stokestransportproblem.hh"

namespace Dune {

template<class Grid, class Scalar>
class TransportProblem : public StokesTransportProblem<Grid, Scalar>
{
public:
    enum {dim=Grid::dimension, numEq=Grid::dimension+2};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::LeafGridView::IntersectionIterator IntersectionIterator;


    virtual FieldVector<Scalar,numEq> initial (const FieldVector<Scalar,dim>& x, const Element& e,
                                               const FieldVector<Scalar,dim>& xi) const
    {
        FieldVector<Scalar,numEq> result(0);

        result[0] = 1.0e-3;

        if (x[0] > 0.4 && x[0] < 0.6 && x[1] > 0.4 && x[1] < 0.6 )
            result[dim] = 1;

        return result;
    }

    virtual FieldVector<Scalar,numEq> q(const FieldVector<Scalar,dim>& x, const Element& e,
                                        const FieldVector<Scalar,dim>& xi) const
    {
        FieldVector<Scalar,numEq> result(0);

        return result;
    }

    virtual BoundaryConditions::Flags bctype (const FieldVector<Scalar,dim>& x, const Element& e,
                                              const IntersectionIterator& intersectionIt,
                                              const FieldVector<Scalar,dim>& xi) const
    {
        return BoundaryConditions::dirichlet;
    }

    virtual FieldVector<Scalar,dim+1> g(const FieldVector<Scalar,dim>& x, const Element& e,
                                        const IntersectionIterator& intersectionIt,
                                        const FieldVector<Scalar,dim>& xi) const
    {
        FieldVector<Scalar,dim+1> result(0);

        result[0] = 1.0e-3;

        return result;
    }

    virtual FieldVector<Scalar,dim+1> J(const FieldVector<Scalar,dim>& x, const Element& e,
                                        const IntersectionIterator& intersectionIt,
                                        const FieldVector<Scalar,dim>& xi)
    {
        FieldVector<Scalar,dim+1> result(0);
        return result;
    }

    virtual FieldMatrix<Scalar,dim,dim> D (const FieldVector<Scalar,dim>& x, const Element& e,
                                           const FieldVector<Scalar,dim>& xi) const
    {
        FieldMatrix<Scalar,dim,dim> res(0);

        for (int kx=0; kx<dim; kx++)
            for (int ky=0; ky<dim; ky++)
                if (kx == ky)
                    res[kx][ky] = 0;

        return res;
    }
    //------------------additional----------------
    virtual Scalar Qg(const FieldVector<Scalar,dim>& x, const Element& e,
                      const FieldVector<Scalar,dim>& xi) const
    {
        Scalar result = 0;
        return result;
    }

    virtual Scalar gravity() const
    {
        return 0;//-9.8;
    }

    virtual Scalar density(const FieldVector<Scalar,dim>& x, const Element& e,
                           const FieldVector<Scalar,dim>& xi) const
    {
        Scalar result = 1;
        return result;
    }
    //------------------additional----------------


    /*
      virtual FieldMatrix<Scalar, dim, dim> velocityGradient(const FieldVector<Scalar,dim>& x) const
      {
      FieldMatrix<Scalar, dim, dim> result(0);

      return result;
      }
    */
    virtual Fluid& gasPhase () const
    {
        return gasPhase_;
    }

    virtual MultiComp& multicomp () const
    {
        return multicomp_;
    }

    TransportProblem()
    {
        pi = 4.0*atan(1.0);
    }

    TransportProblem(Fluid& gasPhase, MultiComp& multicomp = *(new CWaterAir))
        :
        StokesTransportProblem<Grid,Scalar>(gasPhase,multicomp),
        gasPhase_(gasPhase),
        multicomp_(multicomp)
    {
        pi = 4.0*atan(1.0);
    }


protected:
    double pi;
    Fluid& gasPhase_;
    MultiComp& multicomp_;
};

}
#endif
