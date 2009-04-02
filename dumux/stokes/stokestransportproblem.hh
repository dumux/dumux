#ifndef DUNE_STOKESTRANSPORTPROBLEM_HH
#define DUNE_STOKESTRANSPORTPROBLEM_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/istl/bvector.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/grid/utility/intersectiongetter.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dumux/material/twophaserelations.hh>
#include<dumux/material/phaseproperties/phaseproperties_waterair.hh>
#include "dumux/material/multicomponentrelations.hh"

namespace Dune {

template<class Grid, class Scalar> class StokesTransportProblem {

public:
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;
    enum {dim=Grid::dimension, numEq=Grid::dimension+2};

    virtual FieldVector<Scalar,numEq> q(const FieldVector<Scalar,dim>& x, const Element& e,
                                        const FieldVector<Scalar,dim>& xi) const
    {
        DUNE_THROW(NotImplemented, "no q specified, but requested");

        FieldVector<Scalar,numEq> result(0);
        return result;
    }

    virtual BoundaryConditions::Flags bctype (const FieldVector<Scalar,dim>& x, const Element& e,
                                              const IntersectionIterator& intersectionIt,
                                              const FieldVector<Scalar,dim>& xi) const = 0;

    virtual FieldVector<Scalar,numEq> initial(const FieldVector<Scalar,dim>& x, const Element& e,
                                              const FieldVector<Scalar,dim>& xi) const
    {
        DUNE_THROW(NotImplemented, "no q specified, but requested");

        FieldVector<Scalar,numEq> result(0);
        return result;
    }

    virtual FieldVector<Scalar,dim+1> g(const FieldVector<Scalar,dim>& x, const Element& e,
                                        const IntersectionIterator& intersectionIt,
                                        const FieldVector<Scalar,dim>& xi) const = 0;

    virtual FieldVector<Scalar,dim+1> J(const FieldVector<Scalar,dim>& x, const Element& e,
                                        const IntersectionIterator& intersectionIt,
                                        const FieldVector<Scalar,dim>& xi)
    {
        DUNE_THROW(NotImplemented, "no J specified, but requested");

        FieldVector<Scalar,dim+1> result(0);
        return result;
    }

    virtual FieldMatrix<Scalar,dim,dim> D (const FieldVector<Scalar,dim>& x, const Element& e,
                                           const FieldVector<Scalar,dim>& xi) const
    {
        DUNE_THROW(NotImplemented, "Diffusivity ?!");
        FieldMatrix<Scalar,dim,dim> result(0);
        return result;
    }

    virtual FieldVector<Scalar,dim> velocity(const FieldVector<Scalar,dim>& x) const
    {
        DUNE_THROW(NotImplemented, "no exact solution available");

        FieldVector<Scalar,dim> result(0);
        return result;
    }

    virtual FieldVector<Scalar,numEq> velocitypressure(const FieldVector<Scalar,dim>& x) const
    {
        DUNE_THROW(NotImplemented, "no exact solution available");

        FieldVector<Scalar,numEq> result(0);
        return result;
    }

    virtual Scalar pressure(const FieldVector<Scalar,dim>& x) const
    {
        DUNE_THROW(NotImplemented, "no exact solution available");

        return 0;
    }

  virtual Scalar Qg(const FieldVector<Scalar,dim>& x, const Element& e,
				    const FieldVector<Scalar,dim>& xi) const
   {
     DUNE_THROW(NotImplemented, "no exact solution available");

     return 0;
   }


    virtual Scalar partialdensity(const FieldVector<Scalar,dim>& x) const
    {
        DUNE_THROW(NotImplemented, "no exact solution available");

        return 0;
    }

  virtual Scalar density(const FieldVector<Scalar,dim>& x, const Element& e,
          const FieldVector<Scalar,dim>& xi) const
    {
      DUNE_THROW(NotImplemented, "no exact solution available");

      return 0;
    }

  virtual Scalar gravity() const= 0;

    virtual FieldVector<Scalar,numEq> velocitypressuremassfrac(const FieldVector<Scalar,dim>& x) const
    {
        DUNE_THROW(NotImplemented, "no exact solution available");
        FieldVector<Scalar,numEq> result(0);

        return result;
    }

    virtual FieldVector<Scalar,dim+1> velocitymassfrac(const FieldVector<Scalar,dim>& x) const
    {
        DUNE_THROW(NotImplemented, "no exact solution available");
        FieldVector<Scalar,dim+1> result(0);

        return result;
    }
    /*
      virtual FieldMatrix<Scalar, dim, dim> velocityGradient(const FieldVector<Scalar,dim>& x) const
      {
      DUNE_THROW(NotImplemented, "no exact solution available");

      FieldMatrix<Scalar, dim, dim> result(0);
      return result;
      }
    */
    virtual void dirichletIndex(const Dune::FieldVector<Scalar,dim>& x, const Element& e,
                                const IntersectionIterator& intersectionIt,
                                const Dune::FieldVector<Scalar,dim>& xi, Dune::FieldVector<int,numEq>& dirichletIndex) const
    {
        for (int i = 0; i < numEq; i++)
            dirichletIndex[i]=i;
        return;
    }

    virtual Fluid& gasPhase () const
    {
        return gasPhase_;
    }

    virtual MultiComp& multicomp () const
    {
        return multicomp_;
    }

    StokesTransportProblem(Fluid& gasPhase, MultiComp& multicomp = *(new CWaterAir))
        : gasPhase_(gasPhase),
          multicomp_(multicomp)
    {}

    //! always define virtual destructor in abstract base class
    virtual ~StokesTransportProblem()
    {}

protected:
    Fluid& gasPhase_;
    MultiComp& multicomp_;
};

}
#endif
