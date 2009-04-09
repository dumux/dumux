#ifndef DARCYPROBLEM_HH
#define DARCYPROBLEM_HH

#include<dumux/coupled/coupledporousmediaproblem.hh>
#include<dune/common/exceptions.hh>


namespace Dune
{

template<class Grid, class Scalar>
class DarcyProblem  : public CoupledPorousMediaProblem<Grid, Scalar>
{
    enum {dim=Grid::dimension, numEq=1};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<0>::LeafIntersectionIterator IntersectionIterator;

public:
    DarcyProblem ()
    {
        eps_ = 1e-6;
        // CHANGE also in the Stokes problem!
        for (int i=0; i<dim; i++)
            for (int j=0; j<dim; j++)
                if (i==j)
                    permeability_[i][j] = 1e0; // provide intrinsic K / mu
                else
                    permeability_[i][j] = 0;
    }

    const FieldMatrix<Scalar,dim,dim>& K (const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                          const FieldVector<Scalar,dim>& localPos) const
    {
        return permeability_;
    }

    Scalar q (const FieldVector<Scalar,dim>& globalPos, const Element& element,
              const FieldVector<Scalar,dim>& localPos) const
    {
        return 0;
    }

    BoundaryConditions::Flags bctype (const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                      const IntersectionIterator& intersectionIt,
                                      const FieldVector<Scalar,dim>& localPos) const
    {
        if (globalPos[0] > 4 - eps_)
            return BoundaryConditions::dirichlet;

        return BoundaryConditions::neumann;
    }

    // Dirichlet boundary conditions
    Scalar g (const FieldVector<Scalar,dim>& globalPos, const Element& element,
              const IntersectionIterator& intersectionIt,
              const FieldVector<Scalar,dim>& localPos) const
    {
        return 0;
    }

    // Neumann boundary conditions
    Scalar J (const FieldVector<Scalar,dim>& globalPos, const Element& element,
	      const IntersectionIterator& intersectionIt,
	      const FieldVector<Scalar,dim>& localPos) const
    {
        return 0;
    }

    Scalar exact(const FieldVector<Scalar,dim>& globalPos) const
    {
        return (globalPos[0]*globalPos[1]);
    }

    FieldVector<Scalar,dim> exactGrad(const FieldVector<Scalar,dim>& globalPos) const
    {
        FieldVector<Scalar,dim> grad;

        grad[0] = globalPos[1];
        grad[1] = globalPos[0];

        return grad;
    }

private:
    FieldMatrix<Scalar,dim,dim> permeability_;
    Scalar eps_;
};
} // end namespace
#endif
