#ifndef CAPILLARYDARCYPROBLEM_HH
#define CAPILLARYDARCYPROBLEM_HH

#include<dumux/coupled/coupledporousmediaproblem.hh>
#include<dune/common/exceptions.hh>


namespace Dune
{

template<class Grid, class Scalar>
class CapillaryDarcyProblem  : public CoupledPorousMediaProblem<Grid, Scalar>
{
    enum {dim=Grid::dimension, numEq=1};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;

public:
    CapillaryDarcyProblem ()
    {
        eps_ = 1e-10;
        // CHANGE also in the Stokes problem!
        for (int i=0; i<dim; i++)
            for (int j=0; j<dim; j++)
                if (i==j)
                {
                    permTissue_[i][j] = 5e-15; // provide intrinsic K / mu
                    permSmall_[i][j] = 1e-20;
                }
                else
                {
                    permTissue_[i][j] = permSmall_[i][j] = 0;
                }
    }

    const FieldMatrix<Scalar,dim,dim>& K (const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                          const FieldVector<Scalar,dim>& localPos) const
    {
        if (globalPos[0] > 0.2e-3 && globalPos[0] < 0.8e-3)
            return permTissue_;
        else
            return permSmall_;
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
        if (globalPos[1] < -3e-5 + eps_)// || globalPos[0] < 1e-10 || globalPos[0] > 1e-3 - 1e-10)
            return BoundaryConditions::dirichlet;

        return BoundaryConditions::neumann;
    }

    // Dirichlet boundary conditions
    Scalar g (const FieldVector<Scalar,dim>& globalPos, const Element& element,
              const IntersectionIterator& intersectionIt,
              const FieldVector<Scalar,dim>& localPos) const
    {
        if (globalPos[0] < 1e-10)
            return ((-1330.0 - 266.0)/(-3e-5)*globalPos[1] + 266.0);
        else if (globalPos[0] > 1e-3 - 1e-10)
            return ((-1330.0 + 1729.0)/(-3e-5)*globalPos[1] - 1729.0);
        else
            return -1330.0;
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
        return (0);
    }

    FieldVector<Scalar,dim> exactGrad(const FieldVector<Scalar,dim>& globalPos) const
    {
        FieldVector<Scalar,dim> grad(0);

        return grad;
    }

private:
    FieldMatrix<Scalar,dim,dim> permTissue_;
    FieldMatrix<Scalar,dim,dim> permSmall_;
    Scalar eps_;
};
} // end namespace
#endif
