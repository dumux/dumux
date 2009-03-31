#ifndef COUPLEDPOROUSMEDIAPROBLEM_HH
#define COUPLEDPOROUSMEDIAPROBLEM_HH

#include<dune/common/exceptions.hh>

#include<dune/disc/operators/boundaryconditions.hh>

namespace Dune {

template<class Grid, class Scalar>
class CoupledPorousMediaProblem
{
    enum {dim=Grid::dimension, numEq=1};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<0>::LeafIntersectionIterator IntersectionIterator;

public:
    //TODO Permeability can be defined here. This is only preliminary.
    // Should be defined in the soil class
    virtual const FieldMatrix<Scalar,dim,dim>& K (const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                                  const FieldVector<Scalar,dim>& localPos) const = 0;

    //! evaluate source term of the momentum equation
    /*! evaluate source term of the momentum equation at given location
      @param[in]  globalPos    position in global coordinates
      @param[in]  element    entity of codim 0
      @param[in]  localPos   position in reference element
      \return     value of source term
    */
    virtual Scalar q   (const FieldVector<Scalar,dim>& globalPos, const Element& element,
                        const FieldVector<Scalar,dim>& localPos) const
    {
        DUNE_THROW(NotImplemented, "no q specified, but requested");

        FieldVector<Scalar,numEq> result(0);
        return result;
    }

    //! return type of boundary condition at the given global coordinate
    /*! return type of boundary condition at the given global coordinate
      @param[in]  globalPos    position in global coordinates
      \return     boundary condition type given by enum in this class
    */
    virtual BoundaryConditions::Flags bctype (const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                              const IntersectionIterator& intersectionIt,
                                              const FieldVector<Scalar,dim>& localPos) const = 0;

    //! evaluate velocity Dirichlet boundary condition at given position
    /*! evaluate velocity Dirichlet boundary condition at given position
      @param[in]  globalPos    position in global coordinates
      \return     boundary condition value
    */
    virtual Scalar g (const FieldVector<Scalar,dim>& globalPos, const Element& element,
                      const IntersectionIterator& intersectionIt,
                      const FieldVector<Scalar,dim>& localPos) const = 0;

    //! evaluate Neumann boundary condition at given position
    /*! evaluate Neumann boundary condition at given position
      @param[in]  globalPos    position in global coordinates
      \return     boundary condition value
    */
    virtual Scalar J (const FieldVector<Scalar,dim>& globalPos, const Element& element,
                      const IntersectionIterator& intersectionIt,
                      const FieldVector<Scalar,dim>& localPos) const = 0;

    virtual void dirichletIndex(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                const IntersectionIterator& intersectionIt,
                                const FieldVector<Scalar,dim>& localPos, FieldVector<int,numEq>& dirichletIndex) const
    {
        for (int i = 0; i < numEq; i++)
            dirichletIndex[i]=i;
        return;
    }

    virtual Scalar exact(const FieldVector<Scalar,dim>& globalPos) const
    {
        DUNE_THROW(NotImplemented, "no exact solution specified, but requested");
    }

    virtual FieldVector<Scalar,dim> exactGrad(const FieldVector<Scalar,dim>& globalPos) const
    {
        DUNE_THROW(NotImplemented, "no exact gradient specified, but requested");
    }

    CoupledPorousMediaProblem ()
    { }

    virtual ~CoupledPorousMediaProblem()
    { }

    //private:
    //    FieldMatrix<Scalar,dim,dim> permeability_;
};
} // end namespace
#endif
