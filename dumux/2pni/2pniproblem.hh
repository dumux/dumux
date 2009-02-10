// $Id$

#ifndef DUNE_TWOPHASEHEATPROBLEM_HH
#define DUNE_TWOPHASEHEATPROBLEM_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/grid/utility/intersectiongetter.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dumux/material/twophaserelations.hh>
#include<dumux/material/property_baseclasses.hh>
#include<dumux/auxiliary/basicdomain.hh>

/**
 * @file
 * @brief  Base class for defining an instance of the non-isothermal two-phase problem
 * @author Bernd Flemisch, Melanie Darcis
 */

namespace Dune {

/** \todo Please doc me! */

template<class Grid, class Scalar> class TwoPhaseHeatProblem {
    typedef typename Grid::ctype CoordScalar;
    enum {dim=Grid::dimension, numEq=3};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator
            IntersectionIterator;
    typedef BasicDomain<Grid, Scalar> ParentType;
    // the domain traits of the domain
    typedef typename ParentType::DomainTraits DomainTraits;
    typedef typename DomainTraits::LocalPosition LocalPosition;
    typedef typename DomainTraits::GlobalPosition GlobalPosition;

public:
    //! evaluate source term
    /*! evaluate source term at given location
     @param[in]  x    position in global coordinates
     @param[in]  element    entity of codim 0
     @param[in]  xi   position in reference element of e
     \return     value of source term
     */
    virtual FieldVector<Scalar,numEq> q(const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const = 0;

    //! return type of boundary condition at the given global coordinate
    /*! return type of boundary condition at the given global coordinate
     @param[in]  x    position in global coordinates
     \return     boundary condition type given by enum in this class
     */
    //    virtual FieldVector<BoundaryConditions::Flags, numEq> bctype (const GlobalPosition& globalPos, const Element& element,
    //            const IntersectionIterator& isIt,
    //            const LocalPosition& localPos) const = 0;

    virtual FieldVector<BoundaryConditions::Flags, numEq>bctype(
            const GlobalPosition& globalPos, const Element& element,
            const IntersectionIterator& isIt,
            const LocalPosition& localPos) const = 0;

    //! returns index of the primary variable corresponding to the dirichlet boundary condition at the given global coordinate
        /*! returns index of the primary variable corresponding to the dirichlet boundary condition at the given global coordinate
         @param[in]  x    position in global coordinates
         \return     index of the primary variable
         */

    virtual void dirichletIndex(const GlobalPosition& globalPos, const Element& element,
            const IntersectionIterator& isIt,
            const LocalPosition& localPos, FieldVector<int,numEq>& dirichletIdx) const
    {
        for (int equationnumber = 0; equationnumber < numEq; equationnumber++)
            dirichletIdx[equationnumber]=equationnumber;
        return;
    }

    //! evaluate Dirichlet boundary condition at given position
    /*! evaluate Dirichlet boundary condition at given position
     @param[in]  x    position in global coordinates
     \return     boundary condition value
     */
    virtual FieldVector<Scalar,numEq> g(const GlobalPosition& globalPos, const Element& element,
            const IntersectionIterator& isIt,
            const LocalPosition& localPos) const = 0;

    //! evaluate Neumann boundary condition at given position
    /*! evaluate Neumann boundary condition at given position
     @param[in]  x    position in global coordinates
     \return     boundary condition value
     */
    virtual FieldVector<Scalar,numEq> J(const GlobalPosition& globalPos, const Element& element,
            const IntersectionIterator& isIt,
            const LocalPosition& localPos) const = 0;

    //! evaluate initial condition at given position
    /*! evaluate initial boundary condition at given position
     @param[in]  x    position in global coordinates
     \return     boundary condition value
     */
    virtual FieldVector<Scalar,numEq> initial(const GlobalPosition& globalPos,
            const Element& element, const LocalPosition& localPos) const = 0;

    virtual FieldVector<Scalar,dim> gravity() const = 0;

    //! properties of the wetting (liquid) phase
    /*! properties of the wetting (liquid) phase
      \return    wetting phase
     */
    virtual Fluid& wettingPhase () const
    {
        return wettingPhase_;
    }

    //! properties of the nonwetting (liquid) phase
    /*! properties of the nonwetting (liquid) phase
      \return    nonwetting phase
     */
    virtual Fluid& nonwettingPhase () const
    {
        return nonwettingPhase_;
    }

    //! properties of the soil
    /*! properties of the soil
      \return    soil
     */
    virtual Matrix2p<Grid, Scalar>& soil () const
    {
        return soil_;
    }

    //! object for definition of material law
    /*! object for definition of material law (e.g. Brooks-Corey, Van Genuchten, ...)
      \return    material law
     */
    virtual TwoPhaseRelations<Grid, Scalar>& materialLaw () const
    {
        return materialLaw_;
    }

    TwoPhaseHeatProblem(Fluid& liq1, Fluid& liq2, Matrix2p<Grid, Scalar>& soil,
            TwoPhaseRelations<Grid,Scalar>& materialLaw = *(new TwoPhaseRelations<Grid,Scalar>))
    : wettingPhase_(liq1), nonwettingPhase_(liq2), soil_(soil),
      materialLaw_(materialLaw)
      {     }

    //! always define virtual destructor in abstract base class
    virtual ~TwoPhaseHeatProblem() {
    }

protected:
    Fluid& wettingPhase_;
    Fluid& nonwettingPhase_;
    Matrix2p<Grid, Scalar>& soil_;
    TwoPhaseRelations<Grid, Scalar>& materialLaw_;
};

}
#endif
