// $Id$

#ifndef TRANSPORTPROBLEM2P2C_HH
#define TRANSPORTPROBLEM2P2C_HH

#include <iostream>
#include <iomanip>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/referenceelements.hh>
#include <dune/disc/operators/boundaryconditions.hh>

#include <dumux/material/property_baseclasses.hh>
#include <dumux/material/twophaserelations.hh>
#include <dumux/operators/boundaryconditions2p2c.hh>
#include <dumux/fractionalflow/variableclass2p2c.hh>



namespace Dune
{
//! Base class for the definition of 2p2c problems
/** This base class defines all boundary and initial functions which are needed
 * for a decoupled 2p2c computation.
 */
template<class Grid, class Scalar>
class TransportProblem2p2c
{
    enum {dim=Grid::dimension};
    typedef typename Grid::Traits::template Codim<0>::Entity Entity;

public:
    //! Type of concentration boundary condition.
    /**    either the concentration or the saturation have to be defined
     * on boundaries with dirichlet pressure BCs.
     * @param globalPos global coordinates
     * @param element reference to the cell for which the function is to be evaluated
     * @param localPos local coordinates inside element
     */
    virtual BoundaryConditions2p2c::Flags bc_type (const FieldVector<Scalar,dim>& globalPos, const Entity& element,
                                                   const FieldVector<Scalar,dim>& localPos) const = 0;

    //! Type of concentration initisl condition.
    /**    either the concentration or the saturation have to be defined
     * as initial condition.
     * @param globalPos global coordinates
     * @param element reference to the cell for which the function is to be evaluated
     * @param localPos local coordinates inside element
     */
    virtual BoundaryConditions2p2c::Flags initcond_type (const FieldVector<Scalar,dim>& globalPos, const Entity& element,
                                                         const FieldVector<Scalar,dim>& localPos) const = 0;

    //! Type of pressure boundary condition.
    /**    Pressure (dirichlet) or flux (neumann) have to be defined on boundaries.
     * @param globalPos global coordinates
     * @param element reference to the cell for which the function is to be evaluated
     * @param localPos local coordinates inside element
     */
    virtual BoundaryConditions::Flags press_bc_type (const FieldVector<Scalar,dim>& globalPos, const Entity& element,
                                                     const FieldVector<Scalar,dim>& localPos) const = 0;

    //! Permeability tensor \f$ [m^2] \f$
    /**
     * @param globalPos global coordinates
     * @param element reference to the cell for which the function is to be evaluated
     * @param localPos local coordinates inside element
     */
    //    virtual const FieldMatrix<DT,n,n>& K (const FieldVector<DT,n>& globalPos, const Entity& element, const FieldVector<DT,n>& localPos) const
    //    {
    //        return soil.K(globalPos, element, localPos);
    //    }

    //! Feed concentration boundary condition
    /** Feed concentration is the (global) mass fraction of component 1 in the mixture
     * @param globalPos global coordinates
     * @param element reference to the cell for which the function is to be evaluated
     * @param localPos local coordinates inside element
     */
    virtual Scalar dirichletConcentration (const FieldVector<Scalar,dim>& globalPos, const Entity& element,
                                           const FieldVector<Scalar,dim>& localPos) const = 0;

    //! Saturation boundary condition
    /** @param globalPos global coordinates
     * @param element reference to the cell for which the function is to be evaluated
     * @param localPos local coordinates inside element
     */
    virtual Scalar dirichletSat (const FieldVector<Scalar,dim>& globalPos, const Entity& element,
                                 const FieldVector<Scalar,dim>& localPos) const = 0;

    //! Pressure (dirichlet) boundary condition
    /** @param globalPos global coordinates
     * @param element reference to the cell for which the function is to be evaluated
     * @param localPos local coordinates inside element
     */
    virtual Scalar dirichlet (const FieldVector<Scalar,dim>& globalPos, const Entity& element,
                              const FieldVector<Scalar,dim>& localPos) const = 0;

    //! Flux (neumann) boundary condition
    /** @param globalPos global coordinates
     * @param element reference to the cell for which the function is to be evaluated
     * @param localPos local coordinates inside element
     */
    virtual FieldVector<Scalar,2> neumann (const FieldVector<Scalar,dim>& globalPos, const Entity& element,
                                           const FieldVector<Scalar,dim>& localPos) const = 0;

    //! Source of components
    /** Describes the source of the components per unit area
     * @param globalPos global coordinates
     * @param element reference to the cell for which the function is to be evaluated
     * @param localPos local coordinates inside element
     */
    virtual FieldVector<Scalar,2> source (const FieldVector<Scalar,dim>& globalPos, const Entity& element,
                                          const FieldVector<Scalar,dim>& localPos) const = 0;

    //! Saturation initial condition
    /** @param globalPos global coordinates
     * @param element reference to the cell for which the function is to be evaluated
     * @param localPos local coordinates inside element
     */
    virtual Scalar initSat (const FieldVector<Scalar,dim>& globalPos, const Entity& element,
                            const FieldVector<Scalar,dim>& localPos) const = 0;

    //! Feed concentration initial condition
    /** @param globalPos global coordinates
     * @param element reference to the cell for which the function is to be evaluated
     * @param localPos local coordinates inside element
     */
    virtual Scalar initConcentration (const FieldVector<Scalar,dim>& globalPos, const Entity& element,
                                      const FieldVector<Scalar,dim>& localPos) const = 0;

    //! gravity vector
    virtual const FieldVector<Scalar,dim> gravity()
    {
        FieldVector<Scalar,dim> gravity_(0.);
        return gravity_;
    }


    //! Constructor
    /**
     *
     */

    TransportProblem2p2c(Dune::VariableClass2p2c<Grid, Scalar>& var, Liquid_GL& liq, Gas_GL& gas, Matrix2p<Grid, Scalar>& s,
                         TwoPhaseRelations<Grid, Scalar>& law = *(new TwoPhaseRelations<Grid,Scalar>),const bool cap = false)
        :variables(var), liquidPhase(liq), gasPhase(gas), soil(s), materialLaw(law), capillary(cap)
    {
    }

    virtual ~TransportProblem2p2c ()
    {
    }

    VariableClass2p2c<Grid, Scalar>& variables;
    Liquid_GL& liquidPhase;
    Gas_GL& gasPhase;
    Matrix2p<Grid, Scalar>& soil;
    TwoPhaseRelations<Grid, Scalar>& materialLaw;
    const bool capillary;
};

}
#endif
