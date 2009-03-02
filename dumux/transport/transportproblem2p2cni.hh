//$Id$
#ifndef TRANSPORTPROBLEM2P2CNI_HH
#define TRANSPORTPROBLEM2P2CNI_HH

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
#include <dumux/fractionalflow/variableclass2p2cni.hh>

namespace Dune
{
//! Base class for the definition of nonisothermal 2p2c problems
/** This base class defines all boundary and initial functions which are needed
 * for a decoupled nonisothermal 2p2c computation.
 */
template<class G, class RT>
class TransportProblem2p2cni
{

    typedef typename G::ctype DT;
    enum {n=G::dimension, m=1, blocksize=2*G::dimension};
    typedef typename G::Traits::template Codim<0>::Entity Entity;

public:
    //! Type of concentration boundary condition.
    /**    either the concentration or the saturation have to be defined
     * on boundaries with dirichlet pressure BCs.
     * @param x global coordinates
     * @param e reference to the cell for which the function is to be evaluated
     * @param xi local coordinates inside e
     */
    virtual BoundaryConditions2p2c::Flags cbctype (const FieldVector<DT,n>& x, const Entity& e,
                                                   const FieldVector<DT,n>& xi) const = 0;

    //! Type of concentration initisl condition.
    /**    either the concentration or the saturation have to be defined
     * as initial condition.
     * @param x global coordinates
     * @param e reference to the cell for which the function is to be evaluated
     * @param xi local coordinates inside e
     */
    virtual BoundaryConditions2p2c::Flags ictype (const FieldVector<DT,n>& x, const Entity& e,
                                                  const FieldVector<DT,n>& xi) const = 0;

    //! Type of pressure boundary condition.
    /**    Pressure (dirichlet) or flux (neumann) have to be defined on boundaries.
     * @param x global coordinates
     * @param e reference to the cell for which the function is to be evaluated
     * @param xi local coordinates inside e
     */
    virtual BoundaryConditions::Flags pbctype (const FieldVector<DT,n>& x, const Entity& e,
                                               const FieldVector<DT,n>& xi) const = 0;

    //! Feed concentration boundary condition
    /** Feed concentration is the (global) mass fraction of component 1 in the mixture
     * @param x global coordinates
     * @param e reference to the cell for which the function is to be evaluated
     * @param xi local coordinates inside e
     */
    virtual RT gZ (const FieldVector<DT,n>& x, const Entity& e,
                   const FieldVector<DT,n>& xi) const = 0;

    //! Saturation boundary condition
    /** @param x global coordinates
     * @param e reference to the cell for which the function is to be evaluated
     * @param xi local coordinates inside e
     */
    virtual RT gS (const FieldVector<DT,n>& x, const Entity& e,
                   const FieldVector<DT,n>& xi) const = 0;

    //! Pressure (dirichlet) boundary condition
    /** @param x global coordinates
     * @param e reference to the cell for which the function is to be evaluated
     * @param xi local coordinates inside e
     */
    virtual RT gPress (const FieldVector<DT,n>& x, const Entity& e,
                       const FieldVector<DT,n>& xi) const = 0;

    //! Temperature boundary condition
    /** @param x global coordinates
     * @param e reference to the cell for which the function is to be evaluated
     * @param xi local coordinates inside e
     */
    virtual RT gT (const FieldVector<DT,n>& x, const Entity& e,
                   const FieldVector<DT,n>& xi) const = 0;

    //! Flux (neumann) boundary condition
    /** @param x global coordinates
     * @param e reference to the cell for which the function is to be evaluated
     * @param xi local coordinates inside e
     */
    virtual FieldVector<RT,3> J (const FieldVector<DT,n>& x, const Entity& e,
                                 const FieldVector<DT,n>& xi) const = 0;

    //! Source of components
    /** Describes the source of the components per unit volume in \f$ \left[ frac{kg}{m^3} \right] \f$
     * @param x global coordinates
     * @param e reference to the cell for which the function is to be evaluated
     * @param xi local coordinates inside e
     */
    virtual FieldVector<RT,2> q (const FieldVector<DT,n>& x, const Entity& e,
                                 const FieldVector<DT,n>& xi) const = 0;

    //! Heat source
    /** Heat flux from outside per unit volume in \f$ \left[ frac{W}{m^3} \right] \f$
     * @param x global coordinates
     * @param e reference to the cell for which the function is to be evaluated
     * @param xi local coordinates inside e
     */
    virtual RT qh (const FieldVector<DT,n>& x, const Entity& e,
                   const FieldVector<DT,n>& xi) const = 0;

    //! Saturation initial condition
    /** @param x global coordinates
     * @param e reference to the cell for which the function is to be evaluated
     * @param xi local coordinates inside e
     */
    virtual RT S0 (const FieldVector<DT,n>& x, const Entity& e,
                   const FieldVector<DT,n>& xi) const = 0;

    //! Feed concentration initial condition
    /** @param x global coordinates
     * @param e reference to the cell for which the function is to be evaluated
     * @param xi local coordinates inside e
     */
    virtual RT Z1_0 (const FieldVector<DT,n>& x, const Entity& e,
                     const FieldVector<DT,n>& xi) const = 0;

    //! Temperature initial condition \f$ \left[ K \right \f$
    /** @param x global coordinates
     * @param e reference to the cell for which the function is to be evaluated
     * @param xi local coordinates inside e
     */
    virtual RT T_0 (const FieldVector<DT,n>& x, const Entity& e,
                    const FieldVector<DT,n>& xi) const = 0;

    //! gravity vector
    const FieldVector<DT,n>& gravity()
    {
        return gravity_;
    }


    //! Constructor
    /**
     *
     */
    TransportProblem2p2cni(Dune::VariableClass2p2cni<G, RT>& var, Liquid_GL& liq, Gas_GL& gas, Matrix2p<G, RT>& s, TwoPhaseRelations<G, RT>& law,
                           const bool cap = false)
        :variables(var), liquidPhase(liq), gasPhase(gas), soil(s), capillary(cap), materialLaw(law)
    {
    }

    virtual ~TransportProblem2p2cni ()
    {
    }

    FieldVector<DT,n> gravity_;
    const bool capillary;
    TwoPhaseRelations<G, RT>& materialLaw;
    Liquid_GL& liquidPhase;
    Gas_GL& gasPhase;
    Matrix2p<G, RT>& soil;
    VariableClass2p2cni<G, RT>& variables;
};

}
#endif
