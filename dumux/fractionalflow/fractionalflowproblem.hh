// $Id$

#ifndef DUNE_FRACTIONALFLOWPROBLEM_HH
#define DUNE_FRACTIONALFLOWPROBLEM_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dumux/material/twophaserelations.hh>

#include <dumux/material/property_baseclasses.hh>
//#include "dumux/fractionalflow/variableclass.hh"

/**
 * @file
 * @brief  Base class for defining an instance of the diffusion problem
 * @author Bernd Flemisch
 */
/**
 * \ingroup diffusion
 * \defgroup diffusionProblems Problems
 */

namespace Dune
{
/*! \ingroup diffusionProblems
 * @brief base class that defines the parameters of a diffusion equation
 *
 * An interface for defining parameters for the stationary diffusion equation
 * \f$ - \text{div}\, (\lambda K \text{grad}\, p ) = q, \f$,
 * \f$p = g\f$ on \f$\Gamma_1\f$, and \f$\lambda K \text{grad}\, p = J\f$
 * on \f$\Gamma_2\f$. Here,
 * \f$p\f$ denotes the pressure, \f$K\f$ the absolute permeability,
 * and \f$\lambda\f$ the total mobility, possibly depending on the
 * saturation.
 *
 *    Template parameters are:
 *
 *    - Grid  a DUNE grid type
 *    - Scalar    type used for return values
 */
template<class Grid, class Scalar, class VC>
class FractionalFlowProblem
{
public:
    enum
    {    dim=Grid::dimension,dimWorld= Grid::dimensionworld,numEq=1};
protected:
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;

public:
    //! evaluate source term for the pressure equation
    /*! evaluate source term for the pressure equation at given location
     @param[in]  globalPos    position in global coordinates
     @param[in]  element    entity of codim 0
     @param[in]  localPos   position in reference element of element
     \return     value of source term
     */
    virtual Scalar sourcePress (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) = 0;

    //! return type of boundary condition for the pressure equation at the given global coordinate
    /*! return type of boundary condition for the pressure equation at the given global coordinate
     @param[in]  globalPos    position in global coordinates
     \return     boundary condition type given by enum in this class
     */
    virtual BoundaryConditions::Flags bctypePress (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const = 0;

    //! return type of boundary condition for the saturation equation at the given global coordinate
    /*! return type of boundary condition for the saturation equation at the given global coordinate
     @param[in]  globalPos    position in global coordinates
     \return     boundary condition type given by enum in this class
     */
    virtual BoundaryConditions::Flags bctypeSat (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const = 0;

    //! evaluate Dirichlet boundary condition for the pressure equation at given position
    /*! evaluate Dirichlet boundary condition for the pressure equation at given position
     @param[in]  globalPos    position in global coordinates
     \return     boundary condition value
     */
    virtual Scalar dirichletPress (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const = 0;

    //! evaluate Dirichlet boundary condition for the saturation equation at given position
    /*! evaluate Dirichlet boundary condition for the saturation equation at given position
     @param[in]  globalPos    position in global coordinates
     \return     boundary condition value
     */
    virtual Scalar dirichletSat (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        return 1;
    }

    //! evaluate Neumann boundary condition for the pressure equation at given position
    /*! evaluate Neumann boundary condition for the pressure equation at given position
     @param[in]  globalPos    position in global coordinates
     \return     boundary condition value
     */
    virtual Scalar neumannPress (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const = 0;

    //! evaluate initial condition for saturation at given position
    /*! evaluate initial condition for saturation at given position
     @param[in]  globalPos    position in global coordinates
     \return    initial condition value
     */
    virtual Scalar initSat (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const = 0;

    virtual Scalar neumannSat (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos, Scalar factor) const
    {
        return 0;
    }

    //! evaluate gravity
    /*! evaluate gravity
     \return     gravity vector
     */
    const FieldVector<Scalar,dimWorld>& gravity() const
    {
        return gravity_;
    }

    //! constructor
    /** @param law implementation of Material laws. Class TwoPhaseRelations or derived.
     *  @param cap flag to include capillary forces.
     */
    FractionalFlowProblem(VC& variables, Fluid& wettingphase, Fluid& nonwettingphase, Matrix2p<Grid, Scalar>& soil, TwoPhaseRelations<Grid, Scalar>& materialLaw = *(new TwoPhaseRelations<Grid,Scalar>), const bool capillarity = false)
    : variables(variables), wettingPhase(wettingphase), nonWettingPhase(nonwettingphase), soil(soil), capillarity(capillarity), materialLaw(materialLaw),gravity_(0)
    {}

    //! always define virtual destructor in abstract base class
    virtual ~FractionalFlowProblem ()
    {}

    //! a class describing relations between two phases and the porous medium
    VC& variables;
    Fluid& wettingPhase;
    Fluid& nonWettingPhase;
    Matrix2p<Grid, Scalar>& soil;
    const bool capillarity;
    TwoPhaseRelations<Grid, Scalar>& materialLaw;
private:
    FieldVector<Scalar,dimWorld> gravity_;

};

}
#endif
