// $Id$

#ifndef DUNE_DIFFUSIONPROBLEM_HH
#define DUNE_DIFFUSIONPROBLEM_HH

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
 *    - RT    type used for return values
 */
template<class GridView, class Scalar, class VC>
class DiffusionProblem
{

protected:
    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
typedef    typename GridView::Grid Grid;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

public:

    //! evaluate source term
    /*! evaluate source term at given location
     @param[in]  globalPos    position in global coordinates
     @param[in]  element    entity of codim 0
     @param[in]  localPos   position in reference element of element
     \return     value of source term
     */
    virtual Scalar sourcePress (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) = 0;

    //! return type of boundary condition at the given global coordinate
    /*! return type of boundary condition at the given global coordinate
     @param[in]  globalPos    position in global coordinates
     \return     boundary condition type given by enum in this class
     */
    virtual BoundaryConditions::Flags bctypePress (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const = 0;

    //! return type of boundary condition at the given global coordinate
    /*! return type of boundary condition at the given global coordinate
     @param[in]  globalPos    position in global coordinates
     \return     boundary condition type given by enum in this class
     */
    virtual BoundaryConditions::Flags bctypeSat (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        return BoundaryConditions::dirichlet;
    }

    //! evaluate Dirichlet boundary condition at given position
    /*! evaluate Dirichlet boundary condition at given position
     @param[in]  globalPos    position in global coordinates
     \return     boundary condition value
     */
    virtual Scalar dirichletPress (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const = 0;

    //! evaluate Dirichlet boundary condition at given position
    /*! evaluate Dirichlet boundary condition at given position
     @param[in]  globalPos    position in global coordinates
     \return     boundary condition value
     */
    virtual Scalar dirichletSat (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        return 1;
    }

    //! evaluate Neumann boundary condition at given position
    /*! evaluate Neumann boundary condition at given position
     @param[in]  globalPos    position in global coordinates
     \return     boundary condition value
     */
    virtual Scalar neumannPress (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const = 0;

    virtual const FieldVector<Scalar,dim>& gravity() const
    {
        return gravity_;
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
    virtual VC& variables ()
    {
        return variables_;
    }
    virtual Fluid& wettingPhase () const
    {
        return wettingPhase_;
    }
    virtual Fluid& nonWettingPhase () const
    {
        return nonWettingPhase_;
    }

    //! constructor
    /** @param law implementation of Material laws. Class TwoPhaseRelations or derived.
     *  @param cap flag to include capillary forces.
     */
    DiffusionProblem(VC& variables,Fluid& wettingPhase, Fluid& nonWettingPhase, Matrix2p<Grid, Scalar>& soil, TwoPhaseRelations<Grid, Scalar>& materialLaw = *(new TwoPhaseRelations<Grid,Scalar>))
    : variables_(variables), materialLaw_(materialLaw),
    wettingPhase_(wettingPhase), nonWettingPhase_(nonWettingPhase), soil_(soil),
    gravity_(0)
    {}

    //! constructor
    /** @param law implementation of Material laws. Class TwoPhaseRelations or derived.
     *  @param cap flag to include capillary forces.
     */
    DiffusionProblem(VC& variables, TwoPhaseRelations<Grid, Scalar>& materialLaw = *(new TwoPhaseRelations<Grid,Scalar>))
    : variables_(variables), materialLaw_(materialLaw),
    wettingPhase_(materialLaw_.wettingPhase), nonWettingPhase_(materialLaw.nonwettingPhase), soil_(materialLaw.soil), gravity_(0)
    {}

    //! always define virtual destructor in abstract base class
    virtual ~DiffusionProblem ()
    {}

    //! a class describing relations between two phases and the porous medium
private:
    VC& variables_;
    TwoPhaseRelations<Grid,Scalar>& materialLaw_;
    Fluid& wettingPhase_;
    Fluid& nonWettingPhase_;
    Matrix2p<Grid, Scalar>& soil_;
    FieldVector<Scalar,dim> gravity_;
};

}
#endif
