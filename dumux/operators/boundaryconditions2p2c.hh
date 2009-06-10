// $Id$

#ifndef DUNE_BOUNDARYCONDITIONS2P2C_HH
#define DUNE_BOUNDARYCONDITIONS2P2C_HH


namespace Dune
{
/** \ingroup decoupled2p2c
 * \brief Defines type of boundary conditions for 2p2c processes
 *
 *  This is to distinguish BC types for 2p2c processes similar to
 * the class Dune::BoundaryConditions which distinguishes between
 * dirichlet, process and neumann.
 * BoundaryConditions for compositional transport can either be
 * defined as saturations or as total concentrations (decoupled method).
 * Each leads via pressure and constitutive relationships to the other
 * and to the phase concentrations.
 */
struct BoundaryConditions2p2c
{
    enum Flags
        {
            saturation=1,
            concentration=2,
        };
};

}

#endif
