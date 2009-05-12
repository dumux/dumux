// $Id$

/*!
 * \file
 *
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase, two-component model.
 */

#ifndef DUMUX_1P2C_VERTEX_DATA_HH
#define DUMUX_1P2C_VERTEX_DATA_HH

#include <dumux/new_models/boxscheme/boxscheme.hh>
#include <dumux/new_models/boxscheme/p1boxtraits.hh>
#include <dumux/new_models/1p2c/1p2ctraits.hh>
#include <dumux/auxiliary/math.hh>

#include <dumux/auxiliary/apis.hh>
#include <dune/common/collectivecommunication.hh>
#include <vector>
#include <iostream>

namespace Dune
{

/*!
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase, two-component model.
 */

template <class OnePTwoCTraits,
          class Problem>
class OnePTwoCVertexData
{
    typedef OnePTwoCTraits Tr;
    typedef typename Problem::DomainTraits::Scalar Scalar;
    typedef typename Problem::DomainTraits::Grid Grid;

    typedef typename Grid::template Codim<0>::Entity Element;

    typedef Dune::FieldVector<Scalar, Tr::numEq>      SolutionVector;
    typedef Dune::FieldVector<Scalar, Tr::numPhases>  PhasesVector;

    typedef Dune::FieldVector<Scalar, Grid::dimensionworld>  GlobalPosition;
    typedef Dune::FieldVector<Scalar, Grid::dimension>       LocalPosition;

public:
    /*!
     * \brief Update all quantities for a given control volume.
     */
    template <class JacobianImp>
    void update(const SolutionVector   &sol,
                const Element          &element,
                int                     vertIdx,
                bool                    isOldSol,
                JacobianImp            &jac)
    {
        const GlobalPosition &global = element.geometry().corner(vertIdx);
        const LocalPosition   &local =
            Problem::DomainTraits::referenceElement(element.type()).position(vertIdx,
                                                                             Grid::dimension);

        double T = 273;
        double p = 1e-5;

        porosity = jac.problem().soil().porosity(global,element,local);
        viscosity = jac.problem().phase().viscosity(T,p);
        tortuosity = jac.problem().soil().tortuosity(global,element,local);
        diffCoeff = jac.problem().phase().diffCoeff();
        molefraction = sol[Tr::transport];
        pressure = sol[Tr::konti];
    }


public:
    Scalar porosity;
    Scalar viscosity;
    Scalar tortuosity;
    Scalar diffCoeff;
    Scalar molefraction;
    Scalar pressure;
};


} // end namepace

#endif
