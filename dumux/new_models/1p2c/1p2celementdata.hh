// $Id$

/*!
 * \file
 *
 * \brief Contains the quantities which are are constant within a
 *        finite element in the two-phase, two-component model.
 *
 * By for the plain 2p2c model everything is given on the finite
 * volumes, so this class is empty.
 */
#ifndef DUMUX_1P2C_ELEMENT_DATA_HH
#define DUMUX_1P2C_ELEMENT_DATA_HH

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
 * \brief This template class contains the quantities which are are
 *        constant within a finite element in the one-phase,
 *        two-component model.
 *
 * By for the plain 1p2c model everything is given on the finite
 * volumes, so this class is empty.
 */
template <class OnePTwoCTraits,
          class ProblemT>
class OnePTwoCElementData
{
};

} // end namepace

#endif
