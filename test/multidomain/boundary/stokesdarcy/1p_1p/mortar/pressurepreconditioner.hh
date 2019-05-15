// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup TODO doc me.
 * \brief TODO doc me.
 */
#ifndef DUMUX_STOKES_DARCY_MORTAR_PRESSURE_PRECONDITIONER_HH
#define DUMUX_STOKES_DARCY_MORTAR_PRESSURE_PRECONDITIONER_HH

#include <dune/istl/preconditioner.hh>

namespace Dumux {

template<class MortarSolutionVector>
class MortarPressurePreconditioner
: public Dune::Preconditioner<MortarSolutionVector, MortarSolutionVector>
{
    using X = MortarSolutionVector;
    using Y = MortarSolutionVector;

 public:
   /*!
    * \brief Prepare the preconditioner.
    */
   virtual void pre (X& x, Y& b)
   {
       // TODO
   }

   /*!
    * \brief Apply one step of the preconditioner to the system A(v)=d.
    */
   virtual void apply (X& v, const Y& d)
   {
       v = d;
   }

   /*!
    * \brief Clean up.
    */
   virtual void post (X& x)
   {
       // TODO
   }

   //! Category of the preconditioner (see SolverCategory::Category)
   virtual Dune::SolverCategory::Category category() const
   {
       return Dune::SolverCategory::sequential;
   }
 };

} // end namespace Dumux

#endif
