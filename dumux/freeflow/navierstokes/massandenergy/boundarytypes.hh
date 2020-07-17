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
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesBoundaryTypes
 */
#ifndef FREEFLOW_NAVIERSTOKES_MASS_AND_ENERGY_BOUNDARY_TYPES_HH
#define FREEFLOW_NAVIERSTOKES_MASS_AND_ENERGY_BOUNDARY_TYPES_HH

#include <dune/common/exceptions.hh>
#include <dumux/common/boundarytypes.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Class to specify the type of a boundary condition for the Navier-Stokes model.
 */
template <int numEq, class Indices>
class NavierStokesMassAndEnergyBoundaryTypes : public BoundaryTypes<numEq>
{
    using ParentType = BoundaryTypes<numEq>;

public:
    using ParentType::ParentType;

    /*!
     * \brief  Prevent setting all boundary conditions to Dirichlet.
     */
    template<class T = void>
    void setAllDirichlet()
    {
        static_assert(AlwaysFalse<T>::value, "Setting all boundary types to Dirichlet not permitted!");
    }

    /*!
     * \brief Set a Dirichlet boundary condition for a single primary
     *        variable.
     *
     * Depending on the discretization, setting the Dirichlet condition
     * will replace the balance equation with index equal to pvIdx.
     *
     * \param pvIdx The index of the primary variable inside a
     *              PrimaryVariables object.
     */
    void setDirichlet(int pvIdx)
    {
        // if (pvIdx == Indices::pressureIdx)
        //     DUNE_THROW(Dune::InvalidStateException, "Setting Dirichlet for pressure in the mass model is not permitted");
        // else
            ParentType::setDirichlet(pvIdx);
    }

};

} // end namespace Dumux

#endif
