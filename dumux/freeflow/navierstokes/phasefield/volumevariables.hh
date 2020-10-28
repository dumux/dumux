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
 * \ingroup NIModel
 * \brief Base class for the model specific class which provides
 *        access to all volume averaged quantities.
 */

#ifndef DUMUX_NAVIER_STOKES_PHASEFIELD_VOLUME_VARIABLES_HH
#define DUMUX_NAVIER_STOKES_PHASEFIELD_VOLUME_VARIABLES_HH

#include <type_traits>
#include <dune/common/std/type_traits.hh>


namespace Dumux {

/*!
 * \ingroup NIModel
 * \brief The isothermal base class
 */
template<class Traits, class Impl>
class NavierStokesPhasefieldVolumeVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;
    static constexpr bool enablePhasefield = Traits::ModelTraits::enablePhasefield();

public:
    using FluidState = typename Traits::FluidState;
    using FluidSystem = typename Traits::FluidSystem;

    /*!
    * \brief Returns the phasefield at a given sub-control volume
    *
    * \param elemSol A vector containing all primary variables connected to the element
    * \param problem The object specifying the problem which ought to
    *                be simulated
    * \param element An element which contains part of the control volume
    * \param scv The sub-control volume
    */
    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    Scalar getPhasefield(const ElementSolution& elemSol,
                          const Problem& problem,
                          const Element& element,
                          const SubControlVolume& scv) const
    {
        if constexpr (enablePhasefield)
            return elemSol[scv.localDofIndex()][Traits::ModelTraits::Indices::phiIdx];
        else
            return 1.0;
    }

protected:
    const Impl &asImp_() const { return *static_cast<const Impl*>(this); }
    Impl &asImp_() { return *static_cast<Impl*>(this); }

    // \}

};

} // end namespace Dumux

#endif
