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
 * \ingroup SpatialParameters
 * \brief Base class for spatial parameters dealing with thermal and chemical non-equilibrium
 *
 */
#ifndef DUMUX_FV_SPATIALPARAMS_NONEQUILIBRIUM_HH
#define DUMUX_FV_SPATIALPARAMS_NONEQUILIBRIUM_HH

#include <dumux/material/spatialparams/fv.hh>

namespace Dumux {

/*!
 * \ingroup SpatialParameters
 * \brief Definition of the spatial parameters for non-equilibrium
 */
template<class GridGeometry, class Scalar, class Implementation>
class FVNonEquilibriumSpatialParams
: public FVSpatialParams<GridGeometry, Scalar, Implementation>
{
    using ParentType = FVSpatialParams<GridGeometry, Scalar, Implementation>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    //! export the types used for interfacial area calculations
    using AwnSurfaceParams = Scalar;
    using AwsSurfaceParams = Scalar;
    using AnsSurfaceParams = Scalar;

    FVNonEquilibriumSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {  }

    /*!
     * \brief Return the characteristic length for the mass transfer.
     *
     * The position is determined based on the coordinate of
     * the vertex belonging to the considered sub control volume.
     */
    template<class ElementSolution>
    const Scalar characteristicLength(const Element & element,
                                      const SubControlVolume &scv,
                                      const ElementSolution &elemSol) const

    { return this->asImp_().characteristicLengthAtPos(scv.dofPosition()); }

    /*!
     * \brief Return the characteristic length for the mass transfer.
     *
     * \param globalPos The position in global coordinates.
     */
    const Scalar characteristicLengthAtPos(const GlobalPosition & globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a characteristicLengthAtPos() method.");
    }

    /*!
     * \brief Return the pre-factor the the energy transfer
     *
     * The position is determined based on the coordinate of
     * the vertex belonging to the considered sub control volume.
     */
    template<class ElementSolution>
    const Scalar factorEnergyTransfer(const Element &element,
                                      const SubControlVolume &scv,
                                      const ElementSolution &elemSol) const
    { return this->asImp_().factorEnergyTransferAtPos(scv.dofPosition()); }

    /*!
     * \brief Return the pre factor the the energy transfer
     *
     * \param globalPos The position in global coordinates.
     */
    const Scalar factorEnergyTransferAtPos(const  GlobalPosition & globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a factorEnergyTransferAtPos() method.");
    }

    /*!
     * \brief Return the pre-factor the the mass transfer
     *
     * The position is determined based on the coordinate of
     * the vertex belonging to the considered sub control volume.
     */
    template<class ElementSolution>
    const Scalar factorMassTransfer(const Element &element,
                                      const SubControlVolume &scv,
                                      const ElementSolution &elemSol) const
    { return this->asImp_().factorMassTransferAtPos(scv.dofPosition()); }


    /*!
     * \brief Return the pre-factor the the mass transfer
     *
     * \param globalPos The position in global coordinates.
     */
    const Scalar factorMassTransferAtPos(const  GlobalPosition & globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a factorMassTransferAtPos() method.");
    }
};

} // end namespace Dumux

#endif
