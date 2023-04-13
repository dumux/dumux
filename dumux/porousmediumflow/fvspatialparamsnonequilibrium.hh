// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup SpatialParameters
 * \brief Base class for spatial parameters dealing with thermal and chemical non-equilibrium
 *
 */
#ifndef DUMUX_POROUS_MEDIUM_FLOW_FV_SPATIAL_PARAMS_NONEQUILIBRIUM_HH
#define DUMUX_POROUS_MEDIUM_FLOW_FV_SPATIAL_PARAMS_NONEQUILIBRIUM_HH

#include "fvspatialparamsmp.hh"

namespace Dumux {

/*!
 * \ingroup SpatialParameters
 * \brief Definition of the spatial parameters for non-equilibrium
 */
template<class GridGeometry, class Scalar, class Implementation>
class FVPorousMediumFlowSpatialParamsNonEquilibrium
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, Implementation>
{
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, Implementation>;
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

    FVPorousMediumFlowSpatialParamsNonEquilibrium(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

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
    const Scalar characteristicLengthAtPos(const GlobalPosition& globalPos) const
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
    const Scalar factorEnergyTransfer(const Element& element,
                                      const SubControlVolume& scv,
                                      const ElementSolution& elemSol) const
    { return this->asImp_().factorEnergyTransferAtPos(scv.dofPosition()); }

    /*!
     * \brief Return the pre factor the the energy transfer
     * \param globalPos The position in global coordinates.
     */
    const Scalar factorEnergyTransferAtPos(const GlobalPosition& globalPos) const
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
    const Scalar factorMassTransfer(const Element& element,
                                    const SubControlVolume& scv,
                                    const ElementSolution& elemSol) const
    { return this->asImp_().factorMassTransferAtPos(scv.dofPosition()); }


    /*!
     * \brief Return the pre-factor the the mass transfer
     * \param globalPos The position in global coordinates.
     */
    const Scalar factorMassTransferAtPos(const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a factorMassTransferAtPos() method.");
    }
};

} // end namespace Dumux

#endif
