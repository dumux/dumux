// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief base class for spatialparameters dealing with thermal and chemical nonequilibrium
 *
 */
#ifndef DUMUX_FV_SPATIALPARAMS_NONEQUILIBRIUM_HH
#define DUMUX_FV_SPATIALPARAMS_NONEQUILIBRIUM_HH

#include <dumux/material/spatialparams/fv.hh>

// material laws for interfacial area
#include <dumux/material/fluidmatrixinteractions/2pia/efftoabslawia.hh>
#include <dumux/material/fluidmatrixinteractions/2pia/awnsurfacepolynomial2ndorder.hh>
#include <dumux/material/fluidmatrixinteractions/2pia/awnsurfacepcmaxfct.hh>
#include <dumux/material/fluidmatrixinteractions/2pia/awnsurfaceexpswpcto3.hh>

namespace Dumux {

/**
 * \brief Definition of the spatial parameters for nonequilibrium
 */
template<class FVGridGeometry, class Scalar, class Implementation>
class FVNonEquilibriumSpatialParams
: public FVSpatialParams<FVGridGeometry,
                         Scalar,
                         Implementation>
{

    using ParentType = FVSpatialParams<FVGridGeometry, Scalar, Implementation>;
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    //! export the types used for interfacial area calculations

    using AwnSurfaceParams = Scalar;
    using AwsSurfaceParams = Scalar;
    using AnsSurfaceParams = Scalar;

    FVNonEquilibriumSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {  }

    /*!\brief Return a reference to the container object for the
     *        parametrization of the surface between wetting and non-Wetting phase.
     *
     *        The position is determined based on the coordinate of
     *        the vertex belonging to the considered sub control volume.
     * \param element     The finite element
     * \param fvGeometry  The finite volume geometry
     * \param scvIdx      The local index of the sub-control volume */
    template<class ElementSolution>
    AwnSurfaceParams aWettingNonWettingSurfaceParams(const Element &element,
                                                   const SubControlVolume &scv,
                                                   const ElementSolution &elemSol) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide a aWettingNonWettingSurfaceParams() method.");
    }

    /*!\brief Return a reference to the container object for the
     *        parametrization of the surface between non-Wetting and solid phase.
     *
     *        The position is determined based on the coordinate of
     *        the vertex belonging to the considered sub control volume.
     * \param element     The finite element
     * \param fvGeometry  The finite volume geometry
     * \param scvIdx      The local index of the sub-control volume */
    template<class ElementSolution>
    AnsSurfaceParams aNonWettingSolidSurfaceParams(const Element &element,
                                                 const SubControlVolume &scv,
                                                 const ElementSolution &elemSol) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide a aNonWettingSolidSurfaceParams() method.");
    }

    /*!\brief Return a reference to the container object for the
     *        parametrization of the surface between non-Wetting and solid phase.
     *
     *        The position is determined based on the coordinate of
     *        the vertex belonging to the considered sub control volume.
     * \param element     The finite element
     * \param fvGeometry  The finite volume geometry
     * \param scvIdx      The local index of the sub-control volume */
    template<class ElementSolution>
    AwsSurfaceParams aWettingSolidSurfaceParams(const Element &element,
                                              const SubControlVolume &scv,
                                              const ElementSolution &elemSol) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide a aWettingSolidSurfaceParams() method.");
    }

    /*!\brief Return the maximum capillary pressure for the given pc-Sw curve
     *
     *        Of course physically there is no such thing as a maximum capillary pressure.
     *        The parametrization (VG/BC) goes to infinity and physically there is only one pressure
     *        for single phase conditions.
     *        Here, this is used for fitting the interfacial area surface: the capillary pressure,
     *        where the interfacial area is zero.
     *        Technically this value is obtained as the capillary pressure of saturation zero.
     *        This value of course only exists for the case of a regularized pc-Sw relation.
     * \param element     The finite element
     * \param fvGeometry  The finite volume geometry
     * \param scvIdx      The local index of the sub-control volume */
    template<class ElementSolution>
    const Scalar pcMax(const Element &element,
                       const SubControlVolume &scv,
                       const ElementSolution &elemSol) const
    {     DUNE_THROW(Dune::InvalidStateException,
        "The spatial parameters do not provide a aNonWettingSolidSurfaceParams() method.");
    }


    /*!\brief Return the characteristic length for the mass transfer.
     *
     *        The position is determined based on the coordinate of
     *        the vertex belonging to the considered sub controle volume.
     * \param element     The finite element
     * \param fvGeometry  The finite volume geometry
     * \param scvIdx      The local index of the sub control volume */
    template<class ElementSolution>
    const Scalar characteristicLength(const Element & element,
                                      const SubControlVolume &scv,
                                      const ElementSolution &elemSol) const

    { return this->asImp_().characteristicLengthAtPos(scv.dofPosition()); }

    /*!\brief Return the characteristic length for the mass transfer.
     * \param globalPos The position in global coordinates.*/
    const Scalar characteristicLengthAtPos(const GlobalPosition & globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a characteristicLengthAtPos() method.");
    }

    /*!\brief Return the pre factor the the energy transfer
     *
     *        The position is determined based on the coordinate of
     *        the vertex belonging to the considered sub controle volume.
     * \param element     The finite element
     * \param fvGeometry  The finite volume geometry
     * \param scvIdx      The local index of the sub control volume */
    template<class ElementSolution>
    const Scalar factorEnergyTransfer(const Element &element,
                                      const SubControlVolume &scv,
                                      const ElementSolution &elemSol) const
    { return this->asImp_().factorEnergyTransferAtPos(scv.dofPosition()); }

    /*!\brief Return the pre factor the the energy transfer
     * \param globalPos The position in global coordinates.*/
    const Scalar factorEnergyTransferAtPos(const  GlobalPosition & globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a factorEnergyTransferAtPos() method.");
    }

    /*!\brief Return the pre factor the the mass transfer
     *
     *        The position is determined based on the coordinate of
     *        the vertex belonging to the considered sub controle volume.
     * \param element     The finite element
     * \param fvGeometry  The finite volume geometry
     * \param scvIdx      The local index of the sub control volume */
    template<class ElementSolution>
    const Scalar factorMassTransfer(const Element &element,
                                      const SubControlVolume &scv,
                                      const ElementSolution &elemSol) const
    { return this->asImp_().factorMassTransferAtPos(scv.dofPosition()); }


    /*!\brief Return the pre factor the the mass transfer
     * \param globalPos The position in global coordinates.*/
    const Scalar factorMassTransferAtPos(const  GlobalPosition & globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a factorMassTransferAtPos() method.");
    }
};

}

#endif // GUARDIAN
