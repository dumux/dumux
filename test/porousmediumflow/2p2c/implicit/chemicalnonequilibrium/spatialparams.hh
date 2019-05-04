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
 * \ingroup TwoPTwoCTests
 * \brief The spatial parameters for the 2p2c chemical nonequilibrium problem.
 */

#ifndef DUMUX_MPNC_COMPARISON_SPATIAL_PARAMS_HH
#define DUMUX_MPNC_COMPARISON_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/spatialparams/fvnonequilibrium.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedlinearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>

// material laws for interfacial area
#include <dumux/material/fluidmatrixinteractions/2pia/efftoabslawia.hh>
#include <dumux/material/fluidmatrixinteractions/2pia/awnsurfacepolynomial2ndorder.hh>
#include <dumux/material/fluidmatrixinteractions/2pia/awnsurfacepcmaxfct.hh>
#include <dumux/material/fluidmatrixinteractions/2pia/awnsurfaceexpswpcto3.hh>
namespace Dumux {

/**
 * \ingroup TwoPTwoCTests
 * \brief The spatial parameters for the 2p2c chemical nonequilibrium problem.
 */
template<class FVGridGeometry, class Scalar>
class TwoPTwoCChemicalNonequilibriumSpatialParams
: public FVNonEquilibriumSpatialParams<FVGridGeometry, Scalar,
                                       TwoPTwoCChemicalNonequilibriumSpatialParams<FVGridGeometry, Scalar>>
{
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVNonEquilibriumSpatialParams<FVGridGeometry, Scalar,
                                                     TwoPTwoCChemicalNonequilibriumSpatialParams<FVGridGeometry, Scalar>>;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    enum {dimWorld=GridView::dimensionworld};

    using EffectiveLaw = RegularizedBrooksCorey<Scalar>;

public:
    //! Export permeability type
    using PermeabilityType = Scalar;
    //! Export the type used for the material law
    using MaterialLaw = EffToAbsLaw<EffectiveLaw>;
    using MaterialLawParams = typename MaterialLaw::Params;
    using EffectiveIALawAwn = AwnSurfacePcMaxFct<Scalar>;
    using AwnSurface = EffToAbsLawIA<EffectiveIALawAwn, MaterialLawParams>;
    using AwnSurfaceParams = typename AwnSurface::Params;

    TwoPTwoCChemicalNonequilibriumSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        // intrinsic permeabilities
        coarseK_ = 1e-11;

        // the porosity
        porosity_ = 0.4;

         // residual saturations
        coarseMaterialParams_.setSwr(0.2);
        coarseMaterialParams_.setSnr(0.1);

        // parameters for the Brooks-Corey law
        coarseMaterialParams_.setPe(1e4);
        coarseMaterialParams_.setLambda(2.0);

        aWettingNonWettingA1_ = getParam<Scalar>("SpatialParams.WettingNonWettingAreaA1");
        aWettingNonWettingA2_ = getParam<Scalar>("SpatialParams.WettingNonWettingAreaA2");
        aWettingNonWettingA3_ = getParam<Scalar>("SpatialParams.WettingNonWettingAreaA3");

        // wetting-non wetting: surface which goes to zero on the edges, but is a polynomial
        aWettingNonWettingSurfaceParams_.setA1(aWettingNonWettingA1_);
        aWettingNonWettingSurfaceParams_.setA2(aWettingNonWettingA2_);
        aWettingNonWettingSurfaceParams_.setA3(aWettingNonWettingA3_);
        // determine maximum capillary pressure for wetting-nonwetting surface
        using TwoPLaw = EffToAbsLaw<RegularizedBrooksCorey<Scalar>>;
        pcMax_ = TwoPLaw::pc(coarseMaterialParams_, /*sw = */0.0);
        aWettingNonWettingSurfaceParams_.setPcMax(pcMax_);
        characteristicLength_ =getParam<Scalar>("SpatialParams.MeanPoreSize");
        factorMassTransfer_ = getParam<Scalar>("SpatialParams.MassTransferFactor");
    }

    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        return coarseK_;
    }

    /*!
     * \brief Defines the porosity \f$[-]\f$ of the soil
     *
     * \param globalPos The global position of the sub-control volume.
     * \return the material parameters object
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        return porosity_;
    }

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-sw, pc-sw, etc.).
     *
     * \param globalPos The global position of the sub-control volume.
     * \return The material parameters object
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    {
        return coarseMaterialParams_;
    }

    /*!\brief Returns a reference to the container object for the
     *        parametrization of the surface between wetting and non-Wetting phase.
     *
     * The position is determined based on the coordinate of
     * the vertex belonging to the considered sub-control volume.
     *
     * \param element The finite element
     * \param scv The sub-control volume
     * \param elemSol The element solution
     */
    template<class ElementSolution>
    const AwnSurfaceParams& aWettingNonWettingSurfaceParams(const Element &element,
                                                            const SubControlVolume &scv,
                                                            const ElementSolution &elemSol) const
    {
        return aWettingNonWettingSurfaceParams_ ;
    }

    /*!\brief Returns the maximum capillary pressure for the given pc-Sw curve
     *
     * Of course physically there is no such thing as a maximum capillary pressure.
     * The parametrization (VG/BC) goes to infinity and physically there is only one pressure
     * for single phase conditions.
     * Here, this is used for fitting the interfacial area surface: the capillary pressure,
     * where the interfacial area is zero.
     * Technically this value is obtained as the capillary pressure of saturation zero.
     * This value of course only exists for the case of a regularized pc-Sw relation.
     * \param element The finite element
     * \param scv The sub-control volume
     * \param elemSol The element solution
     */
    template<class ElementSolution>
    const Scalar pcMax(const Element &element,
                       const SubControlVolume &scv,
                       const ElementSolution &elemSol) const
    { return aWettingNonWettingSurfaceParams_.pcMax() ; }

    /*!
     * \brief Returns the characteristic length for the mass transfer.
     *
     * \param globalPos The position in global coordinates.
     */
    const Scalar characteristicLengthAtPos(const  GlobalPosition & globalPos) const
    { return characteristicLength_ ; }

    /*!
     * \brief Return the pre factor the the energy transfer.
     *
     * \param globalPos The position in global coordinates.
     */
    const Scalar factorMassTransferAtPos(const  GlobalPosition & globalPos) const
    { return factorMassTransfer_; }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \param globalPos The global position of the sub-control volume.
     * \return The wetting phase index
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::H2OIdx; }

private:

    Scalar coarseK_;
    Scalar porosity_;
    MaterialLawParams coarseMaterialParams_;
    static constexpr Scalar eps_ = 1e-6;

    AwnSurfaceParams aWettingNonWettingSurfaceParams_;
    Scalar pcMax_ ;

    // interfacial area parameters
    Scalar aWettingNonWettingA1_ ;
    Scalar aWettingNonWettingA2_ ;
    Scalar aWettingNonWettingA3_ ;

    Scalar factorMassTransfer_ ;
    Scalar characteristicLength_ ;
};

} // end namespace Dumux

#endif
