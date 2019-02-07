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
 * \ingroup BoundaryTests
 * \brief The spatial parameters class for the test problem using the 2pnc cc model
 */
#ifndef DUMUX_CONSERVATION_SPATIAL_PARAMS_HH
#define DUMUX_CONSERVATION_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysomerton.hh>


namespace Dumux
{

/*!
 * \ingroup TwoPNCModel
 * \ingroup ImplicitTestProblems
 *
 * \brief The spatial parameters class for the test problem using the
 *        1p cc model
 */
template<class FVGridGeometry, class Scalar>
class TwoPTwoCSpatialParams
: public FVSpatialParams<FVGridGeometry, Scalar, TwoPTwoCSpatialParams<FVGridGeometry, Scalar>>
{
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using ParentType = FVSpatialParams<FVGridGeometry, Scalar, TwoPTwoCSpatialParams<FVGridGeometry, Scalar>>;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using EffectiveLaw = RegularizedVanGenuchten<Scalar>;

public:
    using MaterialLaw = EffToAbsLaw<EffectiveLaw>;
    using MaterialLawParams = typename MaterialLaw::Params;
    using PermeabilityType = Scalar;

    TwoPTwoCSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
        : ParentType(fvGridGeometry)
    {
        permeability_ = getParam<Scalar>("Darcy.SpatialParams.Permeability");
        alphaBJ_ = getParam<Scalar>("Darcy.SpatialParams.AlphaBeaversJoseph");
        porosity_ = getParam<Scalar>("Darcy.SpatialParams.Porosity");

        // residual saturations
        params_.setSwr(getParam<Scalar>("Darcy.SpatialParams.Swr"));
        params_.setSnr(getParam<Scalar>("Darcy.SpatialParams.Snr"));
        // parameters for the vanGenuchten law
        params_.setVgAlpha(getParam<Scalar>("Darcy.SpatialParams.VgAlpha"));
        params_.setVgn(getParam<Scalar>("Darcy.SpatialParams.VgN"));
        Scalar threshold = 0.01 * (1.0 - params_.swr() - params_.snr());
        params_.setPcLowSw(getParam<Scalar>("Darcy.SpatialParams.PcLowSw", threshold));
        params_.setPcHighSw(getParam<Scalar>("Darcy.SpatialParams.PcHighSw", 1.0-threshold));
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
     *
     * \param globalPos The global position
     * \return the intrinsic permeability
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return permeability_; }

    /*! \brief Defines the porosity in [-].
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }

    /*! \brief Defines the Beavers-Joseph coefficient in [-].
     *
     * \param globalPos The global position
     */
    Scalar beaversJosephCoeffAtPos(const GlobalPosition& globalPos) const
    { return alphaBJ_; }

    /*!
     * \brief Returns the parameter object for the Brooks-Corey material law.
     *        In this test, we use element-wise distributed material parameters.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return the material parameters object
     */
    template<class ElementSolutionVector>
    const MaterialLawParams& materialLawParams(const Element& element,
                                               const SubControlVolume& scv,
                                               const ElementSolutionVector& elemSol) const
    { return params_; }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \return the wetting phase index
     * \param globalPos The global position
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::phase0Idx; }

private:
    Scalar permeability_;
    Scalar alphaBJ_;
    Scalar porosity_;
    MaterialLawParams params_;
    static constexpr Scalar eps_ = 1.0e-7;
};

} // end namespace Dumux

#endif
