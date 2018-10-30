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
 * \ingroup TwoPTests
 * \brief The spatial params for the incompressible 2p test
 */
#ifndef DUMUX_MCWHORTER_PRESSURE_TEST_SPATIAL_PARAMS_HH
#define DUMUX_MCWHORTER_PRESSURE_TEST_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/porousmediumflow/2p/boxmaterialinterfaceparams.hh>

namespace Dumux {

/*!
 * \ingroup TwoPTests
 * \brief The spatial params for the incompressible 2p test
 */
template<class FVGridGeometry, class Scalar>
class TwoPTestSpatialParams
: public FVSpatialParams<FVGridGeometry, Scalar, TwoPTestSpatialParams<FVGridGeometry, Scalar>>
{
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ThisType = TwoPTestSpatialParams<FVGridGeometry, Scalar>;
    using ParentType = FVSpatialParams<FVGridGeometry, Scalar, ThisType>;

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using EffectiveLaw = RegularizedBrooksCorey<Scalar>;

public:
    using MaterialLaw = EffToAbsLaw<EffectiveLaw>;
    using MaterialLawParams = typename MaterialLaw::Params;
    using PermeabilityType = Scalar;

    TwoPTestSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        materialLawParams_.setSwr(0.2);
        materialLawParams_.setSnr(0.2);


        //parameters for Brooks-Corey law

        // entry pressures
        materialLawParams_.setPe(getParam<Scalar>("SpatialParams.EntryPressure",8000));


        materialLawParams_.setLambda(getParam<Scalar>("SpatialParams.Lambda",3));


        perm_ = getParam<Scalar>("SpatialParams.Permeability", 1.01936799e-13 );



        wettingSaturation_.resize(fvGridGeometry->numDofs(),0.0);
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
     *        In this test, we use element-wise distributed permeabilities.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return permeability
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {

        // do not use a less permeable lens in the test with inverted wettability

            return perm_;
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.2; }

    /*!
     * \brief Returns the parameter object for the Brooks-Corey material law
     *
     * \param globalPos The global position
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    {
        return materialLawParams_;
    }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \return the wetting phase index
     * \param globalPos The global position
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    {
        return FluidSystem::phase0Idx;
    }

    Scalar wettingSaturation(const Element& element,
                           const SubControlVolume& scv) const
    {
        return wettingSaturation_[scv.dofIndex()];
    }

    template<class Scvf>
    Scalar wettingSaturation(const Element& element,
                             const SubControlVolume& scv,
                             const Scvf& scvf) const
    {
        return wettingSaturation_[scvf.outsideScvIdx()];
    }

    void setWettingSaturation(const std::vector<Scalar>& sat)
    { wettingSaturation_ = sat;}

private:

    MaterialLawParams materialLawParams_;

    PermeabilityType perm_;


    std::vector<Scalar> wettingSaturation_;

    static constexpr Scalar eps_ = 1e-8;
};

} // end namespace Dumux

#endif
