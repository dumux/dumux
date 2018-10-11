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
#ifndef DUMUX_IMPES_TRANSPORT_TEST_SPATIAL_PARAMS_HH
#define DUMUX_IMPES_TRANSPORT_TEST_SPATIAL_PARAMS_HH

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
class TwoPTransportSpatialParams
: public FVSpatialParams<FVGridGeometry, Scalar, TwoPTransportSpatialParams<FVGridGeometry, Scalar>>
{
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ThisType = TwoPTransportSpatialParams<FVGridGeometry, Scalar>;
    using ParentType = FVSpatialParams<FVGridGeometry, Scalar, ThisType>;

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using EffectiveLaw = RegularizedBrooksCorey<Scalar>;

public:
    using MaterialLaw = EffToAbsLaw<EffectiveLaw>;
    using MaterialLawParams = typename MaterialLaw::Params;
    using PermeabilityType = Scalar;

    TwoPTransportSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        // residual saturations
        materialLawParamsBackground_.setSwr(0.);
        materialLawParamsBackground_.setSnr(0.);

        materialLawParamsLenses_.setSwr(0.);
        materialLawParamsLenses_.setSnr(0.);

        //parameters for Brooks-Corey law

        // entry pressures
        materialLawParamsBackground_.setPe(getParam<Scalar>("SpatialParams.BackgroundEntryPressure", 0.0));
        materialLawParamsLenses_.setPe(getParam<Scalar>("SpatialParams.LenseEntryPressure", 0.0));

        materialLawParamsBackground_.setLambda(getParam<Scalar>("SpatialParams.BackgroundLambda", 3.0));
        materialLawParamsLenses_.setLambda(getParam<Scalar>("SpatialParams.LenseLambda", 2.0));

        permBackground_ = getParam<Scalar>("SpatialParams.BackgroundPermeability", 1e-10);
        permLenses_ = getParam<Scalar>("SpatialParams.LensPermeability", 1e-10);

        lensOneLowerLeft_ = getParam<GlobalPosition>("SpatialParams.LensOneLowerLeft", GlobalPosition(0.0));
        lensOneUpperRight_ = getParam<GlobalPosition>("SpatialParams.LensOneUpperRight", GlobalPosition(0.0));
        lensTwoLowerLeft_ = getParam<GlobalPosition>("SpatialParams.LensTwoLowerLeft", GlobalPosition(0.0));
        lensTwoUpperRight_ = getParam<GlobalPosition>("SpatialParams.LensTwoUpperRight", GlobalPosition(0.0));
        lensThreeLowerLeft_ = getParam<GlobalPosition>("SpatialParams.LensThreeLowerLeft", GlobalPosition(0.0));
        lensThreeUpperRight_ = getParam<GlobalPosition>("SpatialParams.LensThreeUpperRight", GlobalPosition(0.0));
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
        auto globalPos = element.geometry().center();
        // do not use a less permeable lens in the test with inverted wettability
        if (isLensOne(globalPos))
            return permLenses_;
        else if (isLensTwo(globalPos))
            return permLenses_;
        else if (isLensThree(globalPos))
            return permLenses_;
        else
            return permBackground_;
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.4; }

    /*!
     * \brief Returns the parameter object for the Brooks-Corey material law.
     *        In this test, we use element-wise distributed material parameters.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return the material parameters object
     */
    template<class ElementSolution>
    const MaterialLawParams& materialLawParams(const Element& element,
                                               const SubControlVolume& scv,
                                               const ElementSolution& elemSol) const
    {
        auto globalPos = element.geometry().center();
        if (isLensOne(globalPos))
            return materialLawParamsLenses_;
        else if (isLensTwo(globalPos))
            return materialLawParamsLenses_;
        else if (isLensThree(globalPos))
            return materialLawParamsLenses_;
        else
            return materialLawParamsBackground_;
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

    //! velocity field
    template<class ElementVolumeVariables>
    Scalar volumeFlux(const Element &element,
                      const FVElementGeometry& fvGeometry,
                      const ElementVolumeVariables& elemVolVars,
                      const SubControlVolumeFace& scvf) const
    {
//        GlobalPosition vel(0.0);
//        vel[dimWorld-1] = -1.0e-6;
//        return vel*scvf.unitOuterNormal();
        return volumeFlux_[scvf.index()];
    }

    void setVolumeFlux(const std::vector<Scalar>& f)
    { volumeFlux_ = f; }

private:
    bool isLensOne(const GlobalPosition& globalPos) const
    {
        for (int i = 0; i < dimWorld; i++)
        {
            if (globalPos[i] < lensOneLowerLeft_[i] - eps_ || globalPos[i] > lensOneUpperRight_[i] + eps_)
            {
                return false;
            }
        }
        return true;
    }
    bool isLensTwo(const GlobalPosition& globalPos) const
    {
        for (int i = 0; i < dimWorld; i++)
        {
            if (globalPos[i] < lensTwoLowerLeft_[i] - eps_ || globalPos[i] > lensTwoUpperRight_[i] + eps_)
            {
                return false;
            }
        }
        return true;
    }
    bool isLensThree(const GlobalPosition& globalPos) const
    {
        for (int i = 0; i < dimWorld; i++)
        {
            if (globalPos[i] < lensThreeLowerLeft_[i] - eps_ || globalPos[i] > lensThreeUpperRight_[i] + eps_)
            {
                return false;
            }
        }
        return true;
    }

    MaterialLawParams materialLawParamsBackground_;
    MaterialLawParams materialLawParamsLenses_;
    PermeabilityType permBackground_;
    PermeabilityType permLenses_;
    GlobalPosition lensOneLowerLeft_;
    GlobalPosition lensOneUpperRight_;
    GlobalPosition lensTwoLowerLeft_;
    GlobalPosition lensTwoUpperRight_;
    GlobalPosition lensThreeLowerLeft_;
    GlobalPosition lensThreeUpperRight_;

    std::vector<Scalar> volumeFlux_;

    static constexpr Scalar eps_ = 1.5e-7;
};

} // end namespace Dumux

#endif
