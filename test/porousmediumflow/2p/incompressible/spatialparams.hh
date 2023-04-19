// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTests
 * \brief The spatial params for the incompressible 2p test.
 */

#ifndef DUMUX_INCOMPRESSIBLE_TWOP_TEST_SPATIAL_PARAMS_HH
#define DUMUX_INCOMPRESSIBLE_TWOP_TEST_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>
#include <dumux/porousmediumflow/2p/boxmaterialinterfaces.hh>

namespace Dumux {

/*!
 * \ingroup TwoPTests
 * \brief The spatial params for the incompressible 2p test.
 */
template<class GridGeometry, class Scalar>
class TwoPTestSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, TwoPTestSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using ThisType = TwoPTestSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, ThisType>;

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PcKrSw = FluidMatrix::VanGenuchtenDefault<Scalar>;
    using MaterialInterfaces = BoxMaterialInterfaces<GridGeometry, PcKrSw>;

public:
    using PermeabilityType = Scalar;

    TwoPTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , lensPcKrSw_("SpatialParams.Lens")
    , outerPcKrSw_("SpatialParams.Outer")
    {
        lensIsOilWet_ = getParam<bool>("SpatialParams.LensIsOilWet", false);

        lensLowerLeft_ = getParam<GlobalPosition>("SpatialParams.LensLowerLeft");
        lensUpperRight_ = getParam<GlobalPosition>("SpatialParams.LensUpperRight");

        lensK_ = getParam<Scalar>("SpatialParams.Lens.Permeability", 9.05e-12);
        outerK_ = getParam<Scalar>("SpatialParams.Outer.Permeability", 4.6e-10);
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
     *        In this test, we use element-wise distributed permeabilities.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return The permeability
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {

        // do not use a less permeable lens in the test with inverted wettability
        if (isInLens_(element.geometry().center()) && !lensIsOilWet_)
            return lensK_;
        return outerK_;
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
     *
     * In this test, we use element-wise distributed material parameters.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return The material parameters object
     */
    template<class ElementSolution>
    auto fluidMatrixInteraction(const Element& element,
                                const SubControlVolume& scv,
                                const ElementSolution& elemSol) const
    {
        // do not use different parameters in the test with inverted wettability
        if (isInLens_(element.geometry().center()) && !lensIsOilWet_)
            return makeFluidMatrixInteraction(lensPcKrSw_);
        return makeFluidMatrixInteraction(outerPcKrSw_);
    }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \param globalPos The global position
     * \return The wetting phase index
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    {
        if (isInLens_(globalPos) && lensIsOilWet_)
            return FluidSystem::phase1Idx;
        return FluidSystem::phase0Idx;
    }

    //! Updates the map of which material parameters are associated with a nodal dof.
    template<class SolutionVector>
    void updateMaterialInterfaces(const SolutionVector& x)
    {
        if (GridGeometry::discMethod == DiscretizationMethods::box)
            materialInterfaces_ = std::make_unique<MaterialInterfaces>(this->gridGeometry(), *this, x);
    }

    //! Returns the material parameters associated with a nodal dof
    const MaterialInterfaces& materialInterfaces() const
    { return *materialInterfaces_; }

    //! Returns whether or not the lens is oil wet
    bool lensIsOilWet() const { return lensIsOilWet_; }

private:
    bool isInLens_(const GlobalPosition &globalPos) const
    {
        for (int i = 0; i < dimWorld; ++i) {
            if (globalPos[i] < lensLowerLeft_[i] + eps_ || globalPos[i] > lensUpperRight_[i] - eps_)
                return false;
        }
        return true;
    }

    bool lensIsOilWet_;
    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;

    Scalar lensK_;
    Scalar outerK_;

    const PcKrSw lensPcKrSw_;
    const PcKrSw outerPcKrSw_;

    std::unique_ptr<MaterialInterfaces> materialInterfaces_;

    static constexpr Scalar eps_ = 1.5e-7;
};

} // end namespace Dumux

#endif
