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

#ifndef DUMUX_INCOMPRESSIBLE_TWOPBOXDFM_TEST_SPATIAL_PARAMS_HH
#define DUMUX_INCOMPRESSIBLE_TWOPBOXDFM_TEST_SPATIAL_PARAMS_HH

#include <dumux/discretization/method.hh>

#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>

#include <dumux/porousmediumflow/2p/boxmaterialinterfaces.hh>

namespace Dumux {

/*!
 * \ingroup TwoPTests
 * \brief The spatial params for the incompressible 2p test.
 */
template<class GridGeometry, class Scalar>
class TwoPTestSpatialParams
: public FVPorousMediumFlowSpatialParamsMP< GridGeometry, Scalar, TwoPTestSpatialParams<GridGeometry, Scalar> >
{
    using ThisType = TwoPTestSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, ThisType>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    static constexpr int dimWorld = GridView::dimensionworld;

    using PcKrSw = FluidMatrix::BrooksCoreyDefault<Scalar>;
    using MaterialInterfaces = BoxMaterialInterfaces<GridGeometry, PcKrSw>;
public:
    using PermeabilityType = Scalar;

    TwoPTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , pcKrSwMatrix_("SpatialParams.Matrix")
    , pcKrSwFracture_("SpatialParams.Fracture")
    {}

    /*!
     * \brief Returns how much the domain is extruded at a given sub-control volume.
     *        Here, we extrude the fracture scvs by half the aperture
     */
    template<class ElementSolution>
    Scalar extrusionFactor(const Element& element,
                           const SubControlVolume& scv,
                           const ElementSolution& elemSol) const
    {
        // In the box-scheme, we compute fluxes etc element-wise,
        // thus per element we compute only half a fracture !!!
        static const Scalar aHalf = getParam<Scalar>("SpatialParams.FractureAperture")/2.0;
        if (scv.isOnFracture())
            return aHalf;
        return 1.0;
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
        if (scv.isOnFracture())
            return 5e-10;
        else
            return 1e-12;
    }


    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return The porosity
     */
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        if (scv.isOnFracture())
            return 0.6;
        else
            return 0.15;
    }

    /*!
     * \brief  Returns the fluid-matrix interaction law.
     *
     * In this test, we use element-wise distributed material parameters.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     */
    template<class ElementSolution>
    auto fluidMatrixInteraction(const Element& element,
                                const SubControlVolume& scv,
                                const ElementSolution& elemSol) const
    {
        if (scv.isOnFracture())
            return makeFluidMatrixInteraction(pcKrSwFracture_);
        else
            return makeFluidMatrixInteraction(pcKrSwMatrix_);
    }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \param globalPos The global position
     * \return The wetting phase index
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::phase0Idx; }

    //! Updates the map of which material parameters are associated with a nodal dof.
    template<class SolutionVector>
    void updateMaterialInterfaces(const SolutionVector& x)
    { materialInterfaces_ = std::make_unique<MaterialInterfaces>(this->gridGeometry(), *this, x); }

    //! Returns the material parameters associated with a nodal dof
    const MaterialInterfaces& materialInterfaces() const
    { return *materialInterfaces_; }

private:
    PcKrSw pcKrSwMatrix_;
    PcKrSw pcKrSwFracture_;

    std::unique_ptr<MaterialInterfaces> materialInterfaces_;

    static constexpr Scalar eps_ = 1.5e-7;
};

} // end namespace Dumux

#endif
