// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTests
 * \brief The spatial params for the incompressible 1p test.
 */

#ifndef DUMUX_INCOMPRESSIBLE_ONEPBOXDFM_TEST_SPATIAL_PARAMS_HH
#define DUMUX_INCOMPRESSIBLE_ONEPBOXDFM_TEST_SPATIAL_PARAMS_HH

#include <dumux/discretization/method.hh>

#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup TwoPTests
 * \brief The spatial params for the incompressible 1p test.
 */
template<class GridGeometry, class Scalar>
class OnePTestSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP< GridGeometry, Scalar, OnePTestSpatialParams<GridGeometry, Scalar> >
{
    using ThisType = OnePTestSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

public:
    using PermeabilityType = Scalar;

    OnePTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , matrixPermeability_{getParam<Scalar>("SpatialParams.MatrixPermeability")}
    , fracturePermeability_{getParam<Scalar>("SpatialParams.FracturePermeability")}
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
        return scv.isOnFracture() ? aHalf: 1.0;
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
    { return scv.isOnFracture() ? fracturePermeability_ : matrixPermeability_; }

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
    { return scv.isOnFracture() ? 0.6 : 0.15; }

private:
    Scalar matrixPermeability_;
    Scalar fracturePermeability_;
};

} // end namespace Dumux

#endif
