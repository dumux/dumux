// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief A test problem for the 1p model: A pipe system with circular cross-section
 *        and a branching point embedded in a three-dimensional world
 */

#ifndef DUMUX_ONEP_TUBES_TEST_SPATIALPARAMS_HH
#define DUMUX_ONEP_TUBES_TEST_SPATIALPARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief A test problem for the 1p model: A pipe system with circular cross-section
 *        and a branching point embedded in a three-dimensional world
 */
template<class GridGeometry, class Scalar>
class TubesTestSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                                         TubesTestSpatialParams<GridGeometry, Scalar>>
{
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using ThisType = TubesTestSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    TubesTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
        : ParentType(gridGeometry)
    {
        radius_ = 1.0;

        using std::sqrt;
        radiusMain_ = sqrt(sqrt(4.0/sqrt(3.0)));
    }

    /*!
     * \brief Returns the radius of the circular pipe for the current
     * sub-control volume in [m].
     *
     * \param scv The sub-control volume
     */
    Scalar radius(const SubControlVolume &scv) const
    {
        if(scv.center()[2] > 0.5 - eps_)
            return radiusMain_;
        else
            return radius_;
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
     *
     * \param element The element
     * \param scv The sub-control volume
     * \param elemSol The element solution vector
     * \return The intrinsic permeability
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        const Scalar radius = this->radius(scv);
        const Scalar gamma = 2; // quadratic velocity profile (Poiseuille flow)
        return radius*radius/(2*(2+gamma));
    }

    /*!
     * \brief Define the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1; }

    /*!
     * \brief Returns how much the domain is extruded at a given sub-control volume.
     *
     * This means the factor by which a lower-dimensional (1D or 2D)
     * entity needs to be expanded to get a full dimensional cell. The
     * default is 1.0 which means that 1D problems are actually
     * thought as pipes with a cross section of 1 m^2 and 2D problems
     * are assumed to extend 1 m to the back.
     */
    template<class ElementSolution>
    Scalar extrusionFactor(const Element &element,
                           const SubControlVolume &scv,
                           const ElementSolution& elemSol) const
    {
        const auto radius = this->radius(scv);
        return M_PI*radius*radius;
    }

private:
    Scalar radius_, radiusMain_;
    static constexpr Scalar eps_ = 1e-8;
};

} // end namespace Dumux

#endif
