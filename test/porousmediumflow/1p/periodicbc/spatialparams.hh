// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief The spatial parameters for the incompressible test.
 */
#ifndef DUMUX_INCOMPRESSIBLE_ONEP_TEST_SPATIAL_PARAMS_HH
#define DUMUX_INCOMPRESSIBLE_ONEP_TEST_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief The spatial parameters class for the test problem using the
 *        incompressible 1p model.
 */
template<class TypeTag>
class OnePTestSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GetPropType<TypeTag, Properties::GridGeometry>,
                                         GetPropType<TypeTag, Properties::Scalar>,
                                         OnePTestSpatialParams<TypeTag>>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    using ThisType = OnePTestSpatialParams<TypeTag>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using PermeabilityType = Scalar;
    OnePTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        permeability_ = getParam<Scalar>("SpatialParams.Permeability");
        permeabilityLens_ = getParam<Scalar>("SpatialParams.PermeabilityLens");

        lensLowerLeft_ = getParam<GlobalPosition>("SpatialParams.LensLowerLeft");
        lensUpperRight_ = getParam<GlobalPosition>("SpatialParams.LensUpperRight");
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
     *
     * \param element The element
     * \param scv The sub-control volume
     * \param elemSol The element solution vector
     * \return the intrinsic permeability
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        if (isInLens_(scv.dofPosition()))
            return permeabilityLens_;
        else
            return permeability_;
    }

    /*!
     * \brief Defines the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.4; }

private:
    bool isInLens_(const GlobalPosition &globalPos) const
    {
        for (int i = 0; i < dimWorld; ++i) {
            if (globalPos[i] < lensLowerLeft_[i] + eps_ || globalPos[i] > lensUpperRight_[i] - eps_)
                return false;
        }
        return true;
    }

    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;

    Scalar permeability_, permeabilityLens_;

    static constexpr Scalar eps_ = 1.5e-7;
};

} // end namespace Dumux

#endif
