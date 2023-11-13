// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Darcy Brinkman model for a single-domain evaluation of coupled freeflow and porous medium flows
 */

#ifndef DUMUX_BRINKMAN_SPATIAL_PARAMS_HH
#define DUMUX_BRINKMAN_SPATIAL_PARAMS_HH

#include <dumux/freeflow/spatialparams.hh>
#include <dune/common/fmatrix.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief The spatial parameters class for the darcy-brinkman model
 */

template<class GridGeometry, class Scalar>
class BrinkmanSpatialParams
: public FreeFlowSpatialParams<GridGeometry, Scalar, BrinkmanSpatialParams<GridGeometry, Scalar>>
{
    using ParentType = FreeFlowSpatialParams<GridGeometry, Scalar, BrinkmanSpatialParams<GridGeometry, Scalar>>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    using PermeabilityType = DimWorldMatrix;

    BrinkmanSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry),
    permeability_(1.0),
    ffPermeability_(1.0)
    {
        isotropicK_ = getParam<Scalar>("SpatialParams.Permeability");

        permeability_[0][0] = isotropicK_;
        permeability_[1][1] = isotropicK_ * 1e4;

        pmLowerLeft_ = getParam<GlobalPosition>("SpatialParams.PorousMediumLowerLeft");
        pmUpperRight_ = getParam<GlobalPosition>("SpatialParams.PorousMediumUpperRight");

    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
     *
     * \param element The element
     * \param scv The sub control volume
     * \return the intrinsic permeability
     */
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv) const
    { return isPM_(scv.dofPosition()) ? permeability_ : ffPermeability_; }


    Scalar brinkmanEpsilon(const Element& element,
                           const SubControlVolume& scv) const
    { return isPM_(scv.dofPosition()) ? 1.0 : 0.0; }

private:
    bool isPM_(const GlobalPosition &globalPos) const
    {
        for (int i = 0; i < dimWorld; ++i)
        {
            if (globalPos[i] < pmLowerLeft_[i] + eps_
             || globalPos[i] > pmUpperRight_[i] - eps_)
                return false;
        }
        return true;
    }

    Scalar isotropicK_;
    PermeabilityType permeability_;
    PermeabilityType ffPermeability_;

    GlobalPosition pmLowerLeft_;
    GlobalPosition pmUpperRight_;

    static constexpr Scalar eps_ = 1e-7;
};

} // end namespace Dumux

#endif
