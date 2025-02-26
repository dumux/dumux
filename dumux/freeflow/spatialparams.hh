// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \brief Definition of the spatial parameters for the freeflow problems.
 */
#ifndef DUMUX_FREEFLOW_SPATIAL_PARAMS_HH
#define DUMUX_FREEFLOW_SPATIAL_PARAMS_HH

#include <dumux/common/fvspatialparams.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace Detail::BrinkmanSpatialParams {

template<class GlobalPosition>
struct hasInversePermeabilityAtPos
{
    template<class SpatialParams>
    auto operator()(const SpatialParams& a)
    -> decltype(a.inversePermeabilityAtPos(std::declval<GlobalPosition>()))
    {}
};

template<class GlobalPosition>
struct hasBrinkmanEpsilonAtPos
{
    template<class SpatialParams>
    auto operator()(const SpatialParams& a)
    -> decltype(a.brinkmanEpsilonAtPos(std::declval<GlobalPosition>()))
    {}
};

} // end namespace Detail
#endif

/*!
 * \ingroup FreeflowModels
 * \brief Definition of the spatial parameters for the freeflow problems.
 */
template<class GridGeometry, class Scalar, class Implementation>
class FreeFlowSpatialParams
: public FVSpatialParams<GridGeometry, Scalar, Implementation>
{
    using ParentType = FVSpatialParams<GridGeometry, Scalar, Implementation>;

public:
    FreeFlowSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}
};

/*!
 * \ingroup FreeflowModels
 * \brief Definition of the spatial parameters for the darcy-brinkman problems.
 */
template<class GridGeometry, class Scalar, class Implementation>
class BrinkmanSpatialParams
: public FreeFlowSpatialParams<GridGeometry, Scalar, Implementation>
{
    using ParentType = FreeFlowSpatialParams<GridGeometry, Scalar, Implementation>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;

public:
    BrinkmanSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

    decltype(auto) inversePermeability(const Element& element,
                                       const FVElementGeometry& fvGeometry,
                                       const SubControlVolume& scv) const
    {
        static_assert(decltype(isValid(Detail::BrinkmanSpatialParams::hasInversePermeabilityAtPos<GlobalPosition>())(this->asImp_()))::value," \n\n"
        "   Your spatial params class has to either implement\n\n"
        "         const PermeabilityType& inversePermeabilityAtPos(const GlobalPosition& globalPos) const\n\n"
        "   or overload this function\n\n"
        "         const PermeabilityType& inversePermeability(const Element& element,\n"
        "                                                     const FVElementGeometry& fvGeometry,\n"
        "                                                     const SubControlVolume& scv) const\n"
        "    This can be done simply with the invert() function of the DimMatrix type (e.g. permeability.invert()). \n\n");
        return this->asImp_().inversePermeabilityAtPos(scv.dofPosition());
    }

    Scalar brinkmanEpsilon(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const SubControlVolume& scv) const
    {
        static_assert(decltype(isValid(Detail::BrinkmanSpatialParams::hasBrinkmanEpsilonAtPos<GlobalPosition>())(this->asImp_()))::value," \n\n"
        "   Your spatial params class has to either implement\n\n"
        "         const Scalar& brinkmanEpsilonAtPos(const GlobalPosition& globalPos) const\n\n"
        "   or overload this function\n\n"
        "         const Scalar& brinkmanEpsilon(const Element& element,\n"
        "                                       const FVElementGeometry& fvGeometry,\n"
        "                                       const SubControlVolume& scv) const\n\n");
        return this->asImp_().brinkmanEpsilonAtPos(scv.dofPosition());
    }

};

/*!
 * \ingroup FreeflowModels
 * \brief Definition of the spatial parameters for the freeflow problems.
 */
template<class GridGeometry, class Scalar>
class FreeFlowDefaultSpatialParams
: public FreeFlowSpatialParams<GridGeometry, Scalar, FreeFlowDefaultSpatialParams<GridGeometry, Scalar>>
{
    using ParentType = FreeFlowSpatialParams<GridGeometry, Scalar, FreeFlowDefaultSpatialParams<GridGeometry, Scalar>>;
public:
    using ParentType::ParentType;
};

} // end namespace Dumux

#endif
