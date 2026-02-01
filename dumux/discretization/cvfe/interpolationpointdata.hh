// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEDiscretization
 * \brief Classes representing interpolation point data for control-volume finite element schemes
 */
#ifndef DUMUX_CVFE_IP_DATA_HH
#define DUMUX_CVFE_IP_DATA_HH

#include <type_traits>

namespace Dumux::CVFE {

/*!
 * \ingroup CVFEDiscretization
 * \brief An interpolation point related to an element that includes global and local positions
 */
template<class LocalPosition, class GlobalPosition>
class InterpolationPointData
{
public:
    InterpolationPointData(LocalPosition&& localPos, GlobalPosition&& pos) : local_(std::move(localPos)), global_(std::move(pos)) {}
    InterpolationPointData(const LocalPosition& localPos, const GlobalPosition& pos) :  local_(localPos), global_(pos) {}

    //! The global position of the quadrature point
    const GlobalPosition& global() const
    { return global_; }

    //! The local position of the quadrature point
    const LocalPosition& local() const
    { return local_; }


private:
    LocalPosition local_;
    GlobalPosition global_;
};

/*!
 * \ingroup CVFEDiscretization
 * \brief An interpolation point related to a localDof of an element, giving its global and local positions
 */
template<class LocalPosition, class GlobalPosition, class LocalIndex>
class LocalDofInterpolationPointData : public InterpolationPointData<LocalPosition, GlobalPosition>
{
    using ParentType = InterpolationPointData<LocalPosition, GlobalPosition>;
public:
    LocalDofInterpolationPointData(LocalPosition&& localPos, GlobalPosition&& pos, LocalIndex index)
    : ParentType(localPos, pos), localDofIndex_(index) {}
    LocalDofInterpolationPointData(const LocalPosition& localPos, const GlobalPosition& pos, LocalIndex index)
    : ParentType(localPos, pos), localDofIndex_(index) {}

    //! The local index of the corresponding dof
    LocalIndex localDofIndex() const
    { return localDofIndex_; }

private:
    LocalIndex localDofIndex_;
};

/*!
 * \ingroup CVFEDiscretization
 * \brief An interpolation point related to a global position of an element, giving its local positions by a mapping
 */
template<class LocalMapping, class GlobalPosition>
class InterpolationPointDataLocalMapping
{
    using LocalPosition = std::invoke_result_t<LocalMapping, const GlobalPosition&>;

public:
    InterpolationPointDataLocalMapping(LocalMapping&& mapping, GlobalPosition&& pos) : localMapping_(std::move(mapping)), global_(std::move(pos)) {}
    InterpolationPointDataLocalMapping(LocalMapping&& mapping, const GlobalPosition& pos) : localMapping_(std::move(mapping)), global_(pos) {}

    //! The global position of the quadrature point
    const GlobalPosition& global() const
    { return global_; }

    //! The local position of the quadrature point
    const LocalPosition local() const
    { return localMapping_(global_); }


private:
    LocalMapping localMapping_;
    GlobalPosition global_;
};

/*!
 * \ingroup CVFEDiscretization
 * \brief An interpolation point related to a face of an element
 */
template<class BaseClass, class LocalIndex>
class FaceInterpolationPointData : public BaseClass
{
    using GlobalPosition = std::remove_cvref_t<decltype(std::declval<BaseClass>().global())>;
    using LocalPosition = std::remove_cvref_t<decltype(std::declval<BaseClass>().local())>;

public:
    template<class... Args>
    FaceInterpolationPointData(GlobalPosition&& n, LocalIndex index, Args&&... args)
    : BaseClass(std::forward<Args>(args)...), normal_(std::move(n)), scvfIndex_(index) {}

    template<class... Args>
    FaceInterpolationPointData(const GlobalPosition& n, LocalIndex index, Args&&... args)
    : BaseClass(std::forward<Args>(args)...), normal_(n), scvfIndex_(index) {}

    //! The unit outer normal vector at the quadrature point
    const GlobalPosition& unitOuterNormal() const
    { return normal_; }

    //! The local index of an scvf
    LocalIndex scvfIndex() const
    { return scvfIndex_; }

private:
    GlobalPosition normal_;
    LocalIndex scvfIndex_;
};

} // end namespace Dumux::CVFE

#endif
