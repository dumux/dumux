// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEDiscretization
 * \brief Classes representing integration point data for control-volume finite element schemes
 */
#ifndef DUMUX_CVFE_IP_DATA_HH
#define DUMUX_CVFE_IP_DATA_HH

#include <type_traits>

namespace Dumux::CVFE {

/*!
 * \ingroup CVFEDiscretization
 * \brief An integration point related to an element that includes global and local positions
 */
template<class LocalPosition, class GlobalPosition>
class IntegrationPointData
{
public:
    IntegrationPointData(LocalPosition&& localPos, GlobalPosition&& pos) : ipLocal_(std::move(localPos)), ipGlobal_(std::move(pos)) {}
    IntegrationPointData(const LocalPosition& localPos, const GlobalPosition& pos) :  ipLocal_(localPos), ipGlobal_(pos) {}

    //! The global position of the quadrature point
    const GlobalPosition& ipGlobal() const
    { return ipGlobal_; }

    //! The local position of the quadrature point
    const LocalPosition& ipLocal() const
    { return ipLocal_; }


private:
    LocalPosition ipLocal_;
    GlobalPosition ipGlobal_;
};

/*!
 * \ingroup CVFEDiscretization
 * \brief An integration point related to a localDof of an element, giving its global and local positions
 */
template<class LocalPosition, class GlobalPosition, class LocalIndex>
class LocalDofIntegrationPointData : public IntegrationPointData<LocalPosition, GlobalPosition>
{
    using ParentType = IntegrationPointData<LocalPosition, GlobalPosition>;
public:
    LocalDofIntegrationPointData(LocalPosition&& localPos, GlobalPosition&& pos, LocalIndex index)
    : ParentType(localPos, pos), localDofIndex_(index) {}
    LocalDofIntegrationPointData(const LocalPosition& localPos, const GlobalPosition& pos, LocalIndex index)
    : ParentType(localPos, pos), localDofIndex_(index) {}

    //! The local index of the corresponding dof
    LocalIndex localDofIndex() const
    { return localDofIndex_; }

private:
    LocalIndex localDofIndex_;
};

/*!
 * \ingroup CVFEDiscretization
 * \brief An integration point related to a global position of an element, giving its local positions by a mapping
 */
template<class LocalMapping, class GlobalPosition>
class IntegrationPointDataLocalMapping
{
    using LocalPosition = std::invoke_result_t<LocalMapping, const GlobalPosition&>;

public:
    IntegrationPointDataLocalMapping(LocalMapping&& mapping, GlobalPosition&& pos) : localMapping_(std::move(mapping)), ipGlobal_(std::move(pos)) {}
    IntegrationPointDataLocalMapping(LocalMapping&& mapping, const GlobalPosition& pos) : localMapping_(std::move(mapping)), ipGlobal_(pos) {}

    //! The global position of the quadrature point
    const GlobalPosition& ipGlobal() const
    { return ipGlobal_; }

    //! The local position of the quadrature point
    const LocalPosition ipLocal() const
    { return localMapping_(ipGlobal_); }


private:
    LocalMapping localMapping_;
    GlobalPosition ipGlobal_;
};

/*!
 * \ingroup CVFEDiscretization
 * \brief An integration point related to a face of an element
 */
template<class LocalPosition, class GlobalPosition, class LocalIndex>
class FaceIntegrationPointData : public IntegrationPointData<LocalPosition, GlobalPosition>
{
    using ParentType = IntegrationPointData<LocalPosition, GlobalPosition>;
public:
    FaceIntegrationPointData(GlobalPosition&& localPos, GlobalPosition&& pos, GlobalPosition&& n, LocalIndex index)
    : ParentType(localPos, pos), normal_(std::move(n)), scvfIndex_(index) {}
    FaceIntegrationPointData(const GlobalPosition& localPos, const GlobalPosition& pos, const GlobalPosition& n, LocalIndex index)
    : ParentType(localPos, pos), normal_(n), scvfIndex_(index) {}

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
