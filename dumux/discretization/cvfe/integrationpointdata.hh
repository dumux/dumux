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

namespace Dumux::CVFE {

/*!
 * \ingroup CVFEDiscretization
 * \brief An integration point related to an element
 */
template<class GlobalPosition>
class IntegrationPointData
{
public:
    IntegrationPointData(GlobalPosition&& pos) : ipGlobal_(std::move(pos)) {}
    IntegrationPointData(const GlobalPosition& pos) : ipGlobal_(pos) {}

    //! The global position of the quadrature point
    const GlobalPosition& ipGlobal() const
    { return ipGlobal_; }

private:
    GlobalPosition ipGlobal_;
};

/*!
 * \ingroup CVFEDiscretization
 * \brief An integration point related to a face of an element
 */
template<class GlobalPosition>
class FaceIntegrationPointData : public IntegrationPointData<GlobalPosition>
{
    using ParentType = IntegrationPointData<GlobalPosition>;
public:
    FaceIntegrationPointData(GlobalPosition&& pos, GlobalPosition&& n) : ParentType(pos), normal_(std::move(n)) {}
    FaceIntegrationPointData(const GlobalPosition& pos, const GlobalPosition& n) : ParentType(pos), normal_(n) {}

    //! The unit outer normal vector at the quadrature point
    const GlobalPosition& unitOuterNormal() const
    { return normal_; }

private:
    GlobalPosition normal_;
};

/*!
 * \ingroup CVFEDiscretization
 * \brief An integration point related to a face of an element, specific to finite volume schemes
 */
template<class SubControlVolumeFace>
class FVFaceIntegrationPointData : public FaceIntegrationPointData<typename SubControlVolumeFace::GlobalPosition>
{
    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;
    using ParentType = FaceIntegrationPointData<GlobalPosition>;
public:
    FVFaceIntegrationPointData(const auto& scvf)
    : ParentType(scvf.ipGlobal(), scvf.unitOuterNormal()), scvf_(scvf)
    {}

    //! The sub-control volume face
    const SubControlVolumeFace& scvf() const
    { return scvf_; }

    //! Return the boundary flag
    auto boundaryFlag() const
    { return scvf_.boundaryFlag(); }

private:
    const SubControlVolumeFace& scvf_;
};

} // end namespace Dumux::CVFE

#endif
