// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEDiscretization
 * \brief Classes representing dofs on elements for control-volume finite element schemes
 */
#ifndef DUMUX_CVFE_LOCAL_DOF_HH
#define DUMUX_CVFE_LOCAL_DOF_HH

#include <dumux/common/indextraits.hh>

namespace Dumux::CVFE {

/*!
 * \ingroup CVFEDiscretization
 * \brief A local degree of freedom from an element perspective
 */
template<class LocalIndex>
class LocalDof
{
public:
    LocalDof(LocalIndex index) : index_(index) {}
    LocalIndex index() const { return index_; }
private:
    LocalIndex index_;
};

/*!
 * \ingroup CVFEDiscretization
 * \brief A local degree of freedom associated with a (sub-)control volume from an element perspective
 */
template<class FVElementGeometry>
class FVLocalDof : public LocalDof<typename IndexTraits<typename FVElementGeometry::GridGeometry::GridView>::LocalIndex>
{
    using LocalIndex = typename IndexTraits<typename FVElementGeometry::GridGeometry::GridView>::LocalIndex;
    using ParentType = LocalDof<LocalIndex>;
public:
    FVLocalDof(LocalIndex dofIndex, const FVElementGeometry& fvGeometry)
    : ParentType(dofIndex), fvGeometry_(fvGeometry) {}

    const typename FVElementGeometry::SubControlVolume& scv() const
    { return fvGeometry_.scv(this->index()); }
private:
    const FVElementGeometry& fvGeometry_;
};

} // end namespace Dumux::CVFE

#endif
