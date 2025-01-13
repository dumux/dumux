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
template<class LocalIndex, class GridIndex>
class LocalDof
{
public:
    LocalDof(LocalIndex indexInElement, GridIndex dofIndex, GridIndex eIdx)
    : indexInElement_(indexInElement), dofIndex_(dofIndex), eIdx_(eIdx) {}
    LocalIndex indexInElement() const { return indexInElement_; }
    GridIndex dofIndex() const { return dofIndex_; }
    GridIndex elementIndex() const { return eIdx_; }
private:
    LocalIndex indexInElement_;
    GridIndex dofIndex_;
    GridIndex eIdx_;
};


} // end namespace Dumux::CVFE

namespace Dumux {

//! range over local dofs
template<class FVElementGeometry>
inline auto localDofs(const FVElementGeometry& fvGeometry)
{
    using LocalIndexType = typename IndexTraits<typename FVElementGeometry::GridGeometry::GridView>::LocalIndex;
    using GridIndexType = typename IndexTraits<typename FVElementGeometry::GridGeometry::GridView>::GridIndex;

    return Dune::transformedRangeView(
        Dune::range(fvGeometry.numScv()),
        [&](const auto i) { return CVFE::LocalDof
        {
            static_cast<LocalIndexType>(i),
            static_cast<GridIndexType>(fvGeometry.scv(i).dofIndex()),
            static_cast<GridIndexType>(fvGeometry.scv(i).elementIndex())
        }; }
    );
}

//! range over control-volume local dofs
template<class FVElementGeometry>
inline auto cvLocalDofs(const FVElementGeometry& fvGeometry)
{
    // Default it that all dofs are cv dofs
    return localDofs(fvGeometry);
}

} // end namespace Dumux
#endif
