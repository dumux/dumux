// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CCMpfaDiscretization
 * \brief Helper class providing functionality to compute the geometry
 *        of the interaction-volume local sub-control volumes involved
 *        in the mpfa-o scheme.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_O_SCV_GEOMETRY_HELPER_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_O_SCV_GEOMETRY_HELPER_HH

#include <array>

#include <dune/common/exceptions.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>

namespace Dumux {

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Helper class providing functionality to compute the geometry
 *        of the interaction-volume local sub-control volumes of mpfa-o type.
 */
template<class LocalScvType>
class CCMpfaOScvGeometryHelper
{
    using ctype = typename LocalScvType::ctype;
    using LocalIndexType = typename LocalScvType::LocalIndexType;

    static constexpr int dim = LocalScvType::myDimension;
    static constexpr int dimWorld = LocalScvType::worldDimension;

    struct MLGTraits : public Dune::MultiLinearGeometryTraits<ctype>
    {
        // we know the number of corners is always (2^(dim) corners (1<<dim))
        template< int mydim, int cdim >
        struct CornerStorage
        { using Type = std::array< typename LocalScvType::GlobalCoordinate, (1<<dim) >; };

        // we know all scvs will have the same geometry type
        template< int d >
        struct hasSingleGeometryType
        {
            static const bool v = true;
            static const unsigned int topologyId = Dune::GeometryTypes::cube(d).id();
        };
    };

public:
    //! export the geometry type of the local scvs
    using ScvGeometry = Dune::MultiLinearGeometry<ctype, dim, dimWorld, MLGTraits>;

    //! returns the geometry of the i-th local scv
    template<class InteractionVolume, class FVElementGeometry>
    static ScvGeometry computeScvGeometry(LocalIndexType ivLocalScvIdx,
                                          const InteractionVolume& iv,
                                          const FVElementGeometry& fvGeometry)
    {
        const auto& scv = iv.localScv(ivLocalScvIdx);

        if (dim == 2)
        {
            const auto& firstGridScvf = fvGeometry.scvf(iv.localScvf(scv.localScvfIndex(0)).gridScvfIndex());
            const auto& secondGridScvf = fvGeometry.scvf(iv.localScvf(scv.localScvfIndex(1)).gridScvfIndex());

            typename MLGTraits::template CornerStorage<dim, dimWorld>::Type corners;
            corners[0] = fvGeometry.scv( scv.gridScvIndex() ).center();
            corners[1] = fvGeometry.facetCorner(firstGridScvf);
            corners[2] = fvGeometry.facetCorner(secondGridScvf);
            corners[3] = fvGeometry.vertexCorner(secondGridScvf);

            using std::swap;
            typename LocalScvType::LocalBasis basis;
            basis[0] = corners[1] - corners[0];
            basis[1] = corners[2] - corners[0];
            if ( !fvGeometry.gridGeometry().mpfaHelper().isRightHandSystem(basis) )
                swap(corners[1], corners[2]);

            return ScvGeometry(Dune::GeometryTypes::cube(ScvGeometry::mydimension), corners);
        }
        else if (dim == 3)
            DUNE_THROW(Dune::NotImplemented, "Mpfa-o local scv geometry computation in 3d");
        else
            DUNE_THROW(Dune::InvalidStateException, "Mpfa only works in 2d or 3d");
    }
};

} // end namespace Dumux

#endif
