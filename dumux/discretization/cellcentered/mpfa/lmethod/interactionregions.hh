// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief Base classes for interaction volume of mpfa methods.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_L_INTERACTIONREGIONS_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_L_INTERACTIONREGIONS_HH

namespace Dumux
{
/*!
 * \ingroup Mpfa
 * \brief Base class for the interaction regions of the mpfa-l method.
 */
template<class TypeTag>
struct InteractionRegion
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using Element = typename GridView::template Codim<0>::Entity;
    using LocalIndexType = typename InteractionVolume::LocalIndexType;
    using GlobalIndexType = typename InteractionVolume::GlobalIndexType;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    LocalIndexType contiFaceLocalIdx;
    std::vector<GlobalIndexType> scvIndices;
    std::vector<GlobalPosition> scvCenters;
    std::vector<GlobalPosition> normal;
    std::vector<GlobalPosition> nu;
    std::vector<Scalar> detX;
    std::vector<Element> elements;

    std::array<GlobalIndexType, 2> globalScvfs;

    // Constructor for dim == 2
    template<class ScvSeedType, class OuterScvSeedType>
    InteractionRegion(const Problem& problem,
                      const FVElementGeometry& fvGeometry,
                      const ScvSeedType& scvSeed,
                      const OuterScvSeedType& outerSeed1,
                      const OuterScvSeedType& outerSeed2,
                      const Element& element1,
                      const Element& element2,
                      const Element& element3)
    {
        contiFaceLocalIdx = scvSeed.contiFaceLocalIdx();
        globalScvfs[0] = scvSeed.globalScvfIndices()[contiFaceLocalIdx];
        globalScvfs[1] = contiFaceLocalIdx == 0 ? outerSeed1.globalScvfIndex() : outerSeed2.globalScvfIndex();

        elements.resize(3);
        elements[0] = element1;
        elements[1] = element2;
        elements[2] = element3;

        // The participating sub control entities
        auto&& scv1 = fvGeometry.scv(scvSeed.globalIndex());
        auto&& scv2 = fvGeometry.scv(outerSeed1.globalIndex());
        auto&& scv3 = fvGeometry.scv(outerSeed2.globalIndex());
        auto&& scvf1 = fvGeometry.scvf(scvSeed.globalScvfIndices()[0]);
        auto&& scvf2 = fvGeometry.scvf(scvSeed.globalScvfIndices()[1]);

        // The necessary coordinates and normals
        GlobalPosition v = scvf1.vertexCorner();
        GlobalPosition f1 = scvf1.facetCorner();
        GlobalPosition f2 = scvf2.facetCorner();

        scvCenters.resize(3);
        scvCenters[0] = scv1.center();
        scvCenters[1] = scv2.center();
        scvCenters[2] = scv3.center();

        normal.resize(2);
        normal[0] = scvf1.unitOuterNormal();
        normal[0] *= scvf1.area();
        normal[1] = scvf2.unitOuterNormal();
        normal[1] *= scvf2.area();

        // indices of the participating scvs
        scvIndices.resize(3);
        scvIndices[0] = scv1.index();
        scvIndices[1] = scv2.index();
        scvIndices[2] = scv3.index();

        // calculate nus and detXs
        static const Dune::FieldMatrix<Scalar, dim, dim> R = {{0.0, 1.0}, {-1.0, 0.0}};
        nu.resize(7);
        R.mv(f2-scvCenters[0], nu[0]);
        R.mv(scvCenters[0]-f1, nu[1]);
        R.mv(f1-scvCenters[1], nu[2]);
        R.mv(scvCenters[1]-v, nu[3]);
        R.mv(v-scvCenters[2], nu[4]);
        R.mv(scvCenters[2]-f2, nu[5]);
        R.mv(v-scvCenters[0], nu[6]);

        detX.resize(3);
        detX[0] = (f1-scvCenters[0])*nu[0];
        detX[1] = (v-scvCenters[1])*nu[2];
        detX[2] = (f2-scvCenters[2])*nu[4];
    }
};
}; //end namespace

#endif