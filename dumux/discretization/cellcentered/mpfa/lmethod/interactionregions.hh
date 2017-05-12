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
    using MpfaHelper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using Element = typename GridView::template Codim<0>::Entity;
    using LocalIndexType = typename InteractionVolume::Traits::LocalIndexType;
    using GlobalIndexType = typename InteractionVolume::Traits::GlobalIndexType;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using LocalBasis = std::array<GlobalPosition, dim>;

public:
    LocalIndexType contiFaceLocalIdx;
    std::vector<GlobalIndexType> scvIndices;
    std::vector<GlobalPosition> scvCenters;
    std::vector<GlobalPosition> normal;
    std::vector<GlobalPosition> nu;
    std::vector<Scalar> detX;
    std::vector<Element> elements;

    std::array<GlobalIndexType, 2> globalScvfIndices;

    // Constructor signature for dim == 2
    template<class ScvSeed, class OuterScvSeed>
    InteractionRegion(const Problem& problem,
                      const FVElementGeometry& fvGeometry,
                      const ScvSeed& scvSeed,
                      const OuterScvSeed& outerSeed1,
                      const OuterScvSeed& outerSeed2,
                      const Element& element1,
                      const Element& element2,
                      const Element& element3)
    {
        contiFaceLocalIdx = scvSeed.contiFaceLocalIdx();
        globalScvfIndices[0] = scvSeed.globalScvfIndices()[contiFaceLocalIdx];
        globalScvfIndices[1] = contiFaceLocalIdx == 0 ? outerSeed1.globalScvfIndex() : outerSeed2.globalScvfIndex();

        elements.resize(3);
        elements[0] = element1;
        elements[1] = element2;
        elements[2] = element3;

        // The participating sub control entities
        const auto& scv1 = fvGeometry.scv(scvSeed.globalIndex());
        const auto& scv2 = fvGeometry.scv(outerSeed1.globalIndex());
        const auto& scv3 = fvGeometry.scv(outerSeed2.globalIndex());
        const auto& scvf1 = fvGeometry.scvf(scvSeed.globalScvfIndices()[0]);
        const auto& scvf2 = fvGeometry.scvf(scvSeed.globalScvfIndices()[1]);

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

        // set up the local basis of the elements
        LocalBasis basis1, basis2, basis3, basis4;
        basis1[0] = f1 - scvCenters[0];
        basis1[1] = f2 - scvCenters[0];
        basis2[0] = v - scvCenters[1];
        basis2[1] = f1 - scvCenters[1];
        basis3[0] = f2 - scvCenters[2];
        basis3[1] = v - scvCenters[2];
        basis4[0] = basis1[0];
        basis4[1] = v - scvCenters[0];

        // calculate nus
        const auto nus1 = MpfaHelper::calculateInnerNormals(basis1);
        const auto nus2 = MpfaHelper::calculateInnerNormals(basis2);
        const auto nus3 = MpfaHelper::calculateInnerNormals(basis3);
        const auto nus4 = MpfaHelper::calculateInnerNormals(basis1);

        nu.resize(7);
        nu[0] = nus1[0];
        nu[1] = nus1[1];
        nu[2] = nus2[0];
        nu[3] = nus2[1];
        nu[4] = nus3[0];
        nu[5] = nus3[1];
        nu[6] = nus4[0];

        detX.resize(3);
        detX[0] = MpfaHelper::calculateDetX(basis1);
        detX[1] = MpfaHelper::calculateDetX(basis2);
        detX[2] = MpfaHelper::calculateDetX(basis3);
    }
};
}; //end namespace

#endif