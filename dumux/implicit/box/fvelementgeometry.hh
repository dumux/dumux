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
 * \brief Represents the finite volume geometry of a single element in
 *        the box scheme.
 */
#ifndef DUMUX_BOX_FV_ELEMENTGEOMETRY_HH
#define DUMUX_BOX_FV_ELEMENTGEOMETRY_HH

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/propertysystem.hh>
#include "properties.hh"

namespace Dumux
{
namespace Properties
{
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(ImplicitUseTwoPointFlux);
}

//! \cond INTERNAL

// Functions to initialize the subcontrol volume data
// structures of BoxFVElementGeometry.
// This is a workaround to be able to use template specialization.
// It is required since it is apparently not possible to
// specialize member functions of template classes in this case.
template <typename BoxFVElementGeometry, int dim>
class _BoxFVElemGeomHelper
{
public:
    template<class GlobalPosition>
    static void fillSubContVolData(BoxFVElementGeometry &fvGeometry,
                                   int numVertices,
                                   GlobalPosition *edgeCoord,
                                   GlobalPosition *faceCoord)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "_BoxFVElemGeomHelper::fillSubContVolData dim = " << dim);
    }
};

template <typename BoxFVElementGeometry>
class _BoxFVElemGeomHelper<BoxFVElementGeometry, 1>
{
public:
    static const int dim = 1;

    template<class GlobalPosition>
    static void fillSubContVolData(BoxFVElementGeometry &fvGeometry,
                                   int numVertices,
                                   GlobalPosition *edgeCoord,
                                   GlobalPosition *faceCoord)
    {
        fvGeometry.subContVol[0].volume = 0.5*fvGeometry.elementVolume;
        fvGeometry.subContVol[1].volume = 0.5*fvGeometry.elementVolume;
    }

    template<class SCVGeometry, class GlobalPosition>
    static void computeGeometries(BoxFVElementGeometry &fvGeometry,
                                   int numVertices,
                                   GlobalPosition *edgeCoord,
                                   GlobalPosition *faceCoord)
    {
        std::vector<GlobalPosition> corners = {fvGeometry.subContVol[0].global, fvGeometry.elementGlobal};
        fvGeometry.subContVolGeometries.push_back(SCVGeometry(Dune::GeometryType(Dune::GeometryType::cube, 1), corners));
        corners = {fvGeometry.elementGlobal, fvGeometry.subContVol[1].global};
        fvGeometry.subContVolGeometries.push_back(SCVGeometry(Dune::GeometryType(Dune::GeometryType::cube, 1), corners));
    }
};

template <typename BoxFVElementGeometry>
class _BoxFVElemGeomHelper<BoxFVElementGeometry, 2>
{
public:
    static const int dim = 2;

    template<class GlobalPosition>
    static void fillSubContVolData(BoxFVElementGeometry &fvGeometry,
                                   int numVertices,
                                   GlobalPosition *edgeCoord,
                                   GlobalPosition *faceCoord)
    {
        switch (numVertices) {
        case 3: // triangle
            for (int scvIdx = 0; scvIdx < fvGeometry.numScv; scvIdx++)
                fvGeometry.subContVol[scvIdx].volume = fvGeometry.elementVolume/3.0;
            break;
        case 4: // quadrilateral
            if (GlobalPosition::dimension == dim) {
                fvGeometry.subContVol[0].volume
                  = fvGeometry.quadrilateralArea(fvGeometry.subContVol[0].global,
                                                 edgeCoord[2],
                                                 fvGeometry.elementGlobal,
                                                 edgeCoord[0]);
                fvGeometry.subContVol[1].volume
                  = fvGeometry.quadrilateralArea(fvGeometry.subContVol[1].global,
                                                 edgeCoord[1],
                                                 fvGeometry.elementGlobal,
                                                 edgeCoord[2]);
                fvGeometry.subContVol[2].volume
                  = fvGeometry.quadrilateralArea(fvGeometry.subContVol[2].global,
                                                 edgeCoord[0],
                                                 fvGeometry.elementGlobal,
                                                 edgeCoord[3]);
                fvGeometry.subContVol[3].volume
                  = fvGeometry.quadrilateralArea(fvGeometry.subContVol[3].global,
                                                 edgeCoord[3],
                                                 fvGeometry.elementGlobal,
                                                 edgeCoord[1]);
            }
            else {
                fvGeometry.subContVol[0].volume
                  = fvGeometry.quadrilateralArea3D(fvGeometry.subContVol[0].global,
                                                   edgeCoord[2],
                                                   fvGeometry.elementGlobal,
                                                   edgeCoord[0]);
                fvGeometry.subContVol[1].volume
                  = fvGeometry.quadrilateralArea3D(fvGeometry.subContVol[1].global,
                                                   edgeCoord[1],
                                                   fvGeometry.elementGlobal,
                                                   edgeCoord[2]);
                fvGeometry.subContVol[2].volume
                  = fvGeometry.quadrilateralArea3D(fvGeometry.subContVol[2].global,
                                                   edgeCoord[0],
                                                   fvGeometry.elementGlobal,
                                                   edgeCoord[3]);
                fvGeometry.subContVol[3].volume
                  = fvGeometry.quadrilateralArea3D(fvGeometry.subContVol[3].global,
                                                   edgeCoord[3],
                                                   fvGeometry.elementGlobal,
                                                   edgeCoord[1]);
            }
            break;
        default:
            DUNE_THROW(Dune::NotImplemented,
                       "_BoxFVElemGeomHelper::fillSubContVolData dim = "
                       << dim << ", numVertices = " << numVertices);
        }
    }

    template<class SCVGeometry, class GlobalPosition>
    static void computeGeometries(BoxFVElementGeometry &fvGeometry,
                                   int numVertices,
                                   GlobalPosition *edgeCoord,
                                   GlobalPosition *faceCoord)
    {
        switch (numVertices)
        {
        case 3: // element is triangle
            {
            std::vector<GlobalPosition> corners = {{fvGeometry.subContVol[0].global, edgeCoord[0], edgeCoord[1], fvGeometry.elementGlobal}};
            fvGeometry.subContVolGeometries.push_back(SCVGeometry(Dune::GeometryType(Dune::GeometryType::cube, 2), corners));
            corners = {{edgeCoord[0], fvGeometry.subContVol[1].global, fvGeometry.elementGlobal, edgeCoord[2]}};
            fvGeometry.subContVolGeometries.push_back(SCVGeometry(Dune::GeometryType(Dune::GeometryType::cube, 2), corners));
            corners = {{edgeCoord[1], fvGeometry.elementGlobal, fvGeometry.subContVol[2].global, edgeCoord[2]}};
            fvGeometry.subContVolGeometries.push_back(SCVGeometry(Dune::GeometryType(Dune::GeometryType::cube, 2), corners));
            break;
            }
        case 4: // element is quadrilateral
            {
            std::vector<GlobalPosition> corners = {{fvGeometry.subContVol[0].global, edgeCoord[2], edgeCoord[0], fvGeometry.elementGlobal}};
            fvGeometry.subContVolGeometries.push_back(SCVGeometry(Dune::GeometryType(Dune::GeometryType::cube, 2), corners));
            corners = {{edgeCoord[2], fvGeometry.subContVol[1].global, fvGeometry.elementGlobal, edgeCoord[1]}};
            fvGeometry.subContVolGeometries.push_back(SCVGeometry(Dune::GeometryType(Dune::GeometryType::cube, 2), corners));
            corners = {{edgeCoord[0], fvGeometry.elementGlobal, fvGeometry.subContVol[2].global, edgeCoord[3]}};
            fvGeometry.subContVolGeometries.push_back(SCVGeometry(Dune::GeometryType(Dune::GeometryType::cube, 2), corners));
            corners = {{fvGeometry.elementGlobal, edgeCoord[1], edgeCoord[3], fvGeometry.subContVol[3].global}};
            fvGeometry.subContVolGeometries.push_back(SCVGeometry(Dune::GeometryType(Dune::GeometryType::cube, 2), corners));
            break;
            }
        default:
            DUNE_THROW(Dune::NotImplemented,
                       "_BoxFVElemGeomHelper::computeGeometries dim = "
                       << dim << ", numVertices = " << numVertices);
        }
    }
};

template <typename BoxFVElementGeometry>
class _BoxFVElemGeomHelper<BoxFVElementGeometry, 3>
{
public:
    static const int dim = 3;

    template<class GlobalPosition>
    static void fillSubContVolData(BoxFVElementGeometry &fvGeometry,
                                   int numVertices,
                                   GlobalPosition *edgeCoord,
                                   GlobalPosition *faceCoord)
    {
        switch (numVertices) {
        case 4: // tetrahedron
            for (int scvIdx = 0; scvIdx < fvGeometry.numScv; scvIdx++)
                fvGeometry.subContVol[scvIdx].volume = fvGeometry.elementVolume/4.0;
            break;
        case 5: // pyramid
            fvGeometry.subContVol[0].volume
              = fvGeometry.hexahedronVolume(fvGeometry.subContVol[0].global,
                                            edgeCoord[2],
                                            faceCoord[0],
                                            edgeCoord[0],
                                            edgeCoord[4],
                                            faceCoord[3],
                                            fvGeometry.elementGlobal,
                                            faceCoord[1]);
              fvGeometry.subContVol[1].volume
                = fvGeometry.hexahedronVolume(fvGeometry.subContVol[1].global,
                                              edgeCoord[1],
                                              faceCoord[0],
                                              edgeCoord[2],
                                              edgeCoord[5],
                                              faceCoord[2],
                                              fvGeometry.elementGlobal,
                                              faceCoord[3]);
            fvGeometry.subContVol[2].volume
              = fvGeometry.hexahedronVolume(fvGeometry.subContVol[2].global,
                                            edgeCoord[0],
                                            faceCoord[0],
                                            edgeCoord[3],
                                            edgeCoord[6],
                                            faceCoord[1],
                                            fvGeometry.elementGlobal,
                                            faceCoord[4]);
            fvGeometry.subContVol[3].volume
              = fvGeometry.hexahedronVolume(fvGeometry.subContVol[3].global,
                                            edgeCoord[3],
                                            faceCoord[0],
                                            edgeCoord[1],
                                            edgeCoord[7],
                                            faceCoord[4],
                                            fvGeometry.elementGlobal,
                                            faceCoord[2]);
            fvGeometry.subContVol[4].volume
              = fvGeometry.elementVolume
              - fvGeometry.subContVol[0].volume
              - fvGeometry.subContVol[1].volume
              - fvGeometry.subContVol[2].volume
              - fvGeometry.subContVol[3].volume;
            break;
        case 6: // prism
            fvGeometry.subContVol[0].volume
              = fvGeometry.hexahedronVolume(fvGeometry.subContVol[0].global,
                                            edgeCoord[3],
                                            faceCoord[3],
                                            edgeCoord[4],
                                            edgeCoord[0],
                                            faceCoord[0],
                                            fvGeometry.elementGlobal,
                                            faceCoord[1]);
            fvGeometry.subContVol[1].volume
              = fvGeometry.hexahedronVolume(fvGeometry.subContVol[1].global,
                                            edgeCoord[5],
                                            faceCoord[3],
                                            edgeCoord[3],
                                            edgeCoord[1],
                                            faceCoord[2],
                                            fvGeometry.elementGlobal,
                                            faceCoord[0]);
            fvGeometry.subContVol[2].volume
              = fvGeometry.hexahedronVolume(fvGeometry.subContVol[2].global,
                                            edgeCoord[4],
                                            faceCoord[3],
                                            edgeCoord[5],
                                            edgeCoord[2],
                                            faceCoord[1],
                                            fvGeometry.elementGlobal,
                                            faceCoord[2]);
            fvGeometry.subContVol[3].volume
              = fvGeometry.hexahedronVolume(edgeCoord[0],
                                            faceCoord[0],
                                            fvGeometry.elementGlobal,
                                            faceCoord[1],
                                            fvGeometry.subContVol[3].global,
                                            edgeCoord[6],
                                            faceCoord[4],
                                            edgeCoord[7]);
            fvGeometry.subContVol[4].volume
              = fvGeometry.hexahedronVolume(edgeCoord[1],
                                            faceCoord[2],
                                            fvGeometry.elementGlobal,
                                            faceCoord[0],
                                            fvGeometry.subContVol[4].global,
                                            edgeCoord[8],
                                            faceCoord[4],
                                            edgeCoord[6]);
            fvGeometry.subContVol[5].volume
              = fvGeometry.hexahedronVolume(edgeCoord[2],
                                            faceCoord[1],
                                            fvGeometry.elementGlobal,
                                            faceCoord[2],
                                            fvGeometry.subContVol[5].global,
                                            edgeCoord[7],
                                            faceCoord[4],
                                            edgeCoord[8]);
            break;
        case 8: // hexahedron
            fvGeometry.subContVol[0].volume
              = fvGeometry.hexahedronVolume(fvGeometry.subContVol[0].global,
                                            edgeCoord[6],
                                            faceCoord[4],
                                            edgeCoord[4],
                                            edgeCoord[0],
                                            faceCoord[2],
                                            fvGeometry.elementGlobal,
                                            faceCoord[0]);
            fvGeometry.subContVol[1].volume
              = fvGeometry.hexahedronVolume(fvGeometry.subContVol[1].global,
                                            edgeCoord[5],
                                            faceCoord[4],
                                            edgeCoord[6],
                                            edgeCoord[1],
                                            faceCoord[1],
                                            fvGeometry.elementGlobal,
                                            faceCoord[2]);
            fvGeometry.subContVol[2].volume
              = fvGeometry.hexahedronVolume(fvGeometry.subContVol[2].global,
                                            edgeCoord[4],
                                            faceCoord[4],
                                            edgeCoord[7],
                                            edgeCoord[2],
                                            faceCoord[0],
                                            fvGeometry.elementGlobal,
                                            faceCoord[3]);
            fvGeometry.subContVol[3].volume
              = fvGeometry.hexahedronVolume(fvGeometry.subContVol[3].global,
                                            edgeCoord[7],
                                            faceCoord[4],
                                            edgeCoord[5],
                                            edgeCoord[3],
                                            faceCoord[3],
                                            fvGeometry.elementGlobal,
                                            faceCoord[1]);
            fvGeometry.subContVol[4].volume
              = fvGeometry.hexahedronVolume(edgeCoord[0],
                                            faceCoord[2],
                                            fvGeometry.elementGlobal,
                                            faceCoord[0],
                                            fvGeometry.subContVol[4].global,
                                            edgeCoord[10],
                                            faceCoord[5],
                                            edgeCoord[8]);
            fvGeometry.subContVol[5].volume
              = fvGeometry.hexahedronVolume(edgeCoord[1],
                                            faceCoord[1],
                                            fvGeometry.elementGlobal,
                                            faceCoord[2],
                                            fvGeometry.subContVol[5].global,
                                            edgeCoord[9],
                                            faceCoord[5],
                                            edgeCoord[10]);
            fvGeometry.subContVol[6].volume
              = fvGeometry.hexahedronVolume(edgeCoord[2],
                                            faceCoord[0],
                                            fvGeometry.elementGlobal,
                                            faceCoord[3],
                                            fvGeometry.subContVol[6].global,
                                            edgeCoord[8],
                                            faceCoord[5],
                                            edgeCoord[11]);
            fvGeometry.subContVol[7].volume
              = fvGeometry.hexahedronVolume(edgeCoord[3],
                                            faceCoord[3],
                                            fvGeometry.elementGlobal,
                                            faceCoord[1],
                                            fvGeometry.subContVol[7].global,
                                            edgeCoord[11],
                                            faceCoord[5],
                                            edgeCoord[9]);
            break;
        default:
            DUNE_THROW(Dune::NotImplemented,
                       "_BoxFVElemGeomHelper::fillSubContVolData dim = "
                       << dim << ", numVertices = " << numVertices);
        }
    }

    template<class SCVGeometry, class GlobalPosition>
    static void computeGeometries(BoxFVElementGeometry &fvGeometry,
                                   int numVertices,
                                   GlobalPosition *edgeCoord,
                                   GlobalPosition *faceCoord)
    {
        switch (numVertices)
        {
        case 4: // element is tetrahedron
            {
            Dune::GeometryType type; type.makeHexahedron();

            std::vector<GlobalPosition>
            corners = {{fvGeometry.subContVol[0].global, edgeCoord[0],
                        edgeCoord[1], faceCoord[0],
                        edgeCoord[3], faceCoord[1],
                        faceCoord[2], fvGeometry.elementGlobal}};
            fvGeometry.subContVolGeometries.push_back(SCVGeometry(type, corners));

            corners = {{edgeCoord[0], fvGeometry.subContVol[1].global,
                        faceCoord[0], edgeCoord[2],
                        faceCoord[1], edgeCoord[4],
                        fvGeometry.elementGlobal, faceCoord[3]}};
            fvGeometry.subContVolGeometries.push_back(SCVGeometry(type, corners));

            corners = {{edgeCoord[1], faceCoord[0],
                        fvGeometry.subContVol[2].global, edgeCoord[2],
                        faceCoord[2], fvGeometry.elementGlobal,
                        edgeCoord[5], faceCoord[3]}};
            fvGeometry.subContVolGeometries.push_back(SCVGeometry(type, corners));

            corners = {{edgeCoord[3], faceCoord[1],
                        faceCoord[2], fvGeometry.elementGlobal,
                        fvGeometry.subContVol[3].global, edgeCoord[4],
                        edgeCoord[5], faceCoord[3]}};
            fvGeometry.subContVolGeometries.push_back(SCVGeometry(type, corners));
            break;
            }
        case 6: // element is prism
            {
            Dune::GeometryType type; type.makeHexahedron();

            std::vector<GlobalPosition>
            corners = {{fvGeometry.subContVol[0].global, edgeCoord[3],
                        edgeCoord[4], faceCoord[3],
                        edgeCoord[0], faceCoord[0],
                        faceCoord[1], fvGeometry.elementGlobal}};
            fvGeometry.subContVolGeometries.push_back(SCVGeometry(type, corners));

            corners = {{edgeCoord[3], fvGeometry.subContVol[1].global,
                        faceCoord[3], edgeCoord[5],
                        faceCoord[0], edgeCoord[1],
                        fvGeometry.elementGlobal, faceCoord[2]}};
            fvGeometry.subContVolGeometries.push_back(SCVGeometry(type, corners));

            corners = {{edgeCoord[4], faceCoord[3],
                        fvGeometry.subContVol[2].global, edgeCoord[5],
                        faceCoord[1], fvGeometry.elementGlobal,
                        edgeCoord[2], faceCoord[2]}};
            fvGeometry.subContVolGeometries.push_back(SCVGeometry(type, corners));

            corners = {{edgeCoord[0], faceCoord[0],
                        faceCoord[1], fvGeometry.elementGlobal,
                        fvGeometry.subContVol[3].global, edgeCoord[6],
                        edgeCoord[7], faceCoord[4]}};
            fvGeometry.subContVolGeometries.push_back(SCVGeometry(type, corners));

            corners = {{faceCoord[0], edgeCoord[1],
                        fvGeometry.elementGlobal, faceCoord[2],
                        edgeCoord[6], fvGeometry.subContVol[4].global,
                        faceCoord[4], edgeCoord[8]}};
            fvGeometry.subContVolGeometries.push_back(SCVGeometry(type, corners));

            corners = {{faceCoord[1], fvGeometry.elementGlobal,
                        edgeCoord[2], faceCoord[2],
                        edgeCoord[7], faceCoord[4],
                        fvGeometry.subContVol[5].global, edgeCoord[8]}};
            fvGeometry.subContVolGeometries.push_back(SCVGeometry(type, corners));
            break;
            }
        case 8: // element is hexahedron
            {
            Dune::GeometryType type; type.makeHexahedron();

            std::vector<GlobalPosition>
            corners = {{fvGeometry.subContVol[0].global, edgeCoord[6],
                        edgeCoord[4], faceCoord[4],
                        edgeCoord[0], faceCoord[2],
                        faceCoord[0], fvGeometry.elementGlobal}};
            fvGeometry.subContVolGeometries.push_back(SCVGeometry(type, corners));

            corners = {{edgeCoord[6], fvGeometry.subContVol[1].global,
                        faceCoord[4], edgeCoord[5],
                        faceCoord[2], edgeCoord[1],
                        fvGeometry.elementGlobal, faceCoord[1]}};
            fvGeometry.subContVolGeometries.push_back(SCVGeometry(type, corners));

            corners = {{edgeCoord[4], faceCoord[4],
                        fvGeometry.subContVol[2].global, edgeCoord[7],
                        faceCoord[0], fvGeometry.elementGlobal,
                        edgeCoord[2], faceCoord[3]}};
            fvGeometry.subContVolGeometries.push_back(SCVGeometry(type, corners));

            corners = {{faceCoord[4], edgeCoord[5],
                        edgeCoord[7], fvGeometry.subContVol[3].global,
                        fvGeometry.elementGlobal, faceCoord[1],
                        faceCoord[3], edgeCoord[3]}};
            fvGeometry.subContVolGeometries.push_back(SCVGeometry(type, corners));

            corners = {{edgeCoord[0], faceCoord[2],
                        faceCoord[0], fvGeometry.elementGlobal,
                        fvGeometry.subContVol[4].global, edgeCoord[10],
                        edgeCoord[8], faceCoord[5]}};
            fvGeometry.subContVolGeometries.push_back(SCVGeometry(type, corners));

            corners = {{faceCoord[2], edgeCoord[1],
                        fvGeometry.elementGlobal, faceCoord[1],
                        edgeCoord[10], fvGeometry.subContVol[5].global,
                        faceCoord[5], edgeCoord[9]}};
            fvGeometry.subContVolGeometries.push_back(SCVGeometry(type, corners));

            corners = {{faceCoord[0], fvGeometry.elementGlobal,
                        edgeCoord[2], faceCoord[3],
                        edgeCoord[8], faceCoord[5],
                        fvGeometry.subContVol[6].global, edgeCoord[11]}};
            fvGeometry.subContVolGeometries.push_back(SCVGeometry(type, corners));

            corners = {{fvGeometry.elementGlobal,faceCoord[1],
                        faceCoord[3], edgeCoord[3],
                        faceCoord[5], edgeCoord[9],
                        edgeCoord[11], fvGeometry.subContVol[7].global}};
            fvGeometry.subContVolGeometries.push_back(SCVGeometry(type, corners));
            break;
            }
        default:
            DUNE_THROW(Dune::NotImplemented,
                       "_BoxFVElemGeomHelper::computeGeometries dim = "
                       << dim << ", numVertices = " << numVertices);
        }
    }
};

//! \endcond

/*!
 * \ingroup BoxModel
 * \brief Represents the finite volume geometry of a single element in
 *        the box scheme.
 *
 * The box scheme is a vertex centered finite volume approach. This
 * means that each vertex corrosponds to a control volume which
 * intersects each of the vertex' neighboring elements. If only
 * looking at a single element of the primary grid (which is what this
 * class does), the element is subdivided into multiple fragments of
 * control volumes called sub-control volumes. Each of the element's
 * vertices corrosponds to exactly one sub-control volume in this
 * scenario.
 *
 * For the box methods the sub-control volumes are constructed by
 * connecting the element's center with each edge of the element.
 */
template<class TypeTag>
class BoxFVElementGeometry
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum{dim = GridView::dimension};
    enum{dimWorld = GridView::dimensionworld};

    typedef BoxFVElementGeometry<TypeTag>   ThisType;

    /** \todo Please doc me! */
    friend class _BoxFVElemGeomHelper<ThisType, dim>;

    typedef _BoxFVElemGeomHelper<ThisType, dim> BoxFVElemGeomHelper;

    enum{maxNC = (dim < 3 ? 4 : 8)};
    enum{maxNE = (dim < 3 ? 4 : 12)};
    enum{maxNF = (dim < 3 ? 1 : 6)};
    enum{maxCOS = (dim < 3 ? 2 : 4)};
    enum{maxBF = (dim < 3 ? 8 : 24)};
    enum{maxNFAP = (dim < 3 ? 4 : 8)};
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename Element::Geometry Geometry;
    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldVector<CoordScalar,dim> LocalPosition;
    typedef typename Geometry::JacobianInverseTransposed JacobianInverseTransposed;
    typedef typename Dune::ReferenceElements<CoordScalar, dim> ReferenceElements;
    typedef typename Dune::ReferenceElement<CoordScalar, dim> ReferenceElement;
    typedef typename Dune::MultiLinearGeometry<CoordScalar, dim, dimWorld> BoxGeometry;

    typedef Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1> LocalFiniteElementCache;
    typedef typename LocalFiniteElementCache::FiniteElementType LocalFiniteElement;
    typedef typename LocalFiniteElement::Traits::LocalBasisType::Traits LocalBasisTraits;
    typedef typename LocalBasisTraits::JacobianType ShapeJacobian;

    Scalar quadrilateralArea(const GlobalPosition& p0, const GlobalPosition& p1,
                             const GlobalPosition& p2, const GlobalPosition& p3)
    {
        return 0.5*std::abs((p3[0] - p1[0])*(p2[1] - p0[1]) - (p3[1] - p1[1])*(p2[0] - p0[0]));
    }

    void crossProduct(GlobalPosition& c, const GlobalPosition& a, const GlobalPosition& b)
    {
        c[0] = a[1]*b[2] - a[2]*b[1];
        c[1] = a[2]*b[0] - a[0]*b[2];
        c[2] = a[0]*b[1] - a[1]*b[0];
    }

    Scalar pyramidVolume (const GlobalPosition& p0, const GlobalPosition& p1,
                          const GlobalPosition& p2, const GlobalPosition& p3,
                          const GlobalPosition& p4)
    {
        GlobalPosition a = p2 - p0;
        GlobalPosition b = p3 - p1;
        GlobalPosition n;
        crossProduct(n, a, b);

        a = p4 - p0;

        return 1.0/6.0*(n*a);
    }

    Scalar prismVolume (const GlobalPosition& p0, const GlobalPosition& p1,
                        const GlobalPosition& p2, const GlobalPosition& p3,
                        const GlobalPosition& p4, const GlobalPosition& p5)
    {
        GlobalPosition a = p4 - p0;
        GlobalPosition b = p1 - p3;
        GlobalPosition m;
        crossProduct(m, a, b);

        a = p1 - p0;
        b = p2 - p0;
        GlobalPosition n;
        crossProduct(n, a, b);

        n += m;

        a = p5 - p0;

        return std::abs(1.0/6.0*(n*a));
    }

    Scalar hexahedronVolume (const GlobalPosition& p0, const GlobalPosition& p1,
                             const GlobalPosition& p2, const GlobalPosition& p3,
                             const GlobalPosition& p4, const GlobalPosition& p5,
                             const GlobalPosition& p6, const GlobalPosition& p7)
    {
        return
            prismVolume(p0,p1,p2,p4,p5,p6)
            + prismVolume(p0,p2,p3,p4,p6,p7);
    }

    void normalOfQuadrilateral3D(GlobalPosition &normal,
                                 const GlobalPosition& p0, const GlobalPosition& p1,
                                 const GlobalPosition& p2, const GlobalPosition& p3)
    {
        GlobalPosition a = p2 - p0;
        GlobalPosition b = p3 - p1;
        crossProduct(normal, a, b);
        normal *= 0.5;
    }

    Scalar quadrilateralArea3D(const GlobalPosition& p0, const GlobalPosition& p1,
                               const GlobalPosition& p2, const GlobalPosition& p3)
    {
        GlobalPosition normal;
        normalOfQuadrilateral3D(normal, p0, p1, p2, p3);
        return normal.two_norm();
    }

    void getFaceIndices(int numVertices, int scvfIdx, int& leftFace, int& rightFace)
    {
        static const int edgeToFaceTet[2][6] = {
            {1, 0, 3, 2, 1, 3},
            {0, 2, 0, 1, 3, 2}
        };
        static const int edgeToFacePyramid[2][8] = {
            {0, 2, 3, 0, 1, 3, 4, 2},
            {1, 0, 0, 4, 3, 2, 1, 4}
        };
        static const int edgeToFacePrism[2][9] = {
            {1, 0, 2, 0, 3, 2, 4, 1, 4},
            {0, 2, 1, 3, 1, 3, 0, 4, 2}
        };
        static const int edgeToFaceHex[2][12] = {
            {0, 2, 3, 1, 4, 1, 2, 4, 0, 5, 5, 3},
            {2, 1, 0, 3, 0, 4, 4, 3, 5, 1, 2, 5}
        };

        switch (numVertices) {
        case 4:
            leftFace = edgeToFaceTet[0][scvfIdx];
            rightFace = edgeToFaceTet[1][scvfIdx];
            break;
        case 5:
            leftFace = edgeToFacePyramid[0][scvfIdx];
            rightFace = edgeToFacePyramid[1][scvfIdx];
            break;
        case 6:
            leftFace = edgeToFacePrism[0][scvfIdx];
            rightFace = edgeToFacePrism[1][scvfIdx];
            break;
        case 8:
            leftFace = edgeToFaceHex[0][scvfIdx];
            rightFace = edgeToFaceHex[1][scvfIdx];
            break;
        default:
            DUNE_THROW(Dune::NotImplemented,
                       "BoxFVElementGeometry :: getFaceIndices for numVertices = "
                       << numVertices);
            break;
        }
    }

    void getEdgeIndices(int numVertices, int scvfIdx, int vIdx, int& leftEdge, int& rightEdge)
    {
        static const int faceAndVertexToLeftEdgeTet[4][4] = {
            { 0, 0, 2, -1},
            { 0, 0, -1, 3},
            { 1, -1, 1, 3},
            {-1, 2, 2, 4}
        };
        static const int faceAndVertexToRightEdgeTet[4][4] = {
            { 1, 2, 1, -1},
            { 3, 4, -1, 4},
            { 3, -1, 5, 5},
            {-1, 4, 5, 5}
        };
        static const int faceAndVertexToLeftEdgePyramid[5][5] = {
            { 0, 2, 3, 1, -1},
            { 0, -1, 0, -1, 4},
            {-1, 1, -1, 1, 5},
            { 2, 2, -1, -1, 4},
            {-1, -1, 3, 3, 7}
        };
        static const int faceAndVertexToRightEdgePyramid[5][5] = {
            { 2, 1, 0, 3, -1},
            { 4, -1, 6, -1, 6},
            {-1, 5, -1, 7, 7},
            { 4, 5, -1, -1, 5},
            {-1, -1, 6, 7, 6}
        };
        static const int faceAndVertexToLeftEdgePrism[5][6] = {
            { 3, 3, -1, 0, 1, -1},
            { 4, -1, 4, 0, -1, 2},
            {-1, 5, 5, -1, 1, 2},
            { 3, 3, 5, -1, -1, -1},
            {-1, -1, -1, 6, 6, 8}
        };
        static const int faceAndVertexToRightEdgePrism[5][6] = {
            { 0, 1, -1, 6, 6, -1},
            { 0, -1, 2, 7, -1, 7},
            {-1, 1, 2, -1, 8, 8},
            { 4, 5, 4, -1, -1, -1},
            {-1, -1, -1, 7, 8, 7}
        };
        static const int faceAndVertexToLeftEdgeHex[6][8] = {
            { 0, -1, 4, -1, 8, -1, 2, -1},
            {-1, 5, -1, 3, -1, 1, -1, 9},
            { 6, 1, -1, -1, 0, 10, -1, -1},
            {-1, -1, 2, 7, -1, -1, 11, 3},
            { 4, 6, 7, 5, -1, -1, -1, -1},
            {-1, -1, -1, -1, 10, 9, 8, 11}
        };
        static const int faceAndVertexToRightEdgeHex[6][8] = {
            { 4, -1, 2, -1, 0, -1, 8, -1},
            {-1, 1, -1, 5, -1, 9, -1, 3},
            { 0, 6, -1, -1, 10, 1, -1, -1},
            {-1, -1, 7, 3, -1, -1, 2, 11},
            { 6, 5, 4, 7, -1, -1, -1, -1},
            {-1, -1, -1, -1, 8, 10, 11, 9}
        };

        switch (numVertices) {
        case 4:
            leftEdge = faceAndVertexToLeftEdgeTet[scvfIdx][vIdx];
            rightEdge = faceAndVertexToRightEdgeTet[scvfIdx][vIdx];
            break;
        case 5:
            leftEdge = faceAndVertexToLeftEdgePyramid[scvfIdx][vIdx];
            rightEdge = faceAndVertexToRightEdgePyramid[scvfIdx][vIdx];
            break;
        case 6:
            leftEdge = faceAndVertexToLeftEdgePrism[scvfIdx][vIdx];
            rightEdge = faceAndVertexToRightEdgePrism[scvfIdx][vIdx];
            break;
        case 8:
            leftEdge = faceAndVertexToLeftEdgeHex[scvfIdx][vIdx];
            rightEdge = faceAndVertexToRightEdgeHex[scvfIdx][vIdx];
            break;
        default:
            DUNE_THROW(Dune::NotImplemented,
                       "BoxFVElementGeometry :: getFaceIndices for numVertices = "
                       << numVertices);
            break;
        }
    }

public:
    int boundaryFaceIndex(const int fIdx, const int vIdxInFace) const
    {
        return (fIdx*maxCOS + vIdxInFace);
    }

    struct SubControlVolume //! FV intersected with element
    {
        LocalPosition local; //!< local position
        GlobalPosition global; //!< global position
        LocalPosition localCenter; //!< local position of scv center
        Scalar volume; //!< volume of scv
        Dune::FieldVector<GlobalPosition, maxNC> grad; //! derivative of shape function associated with the sub control volume
        Dune::FieldVector<GlobalPosition, maxNC> gradCenter; //! derivative of shape function at the center of the sub control volume
        Dune::FieldVector<Scalar, maxNC> shapeValue; //! value of shape function associated with the sub control volume
        bool inner;
    };

    struct SubControlVolumeFace //! interior face of a sub control volume
    {
        int i,j; //!< scvf seperates corner i and j of elem
        LocalPosition ipLocal; //!< integration point in local coords
        GlobalPosition ipGlobal; //!< integration point in global coords
        GlobalPosition normal; //!< normal on face pointing to CV j or outward of the domain with length equal to |scvf|
        Scalar area; //!< area of face
        Dune::FieldVector<GlobalPosition, maxNC> grad; //!< derivatives of shape functions at ip
        Dune::FieldVector<Scalar, maxNC> shapeValue; //!< value of shape functions at ip
        Dune::FieldVector<int, maxNFAP> fapIndices; //!< indices w.r.t.neighbors of the flux approximation points
        unsigned numFap; //!< number of flux approximation points
    };

    typedef SubControlVolumeFace BoundaryFace; //!< compatibility typedef

    LocalPosition elementLocal; //!< local coordinate of element center
    GlobalPosition elementGlobal; //!< global coordinate of element center
    Scalar elementVolume; //!< element volume
    SubControlVolume subContVol[maxNC]; //!< data of the sub control volumes
    std::vector<BoxGeometry> subContVolGeometries; //!< geometries of the subcontrol volumes
    SubControlVolumeFace subContVolFace[maxNE]; //!< data of the sub control volume faces
    BoundaryFace boundaryFace[maxBF]; //!< data of the boundary faces
    int numScv; //!< number of subcontrol volumes
    int numScvf; //!< number of inner-domain subcontrolvolume faces
    int numNeighbors; //!< needed for compatibility with cc models
    std::vector<Element> neighbors; //!< needed for compatibility with cc models

    const LocalFiniteElementCache feCache_;

    void updateInner(const Element& element) //!< needed for compatibility with cc models
    {}

    void update(const GridView& gridView, const Element& element)
    {
        const Geometry& geometry = element.geometry();
        Dune::GeometryType geomType = geometry.type();

        const ReferenceElement &referenceElement = ReferenceElements::general(geomType);

        const LocalFiniteElement
            &localFiniteElement = feCache_.get(geomType);

        elementVolume = geometry.volume();
        elementLocal = referenceElement.position(0,0);
        elementGlobal = geometry.global(elementLocal);

        numScv = referenceElement.size(dim);
        numScvf = referenceElement.size(dim-1);
        numNeighbors = 0;

        // subcontrol volumes:
        for (int scvIdx = 0; scvIdx < numScv; scvIdx++) {
            subContVol[scvIdx].local = referenceElement.position(scvIdx, dim);
            subContVol[scvIdx].global = geometry.global(subContVol[scvIdx].local);
            subContVol[scvIdx].inner = true;
        }

        // edges:
        GlobalPosition *edgeCoordinates = new GlobalPosition[numScvf];
        for (int scvfIdx = 0; scvfIdx < numScvf; scvfIdx++) {
            edgeCoordinates[scvfIdx] = geometry.global(referenceElement.position(scvfIdx, dim-1));
        }

        // faces:
        int elementFaces = (dim<3)?0:referenceElement.size(1);
        GlobalPosition *faceCoordinates = new GlobalPosition[elementFaces];
        for (int fIdx = 0; fIdx < elementFaces; fIdx++) {
            faceCoordinates[fIdx] = geometry.global(referenceElement.position(fIdx, 1));
        }

        // fill sub control volume data use specialization for this
        BoxFVElemGeomHelper::fillSubContVolData(*this, numScv, edgeCoordinates, faceCoordinates);
        subContVolGeometries.clear();
        BoxFVElemGeomHelper::template computeGeometries<BoxGeometry>(*this, numScv, edgeCoordinates, faceCoordinates);

        // fill sub control volume face data:
        for (int scvfIdx = 0; scvfIdx < numScvf; scvfIdx++) { // begin loop over edges / sub control volume faces
            SubControlVolumeFace& scvFace = subContVolFace[scvfIdx];

            int i = referenceElement.subEntity(scvfIdx, dim-1, 0, dim);
            int j = referenceElement.subEntity(scvfIdx, dim-1, 1, dim);
            if (numScvf == 4 && (i == 2 || j == 2))
            {
                using std::swap;
                swap(i, j);
            }
            scvFace.i = i;
            scvFace.j = j;

            // calculate the local integration point and
            // the face normal. note that since dim is a
            // constant which is known at compile time
            // the compiler can optimize away all if
            // cases which don't apply.
            LocalPosition ipLocal;
            if (dim==1) {
                ipLocal = 0.5;
                scvFace.ipLocal = ipLocal;
                scvFace.normal = subContVol[j].global - subContVol[i].global;
                scvFace.normal /= scvFace.normal.two_norm();
            }
            else if (dim==2 && dimWorld ==2) {
                ipLocal = referenceElement.position(scvfIdx, dim-1) + elementLocal;
                ipLocal *= 0.5;
                scvFace.ipLocal = ipLocal;
                GlobalPosition diffVec = elementGlobal - edgeCoordinates[scvfIdx];
                scvFace.normal[0] = diffVec[1];
                scvFace.normal[1] = -diffVec[0];

                diffVec = subContVol[j].global - subContVol[i].global;
                // make sure the normal points to the right direction
                if (scvFace.normal * diffVec < 0)
                    scvFace.normal *= -1.0;
            }
            else if (dim==2 && dimWorld==3) {
                ipLocal = referenceElement.position(scvfIdx, dim-1) + elementLocal;
                ipLocal *= 0.5;
                scvFace.ipLocal = ipLocal;
                // normal in 3d world
                const auto faceVec = elementGlobal - edgeCoordinates[scvfIdx];
                const auto elemVec1 = elementGlobal - geometry.corner(0);
                const auto elemVec2 = elementGlobal - geometry.corner(1);
                GlobalPosition elemNormal;
                this->crossProduct(elemNormal, elemVec1, elemVec2);
                this->crossProduct(scvFace.normal, faceVec, elemNormal);

                const auto diffVec = subContVol[j].global - subContVol[i].global;
                // make sure the normal points to the right direction
                if (scvFace.normal * diffVec < 0)
                    scvFace.normal *= -1.0;
            }
            else if (dim==3) {
                int leftFace;
                int rightFace;
                getFaceIndices(numScv, scvfIdx, leftFace, rightFace);
                ipLocal = referenceElement.position(scvfIdx, dim-1) + elementLocal
                    + referenceElement.position(leftFace, 1)
                    + referenceElement.position(rightFace, 1);
                ipLocal *= 0.25;
                scvFace.ipLocal = ipLocal;
                normalOfQuadrilateral3D(scvFace.normal,
                                        edgeCoordinates[scvfIdx], faceCoordinates[rightFace],
                                        elementGlobal, faceCoordinates[leftFace]);
            }

            scvFace.area = scvFace.normal.two_norm();

            bool useTwoPointFlux
                = GET_PARAM_FROM_GROUP(TypeTag, bool, Implicit, UseTwoPointFlux);

            if (useTwoPointFlux)
            {
                GlobalPosition distVec = subContVol[i].global;
                distVec -= subContVol[j].global;
                distVec /= distVec.two_norm2();

                // gradients using a two-point flux approximation
                scvFace.numFap = 2;
                for (unsigned int fapIdx = 0; fapIdx < scvFace.numFap; fapIdx++)
                {
                    scvFace.grad[fapIdx] = distVec;
                    scvFace.shapeValue[fapIdx] = 0.5;
                }
                scvFace.grad[1] *= -1.0;

                scvFace.fapIndices[0] = scvFace.i;
                scvFace.fapIndices[1] = scvFace.j;
            }
            else
            {
                // get the global integration point and the Jacobian inverse
                scvFace.ipGlobal = geometry.global(ipLocal);
                JacobianInverseTransposed jacInvT =
                    geometry.jacobianInverseTransposed(ipLocal);

                // calculate the shape function gradients
                //typedef Dune::FieldVector< Dune::FieldVector< CoordScalar, dim >, 1 > ShapeJacobian;
                typedef Dune::FieldVector< Scalar, 1 > ShapeValue;
                std::vector<ShapeJacobian> localJac;
                std::vector<ShapeValue>    shapeVal;
                localFiniteElement.localBasis().evaluateJacobian(scvFace.ipLocal, localJac);
                localFiniteElement.localBasis().evaluateFunction(scvFace.ipLocal, shapeVal);
                scvFace.numFap = numScv;
                for (int scvIdx = 0; scvIdx < numScv; scvIdx++) {
                    jacInvT.mv(localJac[scvIdx][0], scvFace.grad[scvIdx]);
                    scvFace.shapeValue[scvIdx] = Scalar(shapeVal[scvIdx]);
                    scvFace.fapIndices[scvIdx] = scvIdx;
                }
            }
        } // end loop over edges / sub control volume faces

        // fill boundary face data:
        for (const auto& intersection : intersections(gridView, element))
            if (intersection.boundary())
            {
                int fIdx = intersection.indexInInside();
                int numVerticesOfFace = referenceElement.size(fIdx, 1, dim);
                for (int vIdxInFace = 0; vIdxInFace < numVerticesOfFace; vIdxInFace++)
                {
                    int scvIdx = referenceElement.subEntity(fIdx, 1, vIdxInFace, dim);
                    SubControlVolume& subControlVolume = subContVol[scvIdx];

                    int bfIdx = boundaryFaceIndex(fIdx, vIdxInFace);
                    SubControlVolumeFace& bFace = boundaryFace[bfIdx];

                    subControlVolume.inner = false;
                    switch ((short) dim) {
                    case 1:
                        bFace.ipLocal = referenceElement.position(scvIdx, dim);
                        bFace.area = 1.0;
                        break;
                    case 2:
                        bFace.ipLocal = referenceElement.position(scvIdx, dim)
                                      + referenceElement.position(fIdx, 1);
                        bFace.ipLocal *= 0.5;
                        bFace.area = 0.5*intersection.geometry().volume();
                        break;
                    case 3:
                        int leftEdge;
                        int rightEdge;
                        getEdgeIndices(numScv, fIdx, scvIdx, leftEdge, rightEdge);
                        bFace.ipLocal = referenceElement.position(scvIdx, dim)
                                      + referenceElement.position(fIdx, 1)
                                      + referenceElement.position(leftEdge, dim-1)
                                      + referenceElement.position(rightEdge, dim-1);
                        bFace.ipLocal *= 0.25;
                        bFace.area = quadrilateralArea3D(subControlVolume.global,
                                                         edgeCoordinates[rightEdge],
                                                         faceCoordinates[fIdx],
                                                         edgeCoordinates[leftEdge]);
                        break;
                    default:
                        DUNE_THROW(Dune::NotImplemented, "BoxFVElementGeometry for dim = " << dim);
                    }
                    bFace.ipGlobal = geometry.global(bFace.ipLocal);
                    bFace.i = scvIdx;
                    bFace.j = scvIdx;

                    bFace.normal = intersection.centerUnitOuterNormal();
                    bFace.normal *= bFace.area;

                    typedef Dune::FieldVector< Scalar, 1 > ShapeValue;
                    std::vector<ShapeJacobian> localJac;
                    std::vector<ShapeValue>    shapeVal;
                    localFiniteElement.localBasis().evaluateJacobian(bFace.ipLocal, localJac);
                    localFiniteElement.localBasis().evaluateFunction(bFace.ipLocal, shapeVal);

                    JacobianInverseTransposed jacInvT =
                        geometry.jacobianInverseTransposed(bFace.ipLocal);
                    bFace.numFap = numScv;
                    for (int scvIdx2 = 0; scvIdx2 < numScv; scvIdx2++)
                    {
                        jacInvT.mv(localJac[scvIdx2][0], bFace.grad[scvIdx2]);
                        bFace.shapeValue[scvIdx2] = Scalar(shapeVal[scvIdx2]);
                        bFace.fapIndices[scvIdx2] = scvIdx2;
                    }
                }
            }


        bool evalGradientsAtSCVCenter = GET_PROP_VALUE(TypeTag, EvalGradientsAtSCVCenter);
        if(evalGradientsAtSCVCenter)
        {
            // calculate gradients at the center of the scv
            for (int scvIdx = 0; scvIdx < numScv; scvIdx++){
                SubControlVolume& subControlVolume = subContVol[scvIdx];

                if (dim == 2)
                {
                    switch (scvIdx)
                    {
                    case 0:
                        if (numScv == 4) {
                            subControlVolume.localCenter[0] = 0.25;
                            subControlVolume.localCenter[1] = 0.25;
                        }
                        else {
                            subControlVolume.localCenter[0] = 1.0/6.0;
                            subControlVolume.localCenter[1] = 1.0/6.0;
                        }
                        break;
                    case 1:
                        if (numScv == 4) {
                            subControlVolume.localCenter[0] = 0.75;
                            subControlVolume.localCenter[1] = 0.25;
                        }
                        else {
                            subControlVolume.localCenter[0] = 4.0/6.0;
                            subControlVolume.localCenter[1] = 1.0/6.0;
                        }
                        break;
                    case 2:
                        if (numScv == 4) {
                            subControlVolume.localCenter[0] = 0.25;
                            subControlVolume.localCenter[1] = 0.75;
                        }
                        else {
                            subControlVolume.localCenter[0] = 1.0/6.0;
                            subControlVolume.localCenter[1] = 4.0/6.0;
                        }
                        break;
                    case 3:
                        subControlVolume.localCenter[0] = 0.75;
                        subControlVolume.localCenter[1] = 0.75;
                        break;
                    }
                }

                else if (dim == 3)
                {
                    switch (scvIdx)
                    {
                    case 0:
                        if (numScv == 8) {
                            subControlVolume.localCenter[0] = 0.25;
                            subControlVolume.localCenter[1] = 0.25;
                            subControlVolume.localCenter[2] = 0.25;
                        }
                        else if (numScv == 4) {
                            subControlVolume.localCenter[0] = 3.0/16.0;
                            subControlVolume.localCenter[1] = 3.0/16.0;
                            subControlVolume.localCenter[2] = 3.0/16.0;
                        }
                        break;
                    case 1:
                        if (numScv == 8) {
                            subControlVolume.localCenter[0] = 0.75;
                            subControlVolume.localCenter[1] = 0.25;
                            subControlVolume.localCenter[2] = 0.25;
                        }
                        else if (numScv == 4) {
                            subControlVolume.localCenter[0] = 7.0/16.0;
                            subControlVolume.localCenter[1] = 3.0/16.0;
                            subControlVolume.localCenter[2] = 3.0/16.0;
                        }
                        break;
                    case 2:
                        if (numScv == 8) {
                            subControlVolume.localCenter[0] = 0.25;
                            subControlVolume.localCenter[1] = 0.75;
                            subControlVolume.localCenter[2] = 0.25;
                        }
                        else if (numScv == 4) {
                            subControlVolume.localCenter[0] = 3.0/16.0;
                            subControlVolume.localCenter[1] = 7.0/16.0;
                            subControlVolume.localCenter[2] = 3.0/16.0;
                        }
                        break;
                    case 3:
                        if (numScv == 8) {
                            subControlVolume.localCenter[0] = 0.75;
                            subControlVolume.localCenter[1] = 0.75;
                            subControlVolume.localCenter[2] = 0.25;
                        }
                        else if (numScv == 4) {
                            subControlVolume.localCenter[0] = 3.0/16.0;
                            subControlVolume.localCenter[1] = 3.0/16.0;
                            subControlVolume.localCenter[2] = 7.0/16.0;
                        }
                        break;
                    case 4:
                        subControlVolume.localCenter[0] = 0.25;
                        subControlVolume.localCenter[1] = 0.25;
                        subControlVolume.localCenter[2] = 0.75;
                        break;
                    case 5:
                        subControlVolume.localCenter[0] = 0.75;
                        subControlVolume.localCenter[1] = 0.25;
                        subControlVolume.localCenter[2] = 0.75;
                        break;
                    case 6:
                        subControlVolume.localCenter[0] = 0.25;
                        subControlVolume.localCenter[1] = 0.75;
                        subControlVolume.localCenter[2] = 0.75;
                        break;
                    case 7:
                        subControlVolume.localCenter[0] = 0.75;
                        subControlVolume.localCenter[1] = 0.75;
                        subControlVolume.localCenter[2] = 0.75;
                        break;
                    }
                }
                std::vector<ShapeJacobian> localJac;
                localFiniteElement.localBasis().evaluateJacobian(subControlVolume.localCenter, localJac);

                JacobianInverseTransposed jacInvT =
                    geometry.jacobianInverseTransposed(subControlVolume.localCenter);
                for (int vIdx = 0; vIdx < numScv; vIdx++)
                    jacInvT.mv(localJac[vIdx][0], subControlVolume.gradCenter[vIdx]);
            }
        }

        delete[] edgeCoordinates;
        delete[] faceCoordinates;
    }
};
}

#endif
