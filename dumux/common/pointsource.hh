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
 * \brief A point source class,
 *        i.e. sources located at a single point in space
 */

#ifndef DUMUX_POINTSOURCE_HH
#define DUMUX_POINTSOURCE_HH

#include <dune/geometry/referenceelements.hh>
#include <dumux/common/boundingboxtree.hh>

namespace Dumux
{

namespace Properties
{

// Property forward declarations
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(PrimaryVariables);

} // end namespace Properties

/*!
 * \ingroup Common
 * \brief A point source base class
 */
template<class TypeTag>
class PointSource
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    static const int dimworld = GridView::dimensionworld;
    typedef typename Dune::FieldVector<Scalar, dimworld> GlobalPosition;

public:
    // Contructor
    PointSource(GlobalPosition pos, PrimaryVariables values)
      : pos_(pos), values_(values) {}

    //! return the source values
    const PrimaryVariables& values() const
    { return values_; }

    //! return the source position
    const GlobalPosition& position() const
    { return pos_; }

    //! Divide values by scalar
    template<class T>
    void divideValues(T&& n)
    { values_ /= std::forward<T>(n); }

    //! Add something to the values
    template<class T>
    void addToValues(T&& n)
    { values_ += std::forward<T>(n); }

private:
    GlobalPosition pos_;
    PrimaryVariables values_;
};


/*!
 * \ingroup Common
 * \brief A helper class calculating an DOF index to point source map
 */
template<class TypeTag>
class PointSourceHelper
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    static const int dim = GridView::dimension;

    typedef typename GridView::Grid::ctype CoordScalar;
    typedef typename Dune::ReferenceElements<CoordScalar, dim> ReferenceElements;
    typedef typename Dune::ReferenceElement<CoordScalar, dim> ReferenceElement;

    typedef PointSource<TypeTag> PointSource;
    typedef BoundingBoxTree<GridView> BoundingBoxTree;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:
    //! calculate a DOF index to point source map from given vector of point sources
    static void computePointSourceMap(const Problem& problem,
                                      const std::shared_ptr<BoundingBoxTree>& boundingBoxTree,
                                      std::vector<PointSource>& sources,
                                      std::map<unsigned int, std::vector<PointSource> >& pointSourceMap)
    {
        for (auto&& source : sources)
        {
            // compute in which elements the point source falls
            std::vector<unsigned int> entities = boundingBoxTree->computeEntityCollisions(source.position());
            // split the source values equally among all concerned entities
            source.divideValues(entities.size());
            // loop over all concernes elements
            for (unsigned int eIdx : entities)
            {
                if(isBox)
                {
                    // check in which subcontrolvolume(s) we are
                    auto element = boundingBoxTree->entity(eIdx);
                    FVElementGeometry fvGeometry;
                    fvGeometry.update(problem.gridView(), element);
                    auto globalPos = source.position();
                    // loop over all sub control volumes and check if the point source is inside
                    std::vector<unsigned int> vertices;
                    for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
                    {
                        auto geometry = fvGeometry.subContVolGeometries[scvIdx];
                        const ReferenceElement &refElement = ReferenceElements::general(geometry.type());
                        if (refElement.checkInside(geometry.local(globalPos)))
                            vertices.push_back(problem.model().dofMapper().subIndex(element, scvIdx, dofCodim));
                    }
                    auto sourceValues = source.values();
                    sourceValues /= vertices.size();
                    for (unsigned int vIdx : vertices)
                    {
                        // add the pointsource to the DOF map
                        if (pointSourceMap.count(vIdx))
                            pointSourceMap.at(vIdx).push_back(PointSource(source.position(), sourceValues));
                        else
                            pointSourceMap.insert({vIdx, {PointSource(source.position(), sourceValues)}});
                    }
                }
                else
                {
                    // add the pointsource to the DOF map
                    if (pointSourceMap.count(eIdx))
                        pointSourceMap.at(eIdx).push_back(PointSource(source.position(), source.values()));
                    else
                        pointSourceMap.insert({eIdx, {PointSource(source.position(), source.values())}});
                }
            }
        }
    }
};

} // end namespace Dumux

#endif
