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

#include <dumux/common/boundingboxtree.hh>

namespace Dumux
{

namespace Properties
{

// Property forward declarations
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(PrimaryVariables);
NEW_PROP_TAG(Problem);
NEW_PROP_TAG(FVElementGeometry);
NEW_PROP_TAG(PointSource);
NEW_PROP_TAG(ImplicitIsBox);

} // end namespace Properties

/*!
 * \ingroup Common
 * \brief A point source base class
 */
template<class TypeTag>
class PointSource : public GET_PROP_TYPE(TypeTag, PrimaryVariables)
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    static const int dimworld = GridView::dimensionworld;
    typedef typename Dune::FieldVector<Scalar, dimworld> GlobalPosition;

public:
    // Constructor
    PointSource(GlobalPosition pos, PrimaryVariables values)
      : PrimaryVariables(values), pos_(pos) {}

    //! return the source values
    const PrimaryVariables& values() const
    { return *this; }

    //! return the source position
    const GlobalPosition& position() const
    { return pos_; }

private:
    GlobalPosition pos_;
};


/*!
 * \ingroup Common
 * \brief A helper class calculating a DOF-index to point source map
 */
template<class TypeTag>
class PointSourceHelper
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PointSource) PointSource;

    static const int dim = GridView::dimension;
    static const int dimworld = GridView::dimensionworld;

    typedef BoundingBoxTree<GridView> BoundingBoxTree;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:
    //! calculate a DOF index to point source map from given vector of point sources
    static void computePointSourceMap(const Problem& problem,
                                      const std::shared_ptr<BoundingBoxTree>& boundingBoxTree,
                                      std::vector<PointSource>& sources,
                                      std::map<std::pair<unsigned int, unsigned int>, std::vector<PointSource> >& pointSourceMap)
    {
        for (auto&& source : sources)
        {
            // compute in which elements the point source falls
            std::vector<unsigned int> entities = boundingBoxTree->computeEntityCollisions(source.position());
            // split the source values equally among all concerned entities
            source /= entities.size();
            // loop over all concernes elements
            for (unsigned int eIdx : entities)
            {
                if(isBox)
                {
                    // check in which subcontrolvolume(s) we are
                    const auto element = boundingBoxTree->entity(eIdx);
                    FVElementGeometry fvGeometry;
                    fvGeometry.update(problem.gridView(), element);
                    const auto globalPos = source.position();
                    // loop over all sub control volumes and check if the point source is inside
                    std::vector<unsigned int> scvs;
                    for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
                    {
                        auto geometry = fvGeometry.subContVolGeometries[scvIdx];
                        if (BoundingBoxTreeHelper<dimworld>::pointInGeometry(geometry, globalPos))
                            scvs.push_back(scvIdx);
                    }
                    // for all scvs that where tested positiv add the point sources
                    // to the element/scv to point source map
                    for (unsigned int scvIdx : scvs)
                    {
                        const auto key = std::make_pair(eIdx, scvIdx);
                        if (pointSourceMap.count(key))
                            pointSourceMap.at(key).push_back(source);
                        else
                            pointSourceMap.insert({key, {source}});
                        // split equally on the number of matched scvs
                        pointSourceMap.at(key).back() /= scvs.size();
                    }
                }
                else
                {
                    // add the pointsource to the DOF map
                    const auto key = std::make_pair(eIdx, /*scvIdx=*/ 0);
                    if (pointSourceMap.count(key))
                        pointSourceMap.at(key).push_back(source);
                    else
                        pointSourceMap.insert({key, {source}});
                }
            }
        }
    }
};

} // end namespace Dumux

#endif
