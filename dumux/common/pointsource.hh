// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup Common
 * \brief A point source class,
 *        i.e. sources located at a single point in space
 */

#ifndef DUMUX_POINTSOURCE_HH
#define DUMUX_POINTSOURCE_HH

#include <functional>

#include <dune/common/reservedvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/intersectspointgeometry.hh>
#include <dumux/geometry/intersectingentities.hh>

#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup Common
 * \brief A point source base class
 * \tparam PositionType the position type
 * \tparam ValueType the a vector type storing the source for all equations
 */
template<class PositionType, class ValueType>
class PointSource
{
public:
    //! Export the scalar type
    using Scalar = std::decay_t<decltype(std::declval<ValueType>()[0])>;
    //! Export the position type
    using GlobalPosition = PositionType;
    //! Export the value type
    using Values = ValueType;

    //! Constructor for constant point sources
    PointSource(GlobalPosition pos, Values values)
      : values_(values), pos_(pos), embeddings_(1) {}

    //! Constructor for sol dependent point sources, when there is no
    // value known at the time of initialization
    PointSource(GlobalPosition pos)
      : values_(0.0), pos_(pos), embeddings_(1) {}

    //! Convenience += operator overload modifying only the values
    PointSource& operator+= (Scalar s)
    {
        values_ += s;
        return *this;
    }

    //! Convenience -= operator overload modifying only the values
    PointSource& operator-= (Scalar s)
    {
        values_ -= s;
        return *this;
    }

    //! Convenience *= operator overload modifying only the values
    PointSource& operator*= (Scalar s)
    {
        values_ *= s;
        return *this;
    }

    //! Convenience /= operator overload modifying only the values
    PointSource& operator/= (Scalar s)
    {
        values_ /= s;
        return *this;
    }

    //! Convenience = operator overload modifying only the values
    PointSource& operator= (const Values& values)
    {
        values_ = values;
        return *this;
    }

    //! Convenience = operator overload modifying only the values
    PointSource& operator= (Scalar s)
    {
        values_ = s;
        return *this;
    }

    //! return the source values
    Values values() const
    { return values_; }

    //! return the source position
    const GlobalPosition& position() const
    { return pos_; }

    //! an update function called before adding the value
    // to the local residual in the problem in scvPointSources
    // to be overloaded by derived classes
    template<class Problem, class FVElementGeometry, class ElementVolumeVariables>
    void update(const Problem &problem,
                const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity &element,
                const FVElementGeometry &fvGeometry,
                const ElementVolumeVariables &elemVolVars,
                const typename FVElementGeometry::SubControlVolume &scv)
    {}

    //! set the number of embeddings for this point source
    void setEmbeddings(std::size_t embeddings)
    {
        embeddings_ = embeddings;
    }

    /*!
     * \brief get the number of embeddings for this point source
     * \note A point source might be located on the intersection between several scvs.
     *       If so, there are point sources for every neighboring scv with the same position.
     *       `embeddings` returns the number of neighboring scvs.
     * Example: If I want to inject 1kg/s at a location that is on the inner face of an scv
     *          the point source exists in both scvs. Both have a value of 1kg/s.
     *          We then divide the value by the number of embeddings to not inject 2kg/s but 1kg/s.
     * \note This division is done in the problem.scvPointSources() if this behaviour is not explicitly
     *       changed by e.g. overloading this function in the problem implementation.
     */
    std::size_t embeddings() const
    {
        return embeddings_;
    }

protected:
    Values values_; //!< value of the point source for each equation
private:
    GlobalPosition pos_; //!< position of the point source
    std::size_t embeddings_; //!< how many SCVs the point source is associated with
};

/*!
 * \ingroup Common
 * \brief A point source class with an identifier to attach data
 * \tparam GlobalPosition the position type
 * \tparam SourceValues the a vector type storing the source for all equations
 * \tparam I the ID type
 */
template<class GlobalPosition, class SourceValues, class I>
class IdPointSource : public PointSource<GlobalPosition, SourceValues>
{
    using ParentType = PointSource<GlobalPosition, SourceValues>;
    using Scalar = typename ParentType::Scalar;

public:
    //! export the id type
    using IdType = I;

    //! Constructor for constant point sources
    IdPointSource(GlobalPosition pos, SourceValues values, IdType id)
      :  ParentType(pos, values), id_(id) {}

    //! Constructor for sol dependent point sources, when there is no
    // value known at the time of initialization
    IdPointSource(GlobalPosition pos, IdType id)
      : ParentType(pos, SourceValues(0.0)), id_(id) {}

    //! return the sources identifier
    IdType id() const
    { return id_; }

    //! Convenience = operator overload modifying only the values
    IdPointSource& operator= (const SourceValues& values)
    {
        ParentType::operator=(values);
        return *this;
    }

    //! Convenience = operator overload modifying only the values
    IdPointSource& operator= (Scalar s)
    {
        ParentType::operator=(s);
        return *this;
    }

private:
    IdType id_;
};

/*!
 * \ingroup Common
 * \brief A point source class for time dependent point sources
 */
template<class TypeTag>
class SolDependentPointSource : public PointSource<Dune::FieldVector<typename GetPropType<TypeTag, Properties::GridGeometry>::GridView::ctype,
                                                   GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimensionworld>,
                                                   GetPropType<TypeTag, Properties::NumEqVector>>
{
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using SourceValues = GetPropType<TypeTag, Properties::NumEqVector>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;

    static const int dimworld = GridView::dimensionworld;
    using GlobalPosition = typename Dune::FieldVector<typename GridView::ctype, dimworld>;
    // returns the PointSource values as PrimaryVariables
    using ValueFunction = typename std::function<SourceValues(const Problem &problem,
                                                              const Element &element,
                                                              const FVElementGeometry &fvGeometry,
                                                              const ElementVolumeVariables &elemVolVars,
                                                              const SubControlVolume &scv)>;

    using ParentType = PointSource<GlobalPosition, SourceValues>;
    using Scalar = typename ParentType::Scalar;

public:
    //! Constructor for sol dependent point sources, when there is no
    // value known at the time of initialization
    SolDependentPointSource(GlobalPosition pos,
                            ValueFunction valueFunction)
      : ParentType(pos, SourceValues(0.0)), valueFunction_(valueFunction) {}

    //! an update function called before adding the value
    // to the local residual in the problem in scvPointSources
    // to be overloaded by derived classes
    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const ElementVolumeVariables &elemVolVars,
                const SubControlVolume &scv)
    { this->values_ = valueFunction_(problem, element, fvGeometry, elemVolVars, scv); }

    //! Convenience = operator overload modifying only the values
    SolDependentPointSource& operator= (const SourceValues& values)
    {
        ParentType::operator=(values);
        return *this;
    }

    //! Convenience = operator overload modifying only the values
    SolDependentPointSource& operator= (Scalar s)
    {
        ParentType::operator=(s);
        return *this;
    }

private:
    ValueFunction valueFunction_;
};

/*!
 * \ingroup Common
 * \brief A helper class calculating a sub control volume to point source map
 * This class uses the bounding box tree implementation to indentify in which
 * sub control volume(s) a point source falls.
 */
class BoundingBoxTreePointSourceHelper
{
public:
    //! calculate a DOF index to point source map from given vector of point sources
    template<class GridGeometry, class PointSource, class PointSourceMap>
    static void computePointSourceMap(const GridGeometry& gridGeometry,
                                      const std::vector<PointSource>& sources,
                                      PointSourceMap& pointSourceMap,
                                      const std::string& paramGroup = "")
    {
        const auto& boundingBoxTree = gridGeometry.boundingBoxTree();

        for (const auto& s : sources)
        {
            // compute in which elements the point source falls
            const auto entities = intersectingEntities(s.position(), boundingBoxTree);

            // continue with next point source if no intersection with the grid are found
            if (entities.empty())
                continue;

            // make local copy of point source for the map
            auto source = s;

            // split the source values equally among all concerned entities
            source.setEmbeddings(entities.size()*source.embeddings());
            // loop over all concernes elements
            for (const auto eIdx : entities)
            {
                if constexpr (GridGeometry::discMethod == DiscretizationMethod::box)
                {
                    // check in which subcontrolvolume(s) we are
                    const auto element = boundingBoxTree.entitySet().entity(eIdx);
                    auto fvGeometry = localView(gridGeometry);
                    fvGeometry.bindElement(element);

                    const auto globalPos = source.position();
                    // loop over all sub control volumes and check if the point source is inside
                    constexpr int dim = GridGeometry::GridView::dimension;
                    Dune::ReservedVector<std::size_t, 1<<dim> scvIndices;
                    for (const auto& scv : scvs(fvGeometry))
                        if (intersectsPointGeometry(globalPos, scv.geometry()))
                            scvIndices.push_back(scv.indexInElement());

                    // for all scvs that tested positive add the point sources
                    // to the element/scv to point source map
                    for (const auto scvIdx : scvIndices)
                    {
                        const auto key = std::make_pair(eIdx, scvIdx);
                        if (pointSourceMap.count(key))
                            pointSourceMap.at(key).push_back(source);
                        else
                            pointSourceMap.insert({key, {source}});
                        // split equally on the number of matched scvs
                        auto& s = pointSourceMap.at(key).back();
                        s.setEmbeddings(scvIndices.size()*s.embeddings());
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
