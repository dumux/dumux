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
 * \ingroup Common
 * \brief A point source class,
 *        i.e. sources located at a single point in space
 */

#ifndef DUMUX_POINTSOURCE_HH
#define DUMUX_POINTSOURCE_HH

#include <functional>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/geometry/boundingboxtree.hh>
#include <dumux/common/geometry/intersectspointgeometry.hh>
#include <dumux/common/geometry/intersectingentities.hh>

#include <dumux/discretization/methods.hh>

namespace Dumux
{

// forward declarations
template<class TypeTag>
class PointSourceHelper;

/*!
 * \ingroup Common
 * \brief A point source base class
 */
template<class TypeTag>
class PointSource
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;

    static const int dimworld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;

public:
    //! Constructor for constant point sources
    PointSource(GlobalPosition pos, PrimaryVariables values)
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
    PointSource& operator= (const PrimaryVariables& values)
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
    PrimaryVariables values() const
    { return values_; }

    //! return the source position
    const GlobalPosition& position() const
    { return pos_; }

    //! an update function called before adding the value
    // to the local residual in the problem in scvPointSources
    // to be overloaded by derived classes
    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const ElementVolumeVariables &elemVolVars,
                const SubControlVolume &scv)
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
    PrimaryVariables values_; //!< value of the point source for each equation
private:
    GlobalPosition pos_; //!< position of the point source
    std::size_t embeddings_; //!< how many SCVs the point source is associated with
};

/*!
 * \ingroup Common
 * \brief A point source class with an identifier to attach data
 */
template<class TypeTag, typename IdType>
class IdPointSource : public PointSource<TypeTag>
{
    using ParentType = PointSource<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);

    static const int dimworld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;

public:
    //! Constructor for constant point sources
    IdPointSource(GlobalPosition pos, PrimaryVariables values, IdType id)
      :  ParentType(pos, values), id_(id) {}

    //! Constructor for sol dependent point sources, when there is no
    // value known at the time of initialization
    IdPointSource(GlobalPosition pos, IdType id)
      : ParentType(pos, PrimaryVariables(0.0)), id_(id) {}

    //! return the sources identifier
    IdType id() const
    { return id_; }

    //! Convenience = operator overload modifying only the values
    IdPointSource& operator= (const PrimaryVariables& values)
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
class SolDependentPointSource : public PointSource<TypeTag>
{
    using ParentType = PointSource<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;

    static const int dimworld = GridView::dimensionworld;
    using GlobalPosition = typename Dune::FieldVector<Scalar, dimworld>;
    // returns the PointSource values as PrimaryVariables
    using ValueFunction = typename std::function<PrimaryVariables(const Problem &problem,
                                                                  const Element &element,
                                                                  const FVElementGeometry &fvGeometry,
                                                                  const ElementVolumeVariables &elemVolVars,
                                                                  const SubControlVolume &scv)>;

public:
    //! Constructor for sol dependent point sources, when there is no
    // value known at the time of initialization
    SolDependentPointSource(GlobalPosition pos,
                            ValueFunction valueFunction)
      : ParentType(pos, PrimaryVariables(0.0)), valueFunction_(valueFunction) {}

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
    SolDependentPointSource& operator= (const PrimaryVariables& values)
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
template<class TypeTag>
class BoundingBoxTreePointSourceHelper
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using PointSource = typename GET_PROP_TYPE(TypeTag, PointSource);

    static constexpr int dim = GridView::dimension;
    static constexpr int dimworld = GridView::dimensionworld;

    static constexpr bool isBox = GET_PROP_VALUE(TypeTag, DiscretizationMethod) == DiscretizationMethod::box;
    static constexpr int dofCodim = isBox ? dim : 0;

public:
    //! calculate a DOF index to point source map from given vector of point sources
    template<class PointSourceMap>
    static void computePointSourceMap(const FVGridGeometry& fvGridGeometry,
                                      std::vector<PointSource>& sources,
                                      PointSourceMap& pointSourceMap)
    {
        const auto& boundingBoxTree = fvGridGeometry.boundingBoxTree();

        for (auto&& source : sources)
        {
            // compute in which elements the point source falls
            const auto entities = intersectingEntities(source.position(), boundingBoxTree);
            // split the source values equally among all concerned entities
            source.setEmbeddings(entities.size()*source.embeddings());
            // loop over all concernes elements
            for (unsigned int eIdx : entities)
            {
                if(isBox)
                {
                    // check in which subcontrolvolume(s) we are
                    const auto element = boundingBoxTree.entitySet().entity(eIdx);
                    auto fvGeometry = localView(fvGridGeometry);
                    fvGeometry.bindElement(element);

                    const auto globalPos = source.position();
                    // loop over all sub control volumes and check if the point source is inside
                    std::vector<unsigned int> scvIndices;
                    for (auto&& scv : scvs(fvGeometry))
                    {
                        if (intersectsPointGeometry(globalPos, scv.geometry()))
                            scvIndices.push_back(scv.indexInElement());
                    }
                    // for all scvs that where tested positiv add the point sources
                    // to the element/scv to point source map
                    for (auto scvIdx : scvIndices)
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
