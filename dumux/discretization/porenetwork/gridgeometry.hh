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
 * \ingroup PoreNetworkDiscretization
 * \brief Base class for the finite volume geometry for porenetwork models
 */
#ifndef DUMUX_DISCRETIZATION_PNM_GRID_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_PNM_GRID_GEOMETRY_HH

#include <string>
#include <unordered_map>
#include <functional>

#include <dune/common/exceptions.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/discretization/method.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/defaultmappertraits.hh>
#include <dumux/common/parameters.hh>

#include <dumux/discretization/basegridgeometry.hh>
#include <dumux/discretization/porenetwork/fvelementgeometry.hh>
#include <dumux/discretization/porenetwork/subcontrolvolume.hh>
#include <dumux/discretization/porenetwork/subcontrolvolumeface.hh>
#include <dumux/porenetworkflow/common/throatproperties.hh>
#include <dumux/porenetworkflow/common/poreproperties.hh>

namespace Dumux {

/*!
 * \ingroup PoreNetworkDiscretization
 * \brief Base class for geometry data extraction from the grid data format
 */
template<class Scalar, class GridView>
class DefaultPNMData
{
    using GridIndex = typename IndexTraits<GridView>::GridIndex;
    using SmallLocalIndex = typename IndexTraits<GridView>::SmallLocalIndex;
    using Label = std::int_least8_t;
    using Vertex = typename GridView::template Codim<GridView::dimension>::Entity;
    using Element = typename GridView::template Codim<0>::Entity;

    static const int dim = GridView::dimension;

public:

    template<class GridData>
    void update(const GridView& gridView, const GridData& gridData)
    {
        coordinationNumber_ = gridData.getCoordinationNumbers();

        const auto numThroats = gridView.size(0);
        throatRadius_.resize(numThroats);
        throatLength_.resize(numThroats);
        throatLabel_.resize(numThroats);
        throatCrossSectionalArea_.resize(numThroats);
        throatShapeFactor_.resize(numThroats);

        useSameGeometryForAllPores_ = true;
        useSameShapeForAllThroats_ = true;

        // first check if the same geometry shall be used for all entities ...
        if (hasParamInGroup(gridData.paramGroup(), "Grid.ThroatCrossSectionShape"))
        {
            const auto throatGeometryInput = getParamFromGroup<std::string>(gridData.paramGroup(), "Grid.ThroatCrossSectionShape");
            const auto throatGeometry = Throat::shapeFromString(throatGeometryInput);
            throatGeometry_.resize(1);
            throatGeometry_[0] = throatGeometry;

            std::cout << "Using '" << throatGeometryInput << "' as cross-sectional shape for all throats." << std::endl;
        }
        else // .. otherwise, get the corresponding parameter index from the grid data and resize the respective vector
        {
            std::cout << "Reading shape factors for throats from grid data." << std::endl;
            useSameShapeForAllThroats_ = false;
            throatGeometry_.resize(numThroats);
        }

        // get the vertex parameters
        const auto numPores = gridView.size(dim);
        poreRadius_.resize(numPores);
        poreLabel_.resize(numPores);
        poreVolume_.resize(numPores);

        // first check if the same geometry shall be used for all entities ...
        if (hasParamInGroup(gridData.paramGroup(), "Grid.PoreGeometry"))
        {
            const auto poreGeometryInput = getParamFromGroup<std::string>(gridData.paramGroup(), "Grid.PoreGeometry");
            poreGeometry_.resize(1);
            poreGeometry_[0] = Pore::shapeFromString(poreGeometryInput);;

            std::cout << "Using '" << poreGeometryInput << "' as geometry for all pores." << std::endl;
        }
        else // .. otherwise, get the corresponding parameter index from the grid data and resize the respective vector
        {
            std::cout << "Reading pore shapes from grid data." << std::endl;
            useSameGeometryForAllPores_ = false;
            poreGeometry_.resize(numPores);
        }


        for (const auto& vertex : vertices(gridView))
        {
            static const auto poreRadiusIdx = gridData.parameterIndex("PoreRadius");
            static const auto poreLabelIdx = gridData.parameterIndex("PoreLabel");
            const auto vIdx = gridView.indexSet().index(vertex);
            const auto& params = gridData.parameters(vertex);
            poreRadius_[vIdx] = params[poreRadiusIdx];
            assert(poreRadius_[vIdx] > 0.0);
            poreLabel_[vIdx] = params[poreLabelIdx];

            if (!useSameGeometryForAllPores())
                poreGeometry_[vIdx] = getPoreGeometry_(gridData, vertex);

            poreVolume_[vIdx] = getPoreVolume_(gridData, vertex, vIdx);
        }

        for (const auto& element : elements(gridView))
        {
            const int eIdx = gridView.indexSet().index(element);
            const auto& params = gridData.parameters(element);
            static const auto throatRadiusIdx = gridData.parameterIndex("ThroatRadius");
            static const auto throatLengthIdx = gridData.parameterIndex("ThroatLength");
            static const auto throatLabelIdx = gridData.parameterIndex("ThroatLabel");
            throatRadius_[eIdx] = params[throatRadiusIdx];
            throatLength_[eIdx] = params[throatLengthIdx];
            throatLabel_[eIdx] = params[throatLabelIdx];

            if (!useSameShapeForAllThroats())
            {
                static const auto throatShapeFactorIdx = gridData.parameterIndex("ThroatShapeFactor");
                static const auto throatAreaIdx = gridData.parameterIndex("ThroatCrossSectionalArea");
                throatShapeFactor_[eIdx] = params[throatShapeFactorIdx];
                throatGeometry_[eIdx] = Throat::shape(throatShapeFactor_[eIdx]);
                throatCrossSectionalArea_[eIdx] = params[throatAreaIdx];
            }
            else
            {
                throatCrossSectionalArea_[eIdx] = getThroatCrossSectionalArea_(gridData, element, eIdx);
                throatShapeFactor_[eIdx] = getThroatShapeFactor_(gridData, element, eIdx);
            }

            assert(throatRadius_[eIdx] > 0.0);
            assert(throatLength_[eIdx] > 0.0);
            assert(throatCrossSectionalArea_[eIdx] > 0.0);

            static const bool addThroatVolumeToPoreVolume = getParamFromGroup<bool>(gridData.paramGroup(), "Grid.AddThroatVolumeToPoreVolume", false);
            if (addThroatVolumeToPoreVolume)
            {
                for (int vIdxLocal = 0; vIdxLocal < 2; ++vIdxLocal)
                {
                    const auto vIdx = gridView.indexSet().subIndex(element, vIdxLocal, dim);
                    poreVolume_[vIdx] += 0.5 * throatCrossSectionalArea_[eIdx] * throatLength_[eIdx];
                }
            }
        }

        maybeResizeContainers_();
    }

    //! Returns the pore label (e.g. used for setting BCs)
    Label poreLabel(const GridIndex dofIdxGlobal) const
    { return poreLabel_[dofIdxGlobal]; }

    //! Returns the vector of pore labels
    const std::vector<Label>& poreLabel() const
    { return poreLabel_; }

    //! Returns the radius of the pore
    Scalar poreRadius(const GridIndex dofIdxGlobal) const
    { return poreRadius_[dofIdxGlobal]; }

    //! Returns the vector of pore radii
    const std::vector<Scalar>& poreRadius() const
    { return poreRadius_; }

    //! Returns the volume of the pore
    Scalar poreVolume(const GridIndex dofIdxGlobal) const
    { return poreVolume_[dofIdxGlobal]; }

    //! Returns the vector of pore volumes
    const std::vector<Scalar>& poreVolume() const
    { return poreVolume_; }

    //! Returns the radius of the throat
    Scalar throatRadius(const GridIndex eIdx) const
    { return throatRadius_[eIdx]; }

    //! Returns the vector of throat radii
    const std::vector<Scalar>& throatRadius() const
    { return throatRadius_; }

    //! Returns the length of the throat
    Scalar throatLength(const GridIndex eIdx) const
    { return throatLength_[eIdx]; }

    //! Returns the vector of throat lengths
    const std::vector<Scalar>& throatLength() const
    { return throatLength_; }

    //! Returns an index indicating if a throat is touching the domain boundary
    Label throatLabel(const GridIndex eIdx) const
    { return throatLabel_[eIdx]; }

    //! Returns the vector of throat labels
    const std::vector<Label>& throatLabel() const
    { return throatLabel_; }

    //! Returns the number of throats connected to a pore (coordination number)
    SmallLocalIndex coordinationNumber(const GridIndex dofIdxGlobal) const
    { return coordinationNumber_[dofIdxGlobal]; }

    //! Returns the vector of coordination numbers
    const std::vector<SmallLocalIndex>& coordinationNumber() const
    { return coordinationNumber_; }

    //! the geometry of the pore
    Pore::Shape poreGeometry(const GridIndex vIdx) const
    { return useSameGeometryForAllPores() ? poreGeometry_[0] : poreGeometry_[vIdx]; }

    //! Returns the vector of pore geometries
    const std::vector<Pore::Shape>& poreGeometry() const
    {
        if (useSameGeometryForAllPores())
            DUNE_THROW(Dune::InvalidStateException, "All pores have the same geometry. Create your own vector if needed.");
        else
            return poreGeometry_;
    }

    //! Returns the throat's cross-sectional shape
    Throat::Shape throatCrossSectionShape(const GridIndex eIdx) const
    { return useSameShapeForAllThroats() ? throatGeometry_[0] : throatGeometry_[eIdx]; }

    //! Returns the vector of cross-sectional shapes
    const std::vector<Throat::Shape>& throatCrossSectionShape() const
    {
        if (useSameShapeForAllThroats())
            DUNE_THROW(Dune::InvalidStateException, "All throats have the same shape. Create your own vector if needed.");
        else
            return throatGeometry_;
    }

    //! Returns the throat's cross-sectional area
    Scalar throatCrossSectionalArea(const GridIndex eIdx) const
    { return throatCrossSectionalArea_[eIdx]; }

    //! Returns the vector of throat cross-sectional areas
    const std::vector<Scalar>& throatCrossSectionalArea() const
    { return throatCrossSectionalArea_; }

    //! Returns the throat's shape factor
    [[deprecated("Use throatShapeFactor instead. Will be removed soon.")]]
    Scalar shapeFactor(const GridIndex eIdx) const
    { return throatShapeFactor(eIdx); }

    //! Returns the throat's shape factor
    Scalar throatShapeFactor(const GridIndex eIdx) const
    { return useSameShapeForAllThroats() ? throatShapeFactor_[0] : throatShapeFactor_[eIdx]; }

    //! Returns the vector of throat shape factors
    const std::vector<Scalar>& throatShapeFactor() const
    {
        if (useSameShapeForAllThroats())
            DUNE_THROW(Dune::InvalidStateException, "All throats have the same shape. Create your own vector if needed.");
        else
            return throatShapeFactor_;
    }

    //! Returns whether all pores feature the same shape
    bool useSameGeometryForAllPores() const
    { return useSameGeometryForAllPores_; }

    //! Returns whether all throats feature the same cross-sectional shape
    bool useSameShapeForAllThroats() const
    { return useSameShapeForAllThroats_; }

private:

    //! determine the pore geometry provided as scalar value by the grid file
    template<class GridData>
    Pore::Shape getPoreGeometry_(const GridData& gridData, const Vertex& vertex) const
    {
        static const auto poreGeometryIdx = gridData.parameterIndex("PoreGeometry");
        using T = std::underlying_type_t<Pore::Shape>;
        const auto poreGeometryValue = static_cast<T>(gridData.parameters(vertex)[poreGeometryIdx]);
        return static_cast<Pore::Shape>(poreGeometryValue);
    }

    //! automatically determine the pore volume if not provided by the grid file
    template<class GridData>
    Scalar getPoreVolume_(const GridData& gridData, const Vertex& vertex, const std::size_t vIdx) const
    {
        static const bool gridHasPoreVolume = gridData.gridHasVertexParameter("PoreVolume");

        if (gridHasPoreVolume)
        {
            static const auto poreVolumeIdx = gridData.parameterIndex("PoreVolume");
            return gridData.parameters(vertex)[poreVolumeIdx];
        }
        else
        {
            if (poreGeometry(vIdx) == Pore::Shape::cylinder)
            {
                static const Scalar fixedHeight = getParamFromGroup<Scalar>(gridData.paramGroup(), "Grid.PoreHeight", -1.0);
                const Scalar h = fixedHeight > 0.0 ? fixedHeight : gridData.getParameter(vertex, "PoreHeight");
                return Pore::volume(Pore::Shape::cylinder, poreRadius(vIdx), h);
            }
            else
                return Pore::volume(poreGeometry(vIdx), poreRadius(vIdx));
        }
    }

    //! automatically determine throat cross-sectional area if not provided by the grid file
    template<class GridData>
    Scalar getThroatCrossSectionalArea_(const GridData& gridData, const Element& element, const std::size_t eIdx) const
    {
        static const bool gridHasThroatCrossSectionalArea = gridData.gridHasVertexParameter("ThroatCrossSectionalArea");
        if (gridHasThroatCrossSectionalArea)
        {
            static const auto throatAreaIdx = gridData.parameterIndex("ThroatCrossSectionalArea");
            return gridData.parameters(element)[throatAreaIdx];
        }
        else
        {
            if (const auto shape = throatCrossSectionShape(eIdx); shape == Throat::Shape::rectangle)
            {
                static const auto throatHeight = getParamFromGroup<Scalar>(gridData.paramGroup(), "Grid.ThroatHeight");
                return Throat::totalCrossSectionalAreaForRectangle(throatRadius_[eIdx], throatHeight);
            }
            else
                return Throat::totalCrossSectionalArea(shape, throatRadius_[eIdx]);
        }
    }

    //! automatically determine throat shape factor if not provided by the grid file
    template<class GridData>
    Scalar getThroatShapeFactor_(const GridData& gridData, const Element& element, const std::size_t eIdx) const
    {
        static const bool gridHasThroatShapeFactor = gridData.gridHasVertexParameter("ThroatShapeFactor");
        if (gridHasThroatShapeFactor)
        {
            static const auto throatShapeFactorIdx = gridData.parameterIndex("ThroatShapeFactor");
            return gridData.parameters(element)[throatShapeFactorIdx];
        }
        else
        {
            if (const auto shape = throatCrossSectionShape(eIdx); shape == Throat::Shape::rectangle)
            {
                static const auto throatHeight = getParamFromGroup<Scalar>(gridData.paramGroup(), "Grid.ThroatHeight");
                return Throat::shapeFactorRectangle(throatRadius_[eIdx], throatHeight);
            }
            else if (shape == Throat::Shape::polygon || shape == Throat::Shape::scaleneTriangle)
            {
                static const auto shapeFactor = getParamFromGroup<Scalar>(gridData.paramGroup(), "Grid.ThroatShapeFactor");
                return shapeFactor;
            }
            else
                return Throat::shapeFactor<Scalar>(shape, throatRadius_[eIdx]);
        }
    }

    void maybeResizeContainers_()
    {
        // check if all throat might have the same shape in order to save some memory
        if (!useSameShapeForAllThroats() &&
            std::adjacent_find(throatGeometry_.begin(), throatGeometry_.end(), std::not_equal_to<Throat::Shape>() ) == throatGeometry_.end())
        {
            std::cout << "All throats feature the same shape, resizing containers" << std::endl;
            useSameShapeForAllThroats_ = true;
            const Scalar shapeFactor = throatShapeFactor_[0];
            const auto throatGeometry = throatGeometry_[0];
            throatShapeFactor_.resize(1);
            throatGeometry_.resize(1);
            throatShapeFactor_[0] = shapeFactor;
            throatGeometry_[0] = throatGeometry;
        }

        // check if all throat might have the same shape in order to save some memory
        if (!useSameGeometryForAllPores() &&
            std::adjacent_find(poreGeometry_.begin(), poreGeometry_.end(), std::not_equal_to<Pore::Shape>() ) == poreGeometry_.end())
        {
            std::cout << "All pores feature the same shape, resizing containers" << std::endl;
            useSameGeometryForAllPores_ = true;
            const auto poreGeometry = poreGeometry_[0];
            poreGeometry_.resize(1);
            poreGeometry_[0] = poreGeometry;
        }
    }

    std::vector<Pore::Shape> poreGeometry_;
    std::vector<Throat::Shape> throatGeometry_;
    std::vector<SmallLocalIndex> coordinationNumber_;
    std::vector<Scalar> poreRadius_;
    std::vector<Scalar> poreVolume_;
    std::vector<Label> poreLabel_; // 0:no, 1:general, 2:coupling1, 3:coupling2, 4:inlet, 5:outlet
    std::vector<Label> throatLabel_; // 0:no, 1:general, 2:coupling1, 3:coupling2, 4:inlet, 5:outlet
    std::vector<Scalar> throatRadius_;
    std::vector<Scalar> throatLength_;
    std::vector<Scalar> throatShapeFactor_;
    std::vector<Scalar> throatCrossSectionalArea_;
    bool useSameGeometryForAllPores_;
    bool useSameShapeForAllThroats_;
};

/*!
 * \ingroup PoreNetworkDiscretization
 * \brief The default traits
 * \tparam the grid view type
 */
template<class GridView, class MapperTraits = DefaultMapperTraits<GridView>>
struct PNMDefaultGridGeometryTraits
: public MapperTraits
{
    using SubControlVolume = PNMSubControlVolume<GridView>;
    using SubControlVolumeFace = PNMSubControlVolumeFace<GridView>;

    template<class GridGeometry, bool enableCache>
    using LocalView = PNMFVElementGeometry<GridGeometry, enableCache>;

    using PNMData = DefaultPNMData<typename SubControlVolume::Traits::Scalar, GridView>;
};

/*!
 * \ingroup PoreNetworkDiscretization
 * \brief Base class for the finite volume geometry for porenetwork models
 * \note This class is specialized for versions with and without caching the fv geometries on the grid view
 */
template<class Scalar,
         class GridView,
         bool enableGridGeometryCache = false,
         class Traits = PNMDefaultGridGeometryTraits<GridView> >
class PNMGridGeometry;

/*!
 * \ingroup PoreNetworkDiscretization
 * \brief Base class for the finite volume geometry for porenetwork models
 * \note For caching enabled we store the fv geometries for the whole grid view which is memory intensive but faster
 */
template<class Scalar, class GV, class Traits>
class PNMGridGeometry<Scalar, GV, true, Traits>
: public BaseGridGeometry<GV, Traits>
, public Traits::PNMData
{
    using ThisType = PNMGridGeometry<Scalar, GV, true, Traits>;
    using ParentType = BaseGridGeometry<GV, Traits>;
    using GridIndexType = typename IndexTraits<GV>::GridIndex;
    using LocalIndexType = typename IndexTraits<GV>::LocalIndex;
    using PNMData = typename Traits::PNMData;

    using Element = typename GV::template Codim<0>::Entity;
    using CoordScalar = typename GV::ctype;
    static const int dim = GV::dimension;
    static const int dimWorld = GV::dimensionworld;

public:
    //! export discretization method
    static constexpr DiscretizationMethod discMethod = DiscretizationMethod::box;

    //! export the type of the fv element geometry (the local view type)
    using LocalView = typename Traits::template LocalView<ThisType, true>;
    //! export the type of sub control volume
    using SubControlVolume = typename Traits::SubControlVolume;
    //! export the type of sub control volume
    using SubControlVolumeFace = typename Traits::SubControlVolumeFace;
    //! export dof mapper type
    using DofMapper = typename Traits::VertexMapper;
    //! export the finite element cache type
    using FeCache = Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1>;
    //! export the grid view type
    using GridView = GV;

    //! Constructor
    PNMGridGeometry(const GridView gridView)
    : ParentType(gridView)
    {
        static_assert(GridView::dimension == 1, "Porenetwork model only allow GridView::dimension == 1!");
    }

    //! the vertex mapper is the dofMapper
    //! this is convenience to have better chance to have the same main files for box/tpfa/mpfa...
    const DofMapper& dofMapper() const
    { return this->vertexMapper(); }

    //! The total number of sub control volumes
    std::size_t numScv() const
    {  return numScv_; }

    //! The total number of sun control volume faces
    std::size_t numScvf() const
    { return numScvf_; }

    //! The total number of boundary sub control volume faces
    //! For compatibility reasons with cc methods
    std::size_t numBoundaryScvf() const
    { return 0; }

    //! The total number of degrees of freedom
    std::size_t numDofs() const
    { return this->vertexMapper().size(); }

    //! update all fvElementGeometries (do this again after grid adaption)
    template<class GridData>
    void update(const GridData& gridData)
    {
        ParentType::update();
        PNMData::update(this->gridView(), gridData);

        scvs_.clear();
        scvfs_.clear();

        auto numElements = this->gridView().size(0);
        scvs_.resize(numElements);
        scvfs_.resize(numElements);
        hasBoundaryScvf_.resize(numElements, false);

        boundaryDofIndices_.assign(numDofs(), false);

        numScvf_ = numElements;
        numScv_ = 2*numElements;

        // Build the SCV and SCV faces
        for (const auto& element : elements(this->gridView()))
        {
            // get the element geometry
            auto eIdx = this->elementMapper().index(element);
            auto elementGeometry = element.geometry();

            // construct the sub control volumes
            for (LocalIndexType scvLocalIdx = 0; scvLocalIdx < elementGeometry.corners(); ++scvLocalIdx)
            {
                const auto dofIdxGlobal = this->vertexMapper().subIndex(element, scvLocalIdx, dim);

                // get the fractional volume asssociated with the scv
                const auto volume = this->poreVolume(dofIdxGlobal) / this->coordinationNumber(dofIdxGlobal);

                scvs_[eIdx][scvLocalIdx] = SubControlVolume(dofIdxGlobal,
                                                            scvLocalIdx,
                                                            eIdx,
                                                            elementGeometry.corner(scvLocalIdx),
                                                            volume);

                if (this->poreLabel(dofIdxGlobal) > 0)
                {
                    if (boundaryDofIndices_[dofIdxGlobal])
                        continue;

                    boundaryDofIndices_[dofIdxGlobal] = true;
                    hasBoundaryScvf_[eIdx] = true;
                }
            }

            // construct the inner sub control volume face
            auto unitOuterNormal = elementGeometry.corner(1)-elementGeometry.corner(0);
            unitOuterNormal /= unitOuterNormal.two_norm();
            LocalIndexType scvfLocalIdx = 0;
            scvfs_[eIdx][0] = SubControlVolumeFace(elementGeometry.center(),
                                                   std::move(unitOuterNormal),
                                                   scvfLocalIdx++,
                                                   std::vector<LocalIndexType>({0, 1}));
        }
    }

    //! The finite element cache for creating local FE bases
    const FeCache& feCache() const
    { return feCache_; }

    //! Get the local scvs for an element
    const std::array<SubControlVolume, 2>& scvs(GridIndexType eIdx) const
    { return scvs_[eIdx]; }

    //! Get the local scvfs for an element
    const std::array<SubControlVolumeFace, 1>& scvfs(GridIndexType eIdx) const
    { return scvfs_[eIdx]; }

    //! If a vertex / d.o.f. is on the boundary
    bool dofOnBoundary(GridIndexType dofIdx) const
    { return boundaryDofIndices_[dofIdx]; }

    //! If a vertex / d.o.f. is on a periodic boundary (not implemented)
    bool dofOnPeriodicBoundary(GridIndexType dofIdx) const
    { return false; }

    //! The index of the vertex / d.o.f. on the other side of the periodic boundary
    GridIndexType periodicallyMappedDof(GridIndexType dofIdx) const
    { DUNE_THROW(Dune::NotImplemented, "Periodic boundaries"); }

    //! Returns the map between dofs across periodic boundaries
    std::unordered_map<GridIndexType, GridIndexType> periodicVertexMap() const
    { return std::unordered_map<GridIndexType, GridIndexType>{}; }

    //! Returns whether one of the geometry's scvfs lies on a boundary
    bool hasBoundaryScvf(GridIndexType eIdx) const
    { return hasBoundaryScvf_[eIdx]; }

private:
    const FeCache feCache_;

    std::vector<std::array<SubControlVolume, 2>> scvs_;
    std::vector<std::array<SubControlVolumeFace, 1>> scvfs_;

    std::size_t numScv_;
    std::size_t numScvf_;

    // vertices on the boudary
    std::vector<bool> boundaryDofIndices_;
    std::vector<bool> hasBoundaryScvf_;
};

/*!
 * \ingroup PoreNetworkDiscretization
 * \brief Base class for the finite volume geometry for porenetwork models
 * \note For caching disabled we store only some essential index maps to build up local systems on-demand in
 *       the corresponding FVElementGeometry
 */
template<class Scalar, class GV, class Traits>
class PNMGridGeometry<Scalar, GV, false, Traits>
: public BaseGridGeometry<GV, Traits>
, public Traits::PNMData
{
    using ThisType = PNMGridGeometry<Scalar, GV, false, Traits>;
    using ParentType = BaseGridGeometry<GV, Traits>;
    using GridIndexType = typename IndexTraits<GV>::GridIndex;
    using LocalIndexType = typename IndexTraits<GV>::LocalIndex;
    using PNMData = typename Traits::PNMData;

    static const int dim = GV::dimension;
    static const int dimWorld = GV::dimensionworld;

    using Element = typename GV::template Codim<0>::Entity;
    using CoordScalar = typename GV::ctype;

public:
    //! export discretization method
    static constexpr DiscretizationMethod discMethod = DiscretizationMethod::box;

    //! export the type of the fv element geometry (the local view type)
    using LocalView = typename Traits::template LocalView<ThisType, false>;
    //! export the type of sub control volume
    using SubControlVolume = typename Traits::SubControlVolume;
    //! export the type of sub control volume
    using SubControlVolumeFace = typename Traits::SubControlVolumeFace;
    //! export dof mapper type
    using DofMapper = typename Traits::VertexMapper;
    //! export the finite element cache type
    using FeCache = Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1>;
    //! export the grid view type
    using GridView = GV;

    //! Constructor
    PNMGridGeometry(const GridView gridView)
    : ParentType(gridView)
    {
        static_assert(GridView::dimension == 1, "Porenetwork model only allow GridView::dimension == 1!");
    }

    //! the vertex mapper is the dofMapper
    //! this is convenience to have better chance to have the same main files for box/tpfa/mpfa...
    const DofMapper& dofMapper() const
    { return this->vertexMapper(); }

    //! The total number of sub control volumes
    std::size_t numScv() const
    {  return numScv_; }

    //! The total number of sun control volume faces
    std::size_t numScvf() const
    { return numScvf_; }

    //! The total number of boundary sub control volume faces
    //! For compatibility reasons with cc methods
    std::size_t numBoundaryScvf() const
    { return 0; }

    //! The total number of degrees of freedom
    std::size_t numDofs() const
    { return this->vertexMapper().size(); }

    //! update all fvElementGeometries (do this again after grid adaption)
    template<class GridData>
    void update(const GridData& gridData)
    {
        ParentType::update();
        PNMData::update(this->gridView(), gridData);

        boundaryDofIndices_.assign(numDofs(), false);

        // save global data on the grid's scvs and scvfs
        numScvf_ = this->gridView().size(0);
        numScv_ = 2*numScvf_;

        for (const auto& element : elements(this->gridView()))
        {
            // treat boundaries
            for (LocalIndexType vIdxLocal = 0; vIdxLocal < 2; ++vIdxLocal)
            {
                const auto vIdxGlobal = this->vertexMapper().subIndex(element, vIdxLocal, dim);
                if (this->poreLabel(vIdxGlobal) > 0)
                {
                    if (boundaryDofIndices_[vIdxGlobal])
                        continue;

                    boundaryDofIndices_[vIdxGlobal] = true;
                }
            }
        }
    }

    //! The finite element cache for creating local FE bases
    const FeCache& feCache() const
    { return feCache_; }

    //! If a vertex / d.o.f. is on the boundary
    bool dofOnBoundary(GridIndexType dofIdx) const
    { return boundaryDofIndices_[dofIdx]; }

    //! If a vertex / d.o.f. is on a periodic boundary (not implemented)
    bool dofOnPeriodicBoundary(GridIndexType dofIdx) const
    { return false; }

    //! The index of the vertex / d.o.f. on the other side of the periodic boundary
    GridIndexType periodicallyMappedDof(GridIndexType dofIdx) const
    { DUNE_THROW(Dune::NotImplemented, "Periodic boundaries"); }

    //! Returns the map between dofs across periodic boundaries
    std::unordered_map<GridIndexType, GridIndexType> periodicVertexMap() const
    { return std::unordered_map<GridIndexType, GridIndexType>{}; }

private:

    const FeCache feCache_;

    // Information on the global number of geometries
    std::size_t numScv_;
    std::size_t numScvf_;

    // vertices on the boudary
    std::vector<bool> boundaryDofIndices_;
};

} // end namespace Dumux

#endif
