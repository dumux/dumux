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
 * \ingroup InputOutput
 * \brief Class for grid data attached to dgf or gmsh grid files
 */
#ifndef DUMUX_IO_PORENETWORKGRID_DATA_HH
#define DUMUX_IO_PORENETWORKGRID_DATA_HH

#include <vector>
#include <memory>
#include <type_traits>
#include <random>
#include <unordered_map>

#include <dune/common/exceptions.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/utility/persistentcontainer.hh>
#include <dune/geometry/axisalignedcubegeometry.hh>

// FoamGrid specific includes
#if HAVE_DUNE_FOAMGRID
#include <dune/foamgrid/foamgrid.hh>
#include <dune/foamgrid/dgffoam.hh>
#endif

#include <dumux/common/indextraits.hh>
//#include <dumux/common/geometry/intersectspointgeometry.hh>
//#include <dumux/porenetwork/common/geometry.hh>
#include <dumux/porenetwork/common/throatproperties.hh>


namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief Class for grid data attached to dgf or gmsh grid files
 */
template <class Grid>
class PoreNetworkGridData
{
    static constexpr int dim = Grid::dimension;
    static constexpr int dimWorld = Grid::dimensionworld;
    using Intersection = typename Grid::LeafIntersection;
    using Element = typename Grid::template Codim<0>::Entity;
    using Vertex = typename Grid::template Codim<dim>::Entity;
    using GridView = typename Grid::LeafGridView;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using SmallLocalIndex = typename IndexTraits<GridView>::SmallLocalIndex;

    using Scalar = double;
    using BoundaryList = std::array<int, 2*dimWorld>;

    using PersistentParameterContainer = Dune::PersistentContainer<Grid, std::vector<typename Grid::ctype>>;

public:

    //! constructor for dgf grid data
    PoreNetworkGridData(Dune::GridPtr<Grid> grid, const std::string& paramGroup)
    : dgfGrid_(grid)
    , isDgfData_(true)
    , paramGroup_(paramGroup)
    , numSubregions_(0)
    {
        setParameterIndices_();
    }

    //! constructor for non-dgf grid data
    PoreNetworkGridData(std::shared_ptr<Grid> grid, const std::string& paramGroup)
    : factoryGrid_(grid)
    , isDgfData_(false)
    , paramGroup_(paramGroup)
    {
        numSubregions_ = getParamFromGroup<std::size_t>(paramGroup_, "Grid.NumSubregions", 0);
        setParameterIndices_();
        priorityList_ = getPriorityList_();
        boundaryFaceIndex_ = getBoundaryFacemarkerInput_();
        computeBoundingBox_();
    }

    /*!
     * \brief Call the parameters function of the DGF grid pointer if available for vertex data
     * \note You can only pass vertices that exist on level 0!
     */
    const std::vector<double>& parameters(const Vertex& vertex) const
    {
        if (isDgfData_ && !useCopiedDgfData_)
            return dgfGrid_.parameters(vertex);
        else
        {
            assert(!(*vertexParameters_)[vertex].empty() && "No parameters available. Something might be wrong with your grid file!");
            return (*vertexParameters_)[vertex];
        }
    }

    /*!
     * \brief Call the parameters function of the DGF grid pointer if available for element data
     */
    const std::vector<double>& parameters(const Element& element) const
    {
        if (isDgfData_ && !useCopiedDgfData_)
        {
            if (element.hasFather())
            {
                auto level0Element = element;
                while(level0Element.hasFather())
                    level0Element = level0Element.father();

                return dgfGrid_.parameters(level0Element);
            }
            else
            {
                return dgfGrid_.parameters(element);
            }
        }
        else
        {
            assert(!(*elementParameters_)[element].empty() && "No parameters available. Something might be wrong with your grid file!");
            return (*elementParameters_)[element];
        }
    }

    /*!
     * \brief Call the parameters function of the DGF grid pointer if available
     */
    template <class GridImp, class IntersectionImp>
    const Dune::DGFBoundaryParameter::type& parameters(const Dune::Intersection<GridImp, IntersectionImp>& intersection) const
    {
        if (isDgfData_)
            return dgfGrid_.parameters(intersection);
        else
            DUNE_THROW(Dune::InvalidStateException, "The parameters method is only available if the grid was constructed with a DGF file.");
    }

    /*!
     * \brief Returns the boundary face marker index at given position
     *
     * \param pos The current position
     */
    int boundaryFaceMarkerAtPos(const GlobalPosition& pos) const
    {
        static const auto boundaryFaceMarker = getBoundaryFacemarkerInput_();
        const auto priorityList = getPriorityList_();
        constexpr auto eps = 1e-8; //TODO
        // set the priority which decides the order the vertices on the boundary are indexed
        // by default, vertices on min/max faces in x direcetion have the highest priority, followed by y and z
        for (auto i : priorityList)
        {
            const int idx = [i] ()
            {
                if(i < 2)
                    return 0;
                else if(i < 4)
                    return 1;
                else
                    return 2;
            } ();

            auto isOdd = [] (int n) { return (n % 2 != 0); };
            if(isOdd(i))
            {
                if(pos[idx] > bBoxMax_[idx] - eps)
                   return boundaryFaceMarker[i];
            }
            else
            {
                if(pos[idx] < bBoxMin_[idx] + eps)
                    return boundaryFaceMarker[i];
            }
        }
        return -1;
    }

    /*!
     * \brief Computes and returns the label of a given throat
     *
     * \param element The element (throat)
     */
    auto throatLabel(const Element& element) const
    {
        const auto boundaryIdx = parameterIndex("PoreLabel");
        const auto vertex1 = element.template subEntity<dim>(0);
        const auto vertex2 = element.template subEntity<dim>(1);
        const int marker1 = parameters(vertex1)[boundaryIdx];
        const int marker2 = parameters(vertex2)[boundaryIdx];

        if(marker1 == marker2) // both vertices are inside the domain or on the same boundary face
            return marker1;
        if(marker1 == -1) // vertex1 is inside the domain, vertex2 is on a boundary face
            return marker2;
        if(marker2 == -1) // vertex2 is inside the domain, vertex1 is on a boundary face
            return marker1;

        // vertex1 and vertex2 are on different boundary faces
        if(isDgfData_)
        {
            // when assigning the throat label based on the pore label, we need to specify, which pore label is favored
            const auto priorityList = getParamFromGroup<std::vector<int>>(paramGroup_, "Grid.ThroatLabelPriorityList", std::vector<int>{4,3,2} );
            for(const auto i : priorityList)
            {
                if(marker1 == i)
                    return marker1;
                if(marker2 == i)
                    return marker2;
            }
        }
        else
        {
            // use the priority list to find out which pore label is favored
            for(const auto i : priorityList_)
            {
                if(marker1 == boundaryFaceIndex_[i])
                    return marker1;
                if(marker2 == boundaryFaceIndex_[i])
                    return marker2;
            }
        }

        DUNE_THROW(Dune::InvalidStateException, "Something went wrong with the throat labels");
    }

    /*!
     * \brief Returns the coordination numbers for all pore bodies.
     */
    std::vector<SmallLocalIndex> getCoordinationNumbers() const
    {
        std::vector<SmallLocalIndex>  coordinationNumbers(gridView_().size(dim), 0);

        for (const auto &element : elements(gridView_()))
        {
            for (SmallLocalIndex vIdxLocal = 0; vIdxLocal < 2; ++vIdxLocal)
            {
                const auto vIdxGlobal = gridView_().indexSet().subIndex(element, vIdxLocal, dim);
                coordinationNumbers[vIdxGlobal] +=1;
            }
        }

        if (std::any_of(coordinationNumbers.begin(), coordinationNumbers.end(), [](auto i){ return i == 0; }))
            DUNE_THROW(Dune::InvalidStateException, "One of the pores is not connected to another pore. SanitizeGrid will not help in this case. Check your grid file");

        return coordinationNumbers;
    }


    /*!
     * \brief Assign parameters for generically created grids
     */
    void assignParameters()
    {
        if( isDgfData_)
            DUNE_THROW(Dune::InvalidStateException, "Assigning parameter not possible for dgf grids");

        computeBoundingBox_();

        const auto numVertexParams = numSubregions_ > 0 ? 3 : 2; // add PoreRegionId as additional parameter if there are subregions
        const auto numElementParams = numSubregions_ > 0 ? 4 : 3; // add ThroatRegionId as additional parameter if there are subregions
        vertexParameters_ = makeParamContainer_(*factoryGrid_, numVertexParams, 1/*codim*/);
        elementParameters_ = makeParamContainer_(*factoryGrid_, numElementParams, 0/*codim*/);

        using cytpe = typename GlobalPosition::value_type;
        using InternalBoundingBox = Dune::AxisAlignedCubeGeometry<cytpe, dimWorld, dimWorld>;
        std::vector<InternalBoundingBox> internalBoundingBoxes;

        // divide the network into subregions, if specified
        if (numSubregions_ > 0)
        {
            // get bounding boxes of subregions
            for (int i = 0; i < numSubregions_; ++i)
            {
                auto lowerLeft = getParamFromGroup<GlobalPosition>(paramGroup_, "Grid.LowerLeftSubregion" + std::to_string(i));
                auto upperRight = getParamFromGroup<GlobalPosition>(paramGroup_, "Grid.UpperRightSubregion" + std::to_string(i));
                internalBoundingBoxes.emplace_back(std::move(lowerLeft), std::move(upperRight));
            }

            // assign region ids to elements if they are within a subregion
            for (const auto& element : elements(gridView_()))
            {
                // default value for elements not belonging to a subregion
                setParameter_(element, "ThroatRegionId", -1);
                for (int  i = 0; i < numSubregions_; ++i)
                {
                    const auto& subregion = internalBoundingBoxes[i];
                    if (intersectsPointGeometry(element.geometry().center(), subregion))
                        setParameter_(element, "ThroatRegionId", i);
                }
            }
        }

        // get the maximum possible pore body radii such that pore bodies do not intersect
        // (requires ThroatRegionId, if specified)
        const std::vector<Scalar> maxPoreRadius = getMaxPoreRadii_();
        std::vector<bool> poreRadiusLimited(gridView_().size(dim), false);

        // get a helper function for getting the pore radius of a pore body not belonging to a subregion
        auto defaultPoreRadius = poreRadiusHelper_(-1);

        // get helper functions for pore body radii on subregions
        std::vector<decltype(defaultPoreRadius)> subregionPoreRadius;
        for (int i = 0; i < numSubregions_; ++i)
            subregionPoreRadius.emplace_back(poreRadiusHelper_(i));

        const bool useDeprecatedBehavior = getParamFromGroup<bool>(paramGroup_, "Grid.UseDeprecatedRandomParameterBehavior", false);
        if (useDeprecatedBehavior)
        {
            std::cout << "\n*** WARNING: You are using the deprecated assignment of random pore radii which yields different results compared to the new implementation. Will be removed soon! **** \n" << std::endl;
        }

        // treat the pore body parameters (label, radius and maybe regionId)
        for (const auto& vertex : vertices(gridView_()))
        {
            const auto& pos = vertex.geometry().center();
            const auto vIdxGlobal = gridView_().indexSet().index(vertex);
            setParameter_(vertex, "PoreLabel", boundaryFaceMarkerAtPos(pos));

            // sets the minimum of the given value and the maximum possible pore body radius
            // and keeps track of capped radii
            auto setRadiusAndLogIfCapped = [&](const Scalar value)
            {
                if (value > maxPoreRadius[vIdxGlobal])
                {
                    poreRadiusLimited[vIdxGlobal] = true;
                    setParameter_(vertex, "PoreRadius", maxPoreRadius[vIdxGlobal]);
                }
                else
                    setParameter_(vertex, "PoreRadius", value);
            };

            if (!useDeprecatedBehavior) // TODO remove this if condition
            {
                if (numSubregions_ == 0) // assign radius if no subregions are specified
                    setRadiusAndLogIfCapped(defaultPoreRadius(vertex));
                else // assign region ids and radii to vertices if they are within a subregion
                {
                    // default value for elements not belonging to a subregion
                    setParameter_(vertex, "PoreRegionId", -1);
                    for (int  i = 0; i < numSubregions_; ++i)
                    {
                        const auto& subregion = internalBoundingBoxes[i];
                        if (intersectsPointGeometry(vertex.geometry().center(), subregion))
                            setParameter_(vertex, "PoreRegionId", i);
                    }

                    // set the pore radius
                    if (const auto id = getParameter(vertex, "PoreRegionId"); id >= 0)
                        setRadiusAndLogIfCapped(subregionPoreRadius[id](vertex));
                    else
                        setRadiusAndLogIfCapped(defaultPoreRadius(vertex));
                }
            }
        }

        // treat throat parameters
        auto defaultThroatRadius = throatRadiusHelper_(-1);
        auto defaultThroatLength = throatLengthHelper_(-1);

        // get helper functions for throat radii and lengths on subregions
        std::vector<decltype(defaultThroatRadius)> subregionThroatRadius;
        std::vector<decltype(defaultThroatLength)> subregionThroatLength;
        for (int i = 0; i < numSubregions_; ++i)
        {
            subregionThroatRadius.emplace_back(throatRadiusHelper_(i));
            subregionThroatLength.emplace_back(throatLengthHelper_(i));
        }

        // set throat parameters
        for (const auto& element : elements(gridView_()))
        {
            if (numSubregions_ == 0) // assign values if no subregions are specified
            {
                if (useDeprecatedBehavior) // TODO remove the if condition and the complete block within
                {
                    // This is the old way of iterating over vertices. We now iterate directly over them (whithout going over the elements first)
                    // which changes the order of assiging values to them and therefore the behavior of the random generator
                    typedef typename GridView::template Codim<dim>::Entity PNMVertex;
                    const std::array<PNMVertex, 2> vertices = {element.template subEntity<dim>(0), element.template subEntity<dim>(1)};
                    for (int i = 0; i < 2; ++i)
                    {
                        auto vertex = vertices[i];
                        const Scalar poreRadius = defaultPoreRadius(vertex);
                        const auto vIdxGlobal = gridView_().indexSet().subIndex(element, i, dim);

                        if (poreRadius > maxPoreRadius[vIdxGlobal])
                        {
                            poreRadiusLimited[vIdxGlobal] = true;
                            setParameter_(vertex, "PoreRadius", maxPoreRadius[vIdxGlobal]);
                        }
                        else
                            setParameter_(vertex, "PoreRadius", poreRadius);
                    }
                }

                setParameter_(element, "ThroatRadius", defaultThroatRadius(element));
                setParameter_(element, "ThroatLength", defaultThroatLength(element));
            }
            else // assign values to throats if they are within a subregion
            {
                if (const auto id = getParameter(element, "ThroatRegionId"); id >= 0)
                {
                    setParameter_(element, "ThroatRadius", subregionThroatRadius[id](element));
                    setParameter_(element, "ThroatLength", subregionThroatLength[id](element));
                }
                else
                {
                    setParameter_(element, "ThroatRadius", defaultThroatRadius(element));
                    setParameter_(element, "ThroatLength", defaultThroatLength(element));
                }
            }

            // set the throat label
            setParameter_(element, "ThroatLabel", throatLabel(element));
        }

        const auto numPoreRadiusLimited = std::count(poreRadiusLimited.begin(), poreRadiusLimited.end(), true);
        if (numPoreRadiusLimited > 0)
            std::cout << "*******\nWarning!  " << numPoreRadiusLimited << " out of " << poreRadiusLimited.size()
            << " pore body radii have been capped automatically in order to prevent intersecting pores\n*******" << std::endl;
    }

    void resizeParameterContainers()
    {
        // resize the parameters
        vertexParameters_->resize();
        elementParameters_->resize();
        vertexParameters_->shrinkToFit();
        elementParameters_->shrinkToFit();
    }

    void copyDgfData()
    {
        if (!isDgfData_)
            DUNE_THROW(Dune::InvalidStateException, "copying dgf data only works when a dgf grid is actually used");

        useCopiedDgfData_ = true;
        const auto someVertex = *(vertices(gridView_()).begin());
        const auto someElement = *(elements(gridView_()).begin());
        const auto numVertexParams = dgfGrid_.parameters(someVertex).size();
        const auto numElementParams = dgfGrid_.parameters(someElement).size();
        vertexParameters_ = makeParamContainer_(*dgfGrid_, numVertexParams, 1);
        elementParameters_ = makeParamContainer_(*dgfGrid_, numElementParams, 0);

        for (const auto& element : elements(gridView_()))
        {
            for (int i = 0; i < numElementParams; ++i)
                (*elementParameters_)[element][i] = dgfGrid_.parameters(element)[i];
        }

        for (const auto& vertex : vertices(gridView_()))
        {
            for (int i = 0; i < numVertexParams; ++i)
                (*vertexParameters_)[vertex][i] = dgfGrid_.parameters(vertex)[i];
        }
    }

    /*!
     * \brief Return the index for a given parameter name
     */
    int parameterIndex(const std::string& paramName) const
    {
        // make sure the string is present in the map, throw a Dumux exception otherwise (and not a std one)
        // the [] operator can't be used here due to const correctness
        if (parameterIndex_.count(paramName))
            return parameterIndex_.at(paramName);
        else
        {
            std::string msg;
            if (paramName.find("Throat") != std::string::npos)
                msg = "Make sure to include it in the vector of parameter names ElementParameters = " + paramName + " ... ...";
            else if (paramName.find("Pore") != std::string::npos)
                msg = "Make sure to include it in the vector of parameter names VertexParameters = " + paramName + " ... ...";

            DUNE_THROW(Dumux::ParameterException, paramName << " not set in the input file. \n" << msg);
        }
    }

    /*!
     * \brief Return the parameter group
     */
    const std::string& paramGroup() const
    { return paramGroup_; }

    /*!
     * \brief Return if a given element parameter is provided by the grid
     */
    bool gridHasElementParameter(const std::string& param) const
    {
        return std::any_of(elementParameterNames_.begin(), elementParameterNames_.end(), [&param]( const auto& i ){ return (i == param); });
    }

    /*!
     * \brief Return if a given vertex parameter is provided by the grid
     */
    bool gridHasVertexParameter(const std::string& param) const
    {
        return std::any_of(vertexParameterNames_.begin(), vertexParameterNames_.end(), [&param]( const auto& i ){ return (i == param); });
    }

    /*!
     * \brief Returns the value of an element parameter
     */
    Scalar getParameter(const Element& element, const std::string& param) const
    { return (*elementParameters_)[element][parameterIndex(param)]; }

    /*!
     * \brief Returns the value of an vertex parameter
     */
    Scalar getParameter(const Vertex& vertex, const std::string& param) const
    { return (*vertexParameters_)[vertex][parameterIndex(param)]; }

private:


    void setParameter_(const Element& element, const std::string& param, const Scalar value)
    { (*elementParameters_)[element][parameterIndex(param)] = value; }

    void setParameter_(const Vertex& vertex, const std::string& param, const Scalar value)
    { (*vertexParameters_)[vertex][parameterIndex(param)] = value; }

    void setParameterIndices_()
    {
        using StringVector = std::vector<std::string>;
        const auto defaultElementParameterNames = numSubregions_ > 0 ? StringVector{"ThroatRadius", "ThroatLength", "ThroatRegionId", "ThroatLabel"} : StringVector{"ThroatRadius", "ThroatLength", "ThroatLabel"};
        const auto defaultVertexParameterNames = numSubregions_ > 0 ? StringVector{"PoreRadius", "PoreRegionId", "PoreLabel"} : StringVector{"PoreRadius", "PoreLabel"};
        elementParameterNames_ = std::move(getParamFromGroup<StringVector>(paramGroup_, "Grid.ElementParameters", defaultElementParameterNames));
        vertexParameterNames_ = std::move(getParamFromGroup<StringVector>(paramGroup_, "Grid.VertexParameters", defaultVertexParameterNames));

        // make sure that the number of specified parameters matches with the dgf file
        if (isDgfData_)
        {
            const auto someElement = *(elements(gridView_()).begin());
            const auto someVertex = *(vertices(gridView_()).begin());

            if (elementParameterNames_.size() != dgfGrid_.nofParameters(someElement))
                DUNE_THROW(Dune::InvalidStateException, "Number of user-specified element parameters (" << elementParameterNames_.size()
                            << ") does not match number of element paramters in dgf file (" << dgfGrid_.nofParameters(someElement) << ")");

            if (vertexParameterNames_.size() != dgfGrid_.nofParameters(someVertex))
                DUNE_THROW(Dune::InvalidStateException, "Number of user-specified vertex parameters (" << vertexParameterNames_.size()
                            << ") does not match number of vertex paramters in dgf file (" << dgfGrid_.nofParameters(someVertex) << ")");
        }

        for (int i = 0; i < vertexParameterNames_.size(); ++i)
        {
            std::cout << vertexParameterNames_[i] << " is vertex parameter " << i << std::endl;
            parameterIndex_[vertexParameterNames_[i]] = i;
        }

        for (int i = 0; i < elementParameterNames_.size(); ++i)
        {
            std::cout << elementParameterNames_[i] << " is element parameter " << i << std::endl;
            parameterIndex_[elementParameterNames_[i]] = i;
        }
    }

    // returns the maximum possible pore body radii such that pore bodies do not intersect
    std::vector<Scalar> getMaxPoreRadii_() const
    {
        const auto numVertices = gridView_().size(dim);
        std::vector<Scalar> maxPoreRadius(numVertices, std::numeric_limits<Scalar>::max());

        if (!getParamFromGroup<bool>(paramGroup_, "Grid.CapPoreRadii", true))
            return maxPoreRadius;

        // check for a user-specified fixed throat length
        const Scalar inputThroatLength = getParamFromGroup<Scalar>(paramGroup_, "Grid.ThroatLength", -1.0);
        std::vector<Scalar> subregionInputThroatLengths;

        if (numSubregions_ > 0)
        {
            for (int i = 0; i < numSubregions_; ++i)
            {
                // adapt the parameter group if there are subregions
                const std::string paramGroup = paramGroup_ + ".SubRegion" + std::to_string(i);
                const Scalar subregionInputThroatLength = getParamFromGroup<Scalar>(paramGroup, "Grid.ThroatLength", -1.0);
                subregionInputThroatLengths.push_back(subregionInputThroatLength);
            }
        }

        for (const auto& element : elements(gridView_()))
        {
            // do not cap the pore radius if a user-specified throat length if given
            if (numSubregions_ > 0)
            {
                // check subregions
                const auto subregionId = getParameter(element, "ThroatRegionId");
                if (subregionId >= 0)
                {
                    // throat lies within a subregion
                    if (subregionInputThroatLengths[subregionId] > 0.0)
                        continue;
                }
                else if (inputThroatLength > 0.0) // throat does not lie within a subregion
                    continue;
            }
            else if (inputThroatLength > 0.0) // no subregions
                continue;

            // No fixed throat lengths given, check for max. pore radius.
            // We define this as 1/2 of the length (minus a user specified value) of the shortest pore throat attached to the pore body.
            const Scalar delta = element.geometry().volume();
            static const Scalar minThroatLength = getParamFromGroup<Scalar>(paramGroup_, "Grid.MinThroatLength", 1e-6);
            const Scalar maxRadius = (delta - minThroatLength)/2.0;
            for (int vIdxLocal = 0; vIdxLocal < 2 ; ++vIdxLocal)
            {
                const int vIdxGlobal = gridView_().indexSet().subIndex(element, vIdxLocal, dim);
                maxPoreRadius[vIdxGlobal] = std::min(maxPoreRadius[vIdxGlobal], maxRadius);
            }
        }

        return maxPoreRadius;
    }

    // returns a lambda taking a vertex and returning a radius
    auto poreRadiusHelper_(int subregionId) const
    {
        // adapt the parameter group if there are subregions
        const std::string paramGroup = subregionId < 0 ? paramGroup_ : paramGroup_ + ".Subregion" + std::to_string(subregionId);

        // prepare random number generation for lognormal parameter distribution
        std::mt19937 generator;

        // lognormal random number distribution for pore radii
        std::lognormal_distribution<> poreRadiusDist;

        const Scalar fixedPoreRadius = getParamFromGroup<Scalar>(paramGroup, "Grid.PoreRadius", -1);
        if (fixedPoreRadius < 0.0)
        {
            const auto type = getParamFromGroup<std::string>(paramGroup, "Grid.ParameterType");
            if (type != "lognormal")
                DUNE_THROW(Dune::InvalidStateException, "Unknown parameter type " << type);

            // allow to specify a seed to get reproducible results
            const auto seed = getParamFromGroup<unsigned int>(paramGroup, "Grid.ParameterRandomNumberSeed", std::random_device{}());
            generator.seed(seed);

            // if we use a distribution, get the mean and standard deviation from input file
            const Scalar meanPoreRadius = getParamFromGroup<Scalar>(paramGroup, "Grid.MeanPoreRadius");
            const Scalar stddevPoreRadius = getParamFromGroup<Scalar>(paramGroup, "Grid.StandardDeviationPoreRadius");
            const Scalar variance = stddevPoreRadius*stddevPoreRadius;

            using std::log;
            using std::sqrt;
            const Scalar mu = log(meanPoreRadius/sqrt(1.0 + variance/(meanPoreRadius*meanPoreRadius)));
            const Scalar sigma = sqrt(log(1.0 + variance/(meanPoreRadius*meanPoreRadius)));

            std::lognormal_distribution<>::param_type params(mu, sigma);
            poreRadiusDist.param(params);
        }

        // check if pores for certain labels should be treated in a special way
        const auto poreLabelsToSetFixedRadius = getParamFromGroup<std::vector<int>>(paramGroup, "Grid.PoreLabelsToSetFixedRadius", std::vector<int>{});
        const auto poreLabelsToApplyFactorForRadius = getParamFromGroup<std::vector<int>>(paramGroup, "Grid.PoreLabelsToApplyFactorForRadius", std::vector<int>{});
        const auto poreRadiusForLabel = getParamFromGroup<std::vector<Scalar>>(paramGroup, "Grid.FixedPoreRadiusForLabel", std::vector<Scalar>{});
        const auto poreRadiusFactorForLabel = getParamFromGroup<std::vector<Scalar>>(paramGroup, "Grid.PoreRadiusFactorForLabel", std::vector<Scalar>{});

        auto maybeModifiedRadius = [&, paramGroup, fixedPoreRadius, poreRadiusDist, generator,
                                    poreLabelsToSetFixedRadius, poreLabelsToApplyFactorForRadius,
                                    poreRadiusForLabel, poreRadiusFactorForLabel](const auto& vertex) mutable
        {
            // default radius to return: either a fixed (user-specified) one or a randomly drawn one
            auto radius = [&]() { return fixedPoreRadius > 0.0 ? fixedPoreRadius : poreRadiusDist(generator); };

            // check if pores for certain labels should be treated in a special way
            if (poreLabelsToSetFixedRadius.empty() && poreLabelsToApplyFactorForRadius.empty())
                return radius(); // nothing special to be done

            // get information on how to manipulate the pore radius
            const auto poreLabel = (*vertexParameters_)[vertex][parameterIndex("PoreLabel")];

            // set a fixed radius for a given label
            if (!poreLabelsToSetFixedRadius.empty() || !poreRadiusForLabel.empty())
            {
                if (poreLabelsToSetFixedRadius.size() != poreRadiusForLabel.size())
                    DUNE_THROW(Dumux::ParameterException, "PoreLabelsToSetFixedRadius must be of same size as FixedPoreRadiusForLabel");

                if (const auto it = std::find(poreLabelsToSetFixedRadius.begin(), poreLabelsToSetFixedRadius.end(), poreLabel); it != poreLabelsToSetFixedRadius.end())
                    return poreRadiusForLabel[std::distance(poreLabelsToSetFixedRadius.begin(), it)];
            }

            // multiply the pore radius by a given value for a given label
            if (!poreLabelsToApplyFactorForRadius.empty() || !poreRadiusFactorForLabel.empty())
            {
                if (poreLabelsToApplyFactorForRadius.size() != poreRadiusFactorForLabel.size())
                    DUNE_THROW(Dumux::ParameterException, "PoreLabelsToApplyFactorForRadius must be of same size as PoreRadiusFactorForLabel");

                if (const auto it = std::find(poreLabelsToApplyFactorForRadius.begin(), poreLabelsToApplyFactorForRadius.end(), poreLabel); it != poreLabelsToApplyFactorForRadius.end())
                    return poreRadiusFactorForLabel[std::distance(poreLabelsToApplyFactorForRadius.begin(), it)] * radius();
            }

            // default
            return radius();
        };

        return maybeModifiedRadius;
    }

    // returns a lambda taking a element and returning a radius
    auto throatRadiusHelper_(int subregionId) const
    {
        // adapt the parameter group if there are subregions
        const std::string paramGroup = subregionId < 0 ? paramGroup_ : paramGroup_ + ".Subregion" + std::to_string(subregionId);

        // shape parameter for calculation of throat radius
        const Scalar throatN = getParamFromGroup<Scalar>(paramGroup, "Grid.ThroatRadiusN", -1); //recommended value if set 0.1?

        // check for a user-specified fixed throat radius
        const Scalar fixedThroatRadius = getParamFromGroup<Scalar>(paramGroup, "Grid.ThroatRadius", -1);

        // prepare random number generation for lognormal parameter distribution
        std::mt19937 generatorThroat;

        // lognormal random number distribution for throat radii
        std::lognormal_distribution<> throatRadiusDist;

        //specify lognormal distribution for throat radii
        if (fixedThroatRadius < 0.0)
        {
            const auto type = getParamFromGroup<std::string>(paramGroup, "Grid.ParameterType");
            if (type != "lognormal")
                DUNE_THROW(Dune::InvalidStateException, "Unknown parameter type " << type);

            // allow to specify a seed to get reproducible results
            const auto seed = getParamFromGroup<unsigned int>(paramGroup, "Grid.ParameterRandomNumberSeed", std::random_device{}());
            generatorThroat.seed(seed);

            // if we use a distribution, get the mean and standard deviation from input file
            const Scalar meanThroatRadius = getParamFromGroup<Scalar>(paramGroup, "Grid.MeanThroatRadius", -1);
            const Scalar stddevThroatRadius = getParamFromGroup<Scalar>(paramGroup, "Grid.StandardDeviationThroatRadius");
            const Scalar throatVariance = stddevThroatRadius*stddevThroatRadius;

            using std::log;
            using std::sqrt;
            const Scalar muThroat = log(meanThroatRadius/sqrt(1.0 + throatVariance/(meanThroatRadius*meanThroatRadius)));
            const Scalar sigmaThroat = sqrt(log(1.0 + throatVariance/(meanThroatRadius*meanThroatRadius)));

            std::lognormal_distribution<>::param_type paramsThroat(muThroat, sigmaThroat);
            throatRadiusDist.param(paramsThroat);
        }

        auto getThroatRadius = [&, paramGroup, fixedThroatRadius, throatN, throatRadiusDist, generatorThroat](const Element& element) mutable
        {
            //check if value for throat size distribution is given
            const Scalar meanThroatRadius = getParamFromGroup<Scalar>(paramGroup, "Grid.MeanThroatRadius", -1);

            // the element parameters (throat radius and length)
            if (fixedThroatRadius > 0.0)
                return fixedThroatRadius;
            else if (throatN > 0.0)
            {
                const Scalar delta = element.geometry().volume();
                typedef typename GridView::template Codim<dim>::Entity PNMVertex;
                const std::array<PNMVertex, 2> vertices = {element.template subEntity<dim>(0), element.template subEntity<dim>(1)};

                const Scalar poreRadius0 = parameters(vertices[0])[parameterIndex("PoreRadius")];
                const Scalar poreRadius1 = parameters(vertices[1])[parameterIndex("PoreRadius")];
                return Throat::averagedRadius(poreRadius0, poreRadius1, delta, throatN);
            }
            else if (meanThroatRadius > 0.0)
            {
                return throatRadiusDist(generatorThroat);
            }
            else
                DUNE_THROW(Dune::InvalidStateException, "You must specify a value for throat radius");
        };

        return getThroatRadius;
    }

    // returns a lambda taking a element and returning a length
    auto throatLengthHelper_(int subregionId) const
    {
        auto getThroatLength = [&, subregionId](const Element& element)
        {
            // adapt the parameter group if there are subregions
            const std::string paramGroup = subregionId < 0 ? paramGroup_ : paramGroup_ + ".Subregion" + std::to_string(subregionId);

            static const Scalar inputThroatLength = getParamFromGroup<Scalar>(paramGroup, "Grid.ThroatLength", -1.0);
            if (inputThroatLength > 0.0)
                return inputThroatLength;

            const Scalar delta = element.geometry().volume();

            // decide whether to substract the pore radii from the throat length or not
            static const bool substractRadiiFromThroatLength = getParamFromGroup<bool>(paramGroup, "Grid.SubstractRadiiFromThroatLength", true);

            if (substractRadiiFromThroatLength)
            {
                typedef typename GridView::template Codim<dim>::Entity PNMVertex;
                const std::array<PNMVertex, 2> vertices = {element.template subEntity<dim>(0), element.template subEntity<dim>(1)};
                const Scalar result = delta - getParameter(vertices[0], "PoreRadius") - getParameter(vertices[1], "PoreRadius");
                if (result <= 0.0)
                    DUNE_THROW(Dune::GridError, "Pore radii are so large they intersect! Something went wrong at element " << gridView_().indexSet().index(element));
                else
                    return result;
            }
            else
                return delta;
        };

        return getThroatLength;
    }

    /*!
     * \brief Initializes and returns a container for vertex (codim dim) or element (codim 0) data
     *
     * \param grid The grid
     * \param numParams The number of paramters
     * \param codim The codimension
     */
    auto makeParamContainer_(const Grid& grid, int numParams, int codim) const
    {
        auto parameters = std::make_unique<PersistentParameterContainer>(grid, codim);
        (*parameters).resize();
        for (auto&& v : (*parameters))
            v.resize(numParams);
        return std::move(parameters);
    }

    /*!
     * \brief Returns a list of boundary face indices from user specified input or default values if no input is given
     */
    BoundaryList getBoundaryFacemarkerInput_() const
    {
        BoundaryList boundaryFaceMarker;
        std::fill(boundaryFaceMarker.begin(), boundaryFaceMarker.end(), 0);
        boundaryFaceMarker[0] = 1;
        boundaryFaceMarker[1] = 1;

        if(hasParamInGroup(paramGroup_, "Grid.BoundaryFaceMarker"))
        {
            try {
                boundaryFaceMarker = getParamFromGroup<BoundaryList>(paramGroup_, "Grid.BoundaryFaceMarker");
            }
            catch (Dune::RangeError& e) {
                DUNE_THROW(Dumux::ParameterException, "You must specifiy all boundaries faces: xmin xmax ymin ymax (zmin zmax). \n" << e.what());
            }
            if(std::none_of(boundaryFaceMarker.begin(), boundaryFaceMarker.end(), []( const int i ){ return i == 1; }))
                DUNE_THROW(Dumux::ParameterException, "At least one face must have index 1");
            if(std::any_of(boundaryFaceMarker.begin(), boundaryFaceMarker.end(), []( const int i ){ return (i < 0 || i > 2*dimWorld); }))
                DUNE_THROW(Dumux::ParameterException, "Face indices must range from 0 to " << 2*dimWorld );
        }
        return boundaryFaceMarker;
    }

    /*!
     * \brief Returns a list of boundary face priorities from user specified input or default values if no input is given
     *
     * This essentially determines the index of a node on an edge or corner corner. For instance, a list of {0,1,2} will give highest priority
     * to the "x"-faces and lowest to the "z-faces".
     */
    BoundaryList getPriorityList_() const
    {
        const auto list = [&]()
        {
            BoundaryList priorityList;
            std::iota(priorityList.begin(), priorityList.end(), 0);

            if(hasParamInGroup(paramGroup_, "Grid.PriorityList"))
            {
                try {
                    // priorities can also be set in the input file
                    priorityList = getParamFromGroup<BoundaryList>(paramGroup_, "Grid.PriorityList");
                }
                // make sure that a priority for each direction is set
                catch(Dune::RangeError& e) {
                    DUNE_THROW(Dumux::ParameterException, "You must specifiy priorities for all directions (" << dimWorld << ") \n" << e.what());
                }
                // make sure each direction is only set once
                auto isUnique = [] (auto v)
                {
                    std::sort(v.begin(), v.end());
                    return (std::unique(v.begin(), v.end()) == v.end());
                };
                if(!isUnique(priorityList))
                    DUNE_THROW(Dumux::ParameterException, "You must specifiy priorities for all directions (duplicate directions)");

                //make sure that the directions are correct (ranging from 0 to dimWorld-1)
                if(std::any_of(priorityList.begin(), priorityList.end(), []( const int i ){ return (i < 0 || i >= 2*dimWorld); }))
                    DUNE_THROW(Dumux::ParameterException, "You must specifiy priorities for correct directions (0-" << 2*(dimWorld-1) << ")");
            }
            return priorityList;
        }();
        return list;
    }


    /*!
     * \brief Return the gridView this grid geometry object lives on
     */
    const GridView gridView_() const
    {
        if(isDgfData_)
            return dgfGrid_->leafGridView();
        else
            return factoryGrid_->leafGridView();
    }

    void computeBoundingBox_()
    {
        // calculate the bounding box of the local partition of the grid view
        for (const auto& vertex : vertices(gridView_()))
        {
            for (int i = 0; i < dimWorld; i++)
            {
                using std::min;
                using std::max;
                bBoxMin_[i] = min(bBoxMin_[i], vertex.geometry().corner(0)[i]);
                bBoxMax_[i] = max(bBoxMax_[i], vertex.geometry().corner(0)[i]);
            }
        }
    }


    // dgf grid data
    Dune::GridPtr<Grid> dgfGrid_;

    std::shared_ptr<Grid> factoryGrid_;

    bool isDgfData_ = false;
    bool useCopiedDgfData_ = false;
    std::string paramGroup_;

    std::vector<std::string> vertexParameterNames_;
    std::vector<std::string> elementParameterNames_;

    std::array<unsigned int, 2> vertexParamterIndices_;
    std::array<unsigned int, 3> elementParamterIndices_;

    std::unique_ptr<PersistentParameterContainer> vertexParameters_;
    std::unique_ptr<PersistentParameterContainer> elementParameters_;

    //! the bounding box of the whole domain
    GlobalPosition bBoxMin_ = GlobalPosition(std::numeric_limits<double>::max());
    GlobalPosition bBoxMax_ = GlobalPosition(std::numeric_limits<double>::min());

    BoundaryList priorityList_;
    BoundaryList boundaryFaceIndex_;

    std::size_t numSubregions_;

    std::unordered_map<std::string, int> parameterIndex_;
};

} // namespace Dumux

#endif
