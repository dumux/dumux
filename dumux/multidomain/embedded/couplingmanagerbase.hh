// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedCoupling
 * \brief Coupling manager for low-dimensional domains embedded in the bulk
 *        domain. Point sources on each integration point are computed by an AABB tree.
 */

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGERBASE_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGERBASE_HH

#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <iostream>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <span>
#include <string>
#include <tuple>
#include <utility>
#include <unordered_map>

#include <dune/common/timer.hh>
#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/io/format.hh>
#include <dumux/parallel/parallel_for.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/griddiscretization.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/numericdifferentiation.hh>
#include <dumux/common/parameters.hh>
#include <dumux/assembly/coloring.hh>
#include <dumux/geometry/distance.hh>
#include <dumux/geometry/intersectingentities.hh>
#include <dumux/discretization/method.hh>
#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/common/typetraits/boundary_.hh>
#include <dumux/common/typetraits/periodic.hh>
#include <dumux/multidomain/couplingmanager.hh>
#include <dumux/multidomain/glue.hh>
#include <dumux/multidomain/distributedglue.hh>
#include <dumux/parallel/distributedsums.hh>
#include <dumux/common/pointsources.hh>
#include <dumux/common/concepts/variables_.hh>
#include <dumux/multidomain/embedded/distributedpointsources.hh>
#include <dumux/multidomain/embedded/elementdofs.hh>
#include <dumux/multidomain/embedded/pointsourcecontext.hh>
#include <dumux/multidomain/embedded/pointsourcedata.hh>
#include <dumux/multidomain/embedded/integrationpointsource.hh>
#include <dumux/multidomain/assemblytraits.hh>

namespace Dumux {

//! the default point source traits
template<class MDTraits>
struct DefaultPointSourceTraits
{
private:
    template<std::size_t i> using SubDomainTypeTag = typename MDTraits::template SubDomain<i>::TypeTag;
    template<std::size_t i> using GridGeometry = GetPropType<SubDomainTypeTag<i>, Properties::GridGeometry>;
    template<std::size_t i> using NumEqVector = Dumux::NumEqVector<GetPropType<SubDomainTypeTag<i>, Properties::PrimaryVariables>>;
public:
    //! export the point source type for domain i
    template<std::size_t i>
    using PointSource = IntegrationPointSource<typename GridGeometry<i>::GlobalCoordinate, NumEqVector<i>>;

    //! export the point source helper type  for domain i
    template<std::size_t i>
    using PointSourceHelper = IntegrationPointSourceHelper;

    //! export the point source data type
    using PointSourceData = Dumux::PointSourceData<MDTraits>;
};

/*!
 * \ingroup EmbeddedCoupling
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *        Point sources on each integration point are computed by an AABB tree.
 */
template<class MDTraits, class Implementation, class PSTraits = DefaultPointSourceTraits<MDTraits>>
class EmbeddedCouplingManagerBase
: public CouplingManager<MDTraits>
{
    using ParentType = CouplingManager<MDTraits>;
    using Scalar = typename MDTraits::Scalar;
    static constexpr auto bulkIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto lowDimIdx = typename MDTraits::template SubDomain<1>::Index();
    using SolutionVector = typename MDTraits::SolutionVector;
    using PointSourceData = typename PSTraits::PointSourceData;

    // the sub domain type tags
    template<std::size_t id> using PointSource = typename PSTraits::template PointSource<id>;
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using PrimaryVariables = GetPropType<SubDomainTypeTag<id>, Properties::PrimaryVariables>;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using ElementMapper = typename GridGeometry<id>::ElementMapper;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using ElementSeed = typename GridView<id>::Grid::template Codim<0>::EntitySeed;
    template<std::size_t id> using GridIndex = typename IndexTraits<GridView<id>>::GridIndex;
    template<std::size_t id> using CouplingStencil = std::vector<GridIndex<id>>;

    static constexpr int bulkDim = GridView<bulkIdx>::dimension;
    static constexpr int lowDimDim = GridView<lowDimIdx>::dimension;
    static constexpr int dimWorld = GridView<bulkIdx>::dimensionworld;

    template<std::size_t id>
    static constexpr bool isCVFE()
    { return DiscretizationMethods::isCVFE<typename GridGeometry<id>::DiscretizationMethod>; }

    template<std::size_t id>
    static constexpr bool isSupportedDiscretizationMethod()
    {
        return GridGeometry<id>::discMethod == DiscretizationMethods::cctpfa
            || GridGeometry<id>::discMethod == DiscretizationMethods::ccmpfa
            || isCVFE<id>();
    }

    static_assert(isSupportedDiscretizationMethod<bulkIdx>() && isSupportedDiscretizationMethod<lowDimIdx>(),
                  "The embedded coupling managers support the cell-centered methods and the control-volume finite element methods");

    using GlobalPosition = typename Element<bulkIdx>::Geometry::GlobalCoordinate;
    using GlueType = MultiDomainGlue<GridView<bulkIdx>, GridView<lowDimIdx>, ElementMapper<bulkIdx>, ElementMapper<lowDimIdx>>;
    using BulkGeometryCopy = Dune::MultiLinearGeometry<Scalar, bulkDim, dimWorld>;
    using PieceGeometry = Dune::AffineGeometry<Scalar, lowDimDim, dimWorld>;

    static constexpr std::size_t invalidIndex = std::numeric_limits<std::size_t>::max();

public:
    //! export traits
    using MultiDomainTraits = MDTraits;
    //! export the point source traits
    using PointSourceTraits = PSTraits;
    //! whether the domain is assembled with the variables of the local dofs and the element discretization interface
    template<std::size_t id>
    static constexpr bool usesElementVariables()
    {
        using GridVariables = GetPropType<SubDomainTypeTag<id>, Properties::GridVariables>;
        return Dumux::Concept::GridVariables<GridVariables> && !Dumux::Concept::FVGridVariables<GridVariables>;
    }

    //! export the type of the point source context for domain i, see makePointSources
    template<std::size_t id> using PointSourceContext = EmbeddedCoupling::PointSourceContext<GlobalPosition>;
    //! export stencil types
    template<std::size_t id> using CouplingStencils = std::unordered_map<GridIndex<id>, CouplingStencil<id>>;

    /*!
    * \brief call this after grid adaption
    */
    void updateAfterGridAdaption(std::shared_ptr<const GridGeometry<bulkIdx>> bulkGridGeometry,
                                 std::shared_ptr<const GridGeometry<lowDimIdx>> lowDimGridGeometry)
    {
        glue_ = std::make_shared<GlueType>();
        distributed_ = false;
    }

    /*!
     * \brief Constructor
     */
    EmbeddedCouplingManagerBase(std::shared_ptr<const GridGeometry<bulkIdx>> bulkGridGeometry,
                                std::shared_ptr<const GridGeometry<lowDimIdx>> lowDimGridGeometry)
    {
        updateAfterGridAdaption(bulkGridGeometry, lowDimGridGeometry);
    }

    /*!
     * \brief Methods to be accessed by main
     */
    // \{

    /*!
     * \brief Initialize the coupling: compute the point sources
     * \note On distributed grids this is collective over the processes of the bulk grid.
     */
    void init(std::shared_ptr<Problem<bulkIdx>> bulkProblem,
              std::shared_ptr<Problem<lowDimIdx>> lowDimProblem,
              const SolutionVector& curSol)
    {
        ParentType::updateSolution(curSol);
        this->setSubProblems(std::make_tuple(bulkProblem, lowDimProblem));

        integrationOrder_ = getParam<int>("MixedDimension.IntegrationOrder", 1);
        asImp_().computePointSourceData(integrationOrder_);
    }

    /*!
     * \brief Update the solution the coupling terms are evaluated with
     * \note On distributed grids this is collective over the processes of the bulk grid: the processes
     *       exchange the parts of the coupling data that depend on degrees of freedom of other processes.
     */
    void updateSolution(const SolutionVector& curSol)
    {
        ParentType::updateSolution(curSol);
        if (distributed_)
            exchangeDistributedValues_();
    }

    /*!
     * \brief Set values of the lower-dimensional elements that the coupling terms need
     * \param numValues the number of values per element
     * \param f the function f(lowDimElementIdx, values) writing the values of an element into the
     *        std::span<Scalar> values of size numValues
     * \note Call this before init. The values are read with lowDimElementData(id, i), also on processes
     *       that do not hold the lower-dimensional element of the point source id: on distributed grids they
     *       are sent along with the solution whenever the solution is updated.
     */
    template<class Function>
    void setLowDimElementData(std::size_t numValues, Function&& f)
    {
        if (userElementValues_ != invalidIndex)
            DUNE_THROW(Dune::InvalidStateException, "The data of the lower-dimensional elements can only be set once");
        userElementValues_ = addLowDimElementValues(numValues, std::forward<Function>(f));
    }

    //! The value i set with setLowDimElementData of the lower-dimensional element of the point source id
    Scalar lowDimElementData(std::size_t id, std::size_t i) const
    {
        if (userElementValues_ == invalidIndex)
            DUNE_THROW(Dune::InvalidStateException, "No data of the lower-dimensional elements was set");
        return lowDimElementValue(id, userElementValues_ + i);
    }

    // \}

    /*!
     * \brief Methods to be accessed by the assembly
     */
    // \{

    /*!
     * \brief returns an iterable container of all indices of degrees of freedom of domain j
     *        that couple with / influence the element residual of the given element of domain i
     *
     * \param domainI the domain index of domain i
     * \param element the coupled element of domain í
     * \param domainJ the domain index of domain j
     *
     * \note  The element residual definition depends on the discretization scheme of domain i
     *        box: a container of the residuals of all sub control volumes
     *        cc : the residual of the (sub) control volume
     *        fem: the residual of the element
     * \note  This function has to be implemented by all coupling managers for all combinations of i and j
     */
    template<std::size_t i, std::size_t j>
    const CouplingStencil<j>& couplingStencil(Dune::index_constant<i> domainI,
                                              const Element<i>& element,
                                              Dune::index_constant<j> domainJ) const
    {
        static_assert(i != j, "A domain cannot be coupled to itself!");

        const auto eIdx = gridDiscretization(this->problem(domainI)).elementMapper().index(element);
        if (couplingStencils(domainI).count(eIdx))
            return couplingStencils(domainI).at(eIdx);
        else
            return emptyStencil(domainI);
    }

    /*!
     * \brief evaluates the element residual of a coupled element of domain i which depends on the variables
     *        at the degree of freedom with index dofIdxGlobalJ of domain j
     *
     * \param domainI the domain index of domain i
     * \param localAssemblerI the local assembler assembling the element residual of an element of domain i
     * \param domainJ the domain index of domain j
     * \param dofIdxGlobalJ the index of the degree of freedom of domain j which has an influence on the element residual of domain i
     *
     * \note  we only need to evaluate the source contribution to the residual here as the coupling term is the source
     * \return the element residual
     */
    template<std::size_t i, std::size_t j, class LocalAssemblerI>
    decltype(auto) evalCouplingResidual(Dune::index_constant<i> domainI,
                                        const LocalAssemblerI& localAssemblerI,
                                        Dune::index_constant<j> domainJ,
                                        std::size_t dofIdxGlobalJ)
    {
        static_assert(i != j, "A domain cannot be coupled to itself!");

        if constexpr (requires { localAssemblerI.curElemVars(); })
            return localAssemblerI.evalLocalSourceResidual(localAssemblerI.element(), localAssemblerI.curElemVars());
        else
        {
            typename LocalAssemblerI::LocalResidual::ElementResidualVector residual;

            const auto& element = localAssemblerI.element();
            const auto& fvGeometry = localAssemblerI.fvGeometry();
            static_assert(!Detail::LocalDofs::hasNonCVLocalDofsInterface<std::decay_t<decltype(fvGeometry)>>(),
                          "Local dofs without a control volume need the element variables interface");
            const auto& curElemVolVars = localAssemblerI.curElemVolVars();

            residual.resize(fvGeometry.numScv());
            for (const auto& scv : scvs(fvGeometry))
            {
                auto couplingSource = this->problem(domainI).scvPointSources(element, fvGeometry, curElemVolVars, scv);
                couplingSource += this->problem(domainI).source(element, fvGeometry, curElemVolVars, scv);
                couplingSource *= -GridGeometry<i>::Extrusion::volume(fvGeometry, scv)*curElemVolVars[scv].extrusionFactor();
                residual[scv.indexInElement()] = couplingSource;
            }
            return residual;
        }
    }

    /*!
     * \brief Evaluate the derivatives of the element residual with respect to the coupling values whose
     *        terms come from other processes, see forEachRemoteJacobianEntry
     * \note Does nothing on grids that are not distributed. The local assemblers call this for every element.
     */
    template<std::size_t i, class LocalAssemblerI, class JacobianMatrixDiagBlock, class GridVariables>
    void evalAdditionalDomainDerivatives(Dune::index_constant<i> domainI,
                                         const LocalAssemblerI& localAssemblerI,
                                         const typename LocalAssemblerI::LocalResidual::ElementResidualVector&,
                                         JacobianMatrixDiagBlock&,
                                         GridVariables&)
    { evalRemoteDerivatives(domainI, localAssemblerI); }

    /*!
     * \brief Evaluate the derivatives of the element residual with respect to the coupling values whose
     *        terms come from other processes
     *
     * A coupling value is a sum of the solution at degrees of freedom times weights. The terms of other
     * processes are received in updateSolution and do not change during the assembly, so the local
     * matrices have no entries for their degrees of freedom. This computes the derivative of the source
     * residual of the element with respect to each such value by numeric differentiation.
     */
    template<std::size_t i, class LocalAssemblerI>
    void evalRemoteDerivatives(Dune::index_constant<i> domainI, const LocalAssemblerI& localAssemblerI)
    {
        if (!distributed_)
            return;

        const auto& element = localAssemblerI.element();
        const auto eIdx = gridDiscretization(this->problem(domainI)).elementMapper().index(element);
        const auto& offsets = remoteSumOffsets_[i];
        if (offsets[eIdx] == offsets[eIdx+1])
            return;

        const auto origResidual = localAssemblerI.evalLocalSourceResidual(element);
        auto partialDerivs = origResidual;
        constexpr auto numEq = PrimaryVariables<i>::size();
        static const int numDiffMethod = getParamFromGroup<int>(localAssemblerI.problem().paramGroup(), "Assembly.NumericDifferenceMethod");

        for (std::size_t m = offsets[eIdx]; m < offsets[eIdx+1]; ++m)
        {
            const auto sum = remoteSums_[i][m];
            const auto [kind, s] = sumSamples_[sum];
            Scalar* derivatives = remoteDerivatives_[i].data() + remoteDerivativeOffsets_[i][m];

            const auto differentiate = [&](auto domainJ)
            {
                using Values = PrimaryVariables<domainJ>;
                constexpr auto numComponents = Values::size();
                static const auto eps = this->numericEpsilon(domainJ, this->problem(domainJ).paramGroup());
                const auto value = sampleValue_<Values>(samples_[kind][s], this->curSol(domainJ));
                for (std::size_t c = 0; c < numComponents; ++c)
                {
                    partialDerivs = 0.0;
                    const auto evalResiduals = [&](Scalar deflected)
                    {
                        sumDeflection_ = SumDeflection{this, sum, c, deflected - value[c]};
                        return localAssemblerI.evalLocalSourceResidual(element);
                    };
                    NumericDifferentiation::partialDerivative(evalResiduals, value[c], partialDerivs, origResidual,
                                                              eps(value[c], c), numDiffMethod);
                    sumDeflection_ = SumDeflection{};

                    for (std::size_t k = 0; k < partialDerivs.size(); ++k)
                        for (std::size_t eqIdx = 0; eqIdx < numEq; ++eqIdx)
                            derivatives[(k*numEq + eqIdx)*numComponents + c] = partialDerivs[k][eqIdx];
                }
            };

            if (kind == bulkSample)
                differentiate(bulkIdx);
            else
                differentiate(lowDimIdx);
        }

        zeroConstrainedDerivatives_(domainI, localAssemblerI);
    }

    /*!
     * \brief Receive the columns of the Jacobian entries whose column is a degree of freedom of another process
     * \param globalIndex the function globalIndex(domainIdx, dofIdx) returning a global index (std::int64_t) of
     *        the first primary variable of a degree of freedom of this process; the primary variable pvIdx has
     *        the index globalIndex(domainIdx, dofIdx) + pvIdx
     * \note Collective. Call it again if the global indices change.
     */
    template<class GlobalIndex>
    void updateRemoteColumns(GlobalIndex&& globalIndex)
    {
        if (!distributed_)
            return;

        remoteColumns_ = distributedSums_.template exchangeLists<RemoteColumn>([&](std::size_t sum)
        {
            const auto [kind, s] = sumSamples_[sum];
            const auto& sample = samples_[kind][s];
            std::vector<RemoteColumn> columns(sample.dofs.size());
            for (std::size_t j = 0; j < sample.dofs.size(); ++j)
            {
                const std::int64_t index = kind == bulkSample ? globalIndex(bulkIdx, sample.dofs[j])
                                                              : globalIndex(lowDimIdx, sample.dofs[j]);
                columns[j] = RemoteColumn{index, sample.weights[j]};
            }
            return columns;
        });
    }

    /*!
     * \brief Call f(domainI, dofIdx, eqIdx, column, value) for every Jacobian entry whose column is a
     *        degree of freedom of another process
     *
     * The row is the equation eqIdx of the degree of freedom dofIdx of domain domainI on this process, the
     * column a global index as given to updateRemoteColumns. Entries are given in the same order as long as
     * the point sources do not change; a column can appear several times in a row, the values are to be added.
     * The values are those of the last assembly of the Jacobian.
     */
    template<class Function>
    void forEachRemoteJacobianEntry(Function&& f) const
    {
        if (!distributed_)
            return;

        if (remoteColumns_.size() != distributedSums_.size())
            DUNE_THROW(Dune::InvalidStateException, "Call updateRemoteColumns before forEachRemoteJacobianEntry");

        forEachRemoteDerivative_([&](auto domainI, std::size_t dofIdx, std::size_t eqIdx, std::size_t sum, std::size_t c, Scalar derivative)
        {
            for (const auto& column : remoteColumns_[sum])
                f(domainI, dofIdx, eqIdx, column.index + static_cast<std::int64_t>(c), derivative*column.weight);
        });
    }

    /*!
     * \brief Add the product of the Jacobian entries whose column is a degree of freedom of another process
     *        with x to y
     * \note Collective. Every process uses its values of x of the degrees of freedom it sends terms of.
     */
    void addRemoteJacobianProduct(const SolutionVector& x, SolutionVector& y)
    {
        if (!distributed_)
            return;

        distributedSums_.exchange([&](std::size_t sum, std::span<Scalar> values)
        {
            const auto [kind, s] = sumSamples_[sum];
            if (kind == bulkSample)
                addLocalTerms_(samples_[kind][s], x[bulkIdx], values, PrimaryVariables<bulkIdx>::size());
            else
                addLocalTerms_(samples_[kind][s], x[lowDimIdx], values, PrimaryVariables<lowDimIdx>::size());
        }, remoteProductValues_);

        forEachRemoteDerivative_([&](auto domainI, std::size_t dofIdx, std::size_t eqIdx, std::size_t sum, std::size_t c, Scalar derivative)
        { y[domainI][dofIdx][eqIdx] += derivative*distributedSums_.received(sum, remoteProductValues_)[c]; });
    }

    // \}

    /* \brief Compute integration point point sources and associated data
     *
     * This method uses grid glue to intersect the given grids. Over each intersection
     * we later need to integrate a source term. This method places point sources
     * at each quadrature point and provides the point source with the necessary
     * information to compute integrals (quadrature weight and integration element)
     * \param order The order of the quadrature rule for integration of sources over an intersection
     * \note The second parameter has no effect
     */
    void computePointSourceData(std::size_t order = 1, bool /*verbose*/ = false)
    {
        if (gridView(bulkIdx).comm().size() > 1)
        {
            computeDistributedPointSourceData_(order);
            return;
        }

        // initialize the maps
        // do some logging and profiling
        Dune::Timer watch;
        std::cout << "Initializing the point sources..." << std::endl;

        // clear all internal members like pointsource vectors and stencils
        // initializes the point source id counter
        clear();

        // precompute the dof indices for efficiency for the control-volume finite element methods
        this->precomputeElementDofIndices(bulkIdx);
        this->precomputeElementDofIndices(lowDimIdx);

        const auto& bulkGridGeometry = gridDiscretization(this->problem(bulkIdx));
        const auto& lowDimGridGeometry = gridDiscretization(this->problem(lowDimIdx));

        // intersect the bounding box trees
        glueGrids();

        pointSourceData_.reserve(this->glue().size());
        averageDistanceToBulkCell_.reserve(this->glue().size());
        for (const auto& is : intersections(this->glue()))
        {
            // all inside elements are identical...
            const auto& inside = is.targetEntity(0);
            // get the intersection geometry for integrating over it
            const auto intersectionGeometry = is.geometry();

            // get the Gaussian quadrature rule for the local intersection
            const auto& quad = Dune::QuadratureRules<Scalar, lowDimDim>::rule(intersectionGeometry.type(), order);
            const std::size_t lowDimElementIdx = lowDimGridGeometry.elementMapper().index(inside);

            // iterate over all quadrature points
            for (auto&& qp : quad)
            {
                // compute the coupling stencils
                for (std::size_t outsideIdx = 0; outsideIdx < is.numDomainNeighbors(); ++outsideIdx)
                {
                    const auto& outside = is.domainEntity(outsideIdx);
                    const std::size_t bulkElementIdx = bulkGridGeometry.elementMapper().index(outside);

                    // each quadrature point will be a point source for the sub problem
                    const auto globalPos = intersectionGeometry.global(qp.position());
                    const auto id = idCounter_++;
                    const auto qpweight = qp.weight();
                    const auto ie = intersectionGeometry.integrationElement(qp.position());
                    pointSources(bulkIdx).emplace_back(globalPos, id, qpweight, ie, bulkElementIdx);
                    pointSources(bulkIdx).back().setEmbeddings(is.numDomainNeighbors());
                    pointSources(lowDimIdx).emplace_back(globalPos, id, qpweight, ie, lowDimElementIdx);
                    pointSources(lowDimIdx).back().setEmbeddings(is.numDomainNeighbors());

                    // pre compute additional data used for the evaluation of
                    // the actual solution dependent source term
                    PointSourceData psData;

                    if constexpr (isCVFE<lowDimIdx>())
                    {
                        using ShapeValues = std::vector<Dune::FieldVector<Scalar, 1> >;
                        const auto lowDimGeometry = gridDiscretization(this->problem(lowDimIdx)).element(lowDimElementIdx).geometry();
                        ShapeValues shapeValues;
                        this->getShapeValues(lowDimIdx, gridDiscretization(this->problem(lowDimIdx)), lowDimGeometry, globalPos, shapeValues);
                        psData.addLowDimInterpolation(shapeValues, this->elementDofIndices(lowDimIdx, lowDimElementIdx), lowDimElementIdx);
                    }
                    else
                    {
                        psData.addLowDimInterpolation(lowDimElementIdx);
                    }

                    // add data needed to compute integral over the circle
                    if constexpr (isCVFE<bulkIdx>())
                    {
                        using ShapeValues = std::vector<Dune::FieldVector<Scalar, 1> >;
                        const auto bulkGeometry = gridDiscretization(this->problem(bulkIdx)).element(bulkElementIdx).geometry();
                        ShapeValues shapeValues;
                        this->getShapeValues(bulkIdx, gridDiscretization(this->problem(bulkIdx)), bulkGeometry, globalPos, shapeValues);
                        psData.addBulkInterpolation(shapeValues, this->elementDofIndices(bulkIdx, bulkElementIdx), bulkElementIdx);
                    }
                    else
                    {
                        psData.addBulkInterpolation(bulkElementIdx);
                    }

                    // publish point source data in the global vector
                    this->pointSourceData().emplace_back(std::move(psData));

                    // compute average distance to bulk cell
                    averageDistanceToBulkCell_.push_back(averageDistancePointGeometry(globalPos, outside.geometry()));

                    // export the lowdim coupling stencil
                    // we insert all dofs / elements and make it unique later
                    if (isCVFE<bulkIdx>())
                    {
                        const auto& dofs = this->elementDofIndices(bulkIdx, bulkElementIdx);
                        this->couplingStencils(lowDimIdx)[lowDimElementIdx].insert(this->couplingStencils(lowDimIdx)[lowDimElementIdx].end(),
                                                                                   dofs.begin(), dofs.end());
                    }
                    else
                    {
                        this->couplingStencils(lowDimIdx)[lowDimElementIdx].push_back(bulkElementIdx);
                    }

                    // export the bulk coupling stencil
                    // we insert all dofs / elements and make it unique later
                    if (isCVFE<lowDimIdx>())
                    {
                        const auto& dofs = this->elementDofIndices(lowDimIdx, lowDimElementIdx);
                        this->couplingStencils(bulkIdx)[bulkElementIdx].insert(this->couplingStencils(bulkIdx)[bulkElementIdx].end(),
                                                                               dofs.begin(), dofs.end());

                    }
                    else
                    {
                        this->couplingStencils(bulkIdx)[bulkElementIdx].push_back(lowDimElementIdx);
                    }
                }
            }
        }

        // make stencils unique
        using namespace Dune::Hybrid;
        forEach(integralRange(std::integral_constant<std::size_t, MDTraits::numSubDomains>{}), [&](const auto domainIdx)
        {
            for (auto&& stencil : this->couplingStencils(domainIdx))
            {
                std::sort(stencil.second.begin(), stencil.second.end());
                stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());
            }
        });

        std::cout << "took " << watch.elapsed() << " seconds." << std::endl;
    }

    /*!
     * \brief Methods to be accessed by the subproblems
     */
    // \{

    /*!
     * \brief Return a reference to the pointSource data
     * \note On distributed grids an element index of the data is only valid on the process holding the
     *       element as interior element. Values of the lower-dimensional element that are needed where the
     *       bulk sources are can be passed with setLowDimElementData.
     */
    const PointSourceData& pointSourceData(std::size_t id) const
    { return pointSourceData_[id]; }

    //! Return a reference to the bulk problem
    template<std::size_t id>
    const GridView<id>& gridView(Dune::index_constant<id> domainIdx) const
    { return gridDiscretization(this->problem(domainIdx)).gridView(); }

    //! Return data for a bulk point source with the identifier id
    PrimaryVariables<bulkIdx> bulkPriVars(std::size_t id) const
    {
        if (distributed_)
            return sampleValue_<PrimaryVariables<bulkIdx>>(samples_[bulkSample][idBulkSample_[id]], this->curSol(bulkIdx));
        return pointSourceData_[id].interpolateBulk(this->curSol(bulkIdx));
    }

    //! Return data for a low dim point source with the identifier id
    PrimaryVariables<lowDimIdx> lowDimPriVars(std::size_t id) const
    {
        if (distributed_)
            return sampleValue_<PrimaryVariables<lowDimIdx>>(samples_[lowDimSample][idLowDimSample_[id]], this->curSol(lowDimIdx));
        return pointSourceData_[id].interpolateLowDim(this->curSol(lowDimIdx));
    }

    //! return the average distance to the coupled bulk cell center
    Scalar averageDistance(std::size_t id) const
    { return averageDistanceToBulkCell_[id]; }

    //! Return reference to bulk point sources
    const std::vector<PointSource<bulkIdx>>& bulkPointSources() const
    { return std::get<bulkIdx>(pointSources_); }

    //! Return reference to low dim point sources
    const std::vector<PointSource<lowDimIdx>>& lowDimPointSources() const
    { return std::get<lowDimIdx>(pointSources_); }

    //! Return the point source if domain i
    template<std::size_t i>
    const std::vector<PointSource<i>>& pointSources(Dune::index_constant<i> dom) const
    { return std::get<i>(pointSources_); }

    //! Return reference to bulk coupling stencil member of domain i
    template<std::size_t i>
    const CouplingStencils<i>& couplingStencils(Dune::index_constant<i> dom) const
    { return std::get<i>(couplingStencils_); }

    /*!
     * \brief The contexts of the point sources of domain i
     * \note One context for every point source in an element, see EmbeddedCoupling::PointSourceContext
     */
    template<std::size_t i>
    std::vector<PointSourceContext<i>> pointSourceContexts(Dune::index_constant<i> domainI) const
    {
        std::vector<PointSourceContext<i>> contexts;
        const auto& sources = pointSources(domainI);
        const auto& additional = std::get<i>(additionalPointSourceContexts_);
        contexts.reserve(sources.size() + additional.size());
        for (const auto& source : sources)
            contexts.emplace_back(source.position(), source.id(),
                                  source.quadratureWeight()*source.integrationElement(),
                                  source.elementIndex(), source.embeddings());
        contexts.insert(contexts.end(), additional.begin(), additional.end());
        return contexts;
    }

    /*!
     * \brief The point sources of domain i for a problem using the element discretization interface
     *
     * \param domainI the domain index
     * \param f the evaluation function `f(elemDisc, elemVars, context)` returning the source values of a context
     *        for the problem to return from `pointSources()`. The context is
     *        an `EmbeddedCoupling::PointSourceContext` with the identifier of the point source, to access the
     *        data of the coupling, and its weight. The type of the context depends on the local dof, so the
     *        function takes it as `const auto&`.
     */
    template<std::size_t i, class Function>
    auto makePointSources(Dune::index_constant<i> domainI, Function f) const
    {
        return Experimental::makePointSources(
            gridDiscretization(this->problem(domainI)), pointSourceContexts(domainI), std::move(f)
        );
    }

    //! Return reference to point source data vector member
    const std::vector<PointSourceData>& pointSourceData() const
    { return pointSourceData_; }

    //! Return a reference to an empty stencil
    template<std::size_t i>
    const CouplingStencil<i>& emptyStencil(Dune::index_constant<i> dom) const
    { return std::get<i>(emptyStencil_); }

    /*!
     * \brief The global indices of the degrees of freedom of an element of domain i
     *        for the control-volume finite element methods, in the order of the basis functions
     */
    template<std::size_t i>
    const std::vector<GridIndex<i>>& elementDofIndices(Dune::index_constant<i> dom, GridIndex<i> eIdx) const
    { return std::get<i>(elementDofIndices_)[eIdx]; }

    /*!
     * \brief the solution vector of the subproblem
     * \param domainIdx The domain index
     * \note in case of numeric differentiation the solution vector always carries the deflected solution
     */
    template<std::size_t i>
    auto& curSol(Dune::index_constant<i> domainIdx)
    { return ParentType::curSol(domainIdx); }

    /*!
     * \brief the solution vector of the subproblem
     * \param domainIdx The domain index
     * \note in case of numeric differentiation the solution vector always carries the deflected solution
     */
    template<std::size_t i>
    const auto& curSol(Dune::index_constant<i> domainIdx) const
    { return ParentType::curSol(domainIdx); }

    /*!
     * \brief Compute colors for multithreaded assembly
     */
    void computeColorsForAssembly()
    {
        Dune::Timer timer;

        // start with elements connected within the subdomains
        const auto& bulkGG = gridDiscretization(asImp_().problem(bulkIdx));
        const auto& lowDimGG = gridDiscretization(asImp_().problem(lowDimIdx));

        auto connectedElementsBulk = Detail::computeConnectedElements(bulkGG);
        auto connectedElementsLowDim = Detail::computeConnectedElements(lowDimGG);

        // to the elements from the own domain, we have to add the extended source stencil (non-local couplings)
        // we consider this only in the bulk domain: the elements whose extended source stencil contains a dof
        std::vector<std::vector<std::size_t>> extendedSourceElements(bulkGG.numDofs());
        for (const auto& element : elements(bulkGG.gridView()))
        {
            const auto eIdx = bulkGG.elementMapper().index(element);
            for (auto dofIdx : asImp_().extendedSourceStencil(eIdx))
                extendedSourceElements[dofIdx].push_back(eIdx);
        }

        // and the elements reading a dof, connectedElementsBulk is indexed by vertex for
        // the control-volume finite element methods
        std::vector<std::vector<std::size_t>> dofElementsBulk;
        if constexpr (isCVFE<bulkIdx>())
        {
            dofElementsBulk.resize(bulkGG.numDofs());
            for (const auto& element : elements(bulkGG.gridView()))
            {
                const auto eIdx = bulkGG.elementMapper().index(element);
                for (const auto dofIdx : elementDofIndices(bulkIdx, eIdx))
                    dofElementsBulk[dofIdx].push_back(eIdx);
            }
        }
        const auto& readingElementsBulk = isCVFE<bulkIdx>() ? dofElementsBulk : connectedElementsBulk;

        // then we also need a similar container that tells us
        // which elements _from the other domain_ modify our dofs.
        // That, we can obtain via the coupling stencils
        std::array<std::vector<std::vector<std::size_t>>, 2> connectedElementsCoupling;

        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::index_constant<MDTraits::numSubDomains>{}), [&](const auto i)
        {
            connectedElementsCoupling[i()].resize(gridDiscretization(asImp_().problem(i)).numDofs());
            forEach(integralRange(Dune::index_constant<MDTraits::numSubDomains>{}), [&](const auto j)
            {
                if constexpr (i != j)
                {
                    for (const auto& element : elements(gridDiscretization(asImp_().problem(j)).gridView()))
                    {
                        const auto eIdx = gridDiscretization(asImp_().problem(j)).elementMapper().index(element);
                        for (auto dofIdx : asImp_().couplingStencil(j, element, i))
                        {
                            const auto& elems = connectedElementsCoupling[i()][dofIdx];
                            if (std::find(elems.begin(), elems.end(), eIdx) == elems.end())
                                connectedElementsCoupling[i()][dofIdx].push_back(eIdx);
                        }
                    }
                }
            });
        });

        forEach(integralRange(Dune::index_constant<MDTraits::numSubDomains>{}), [&](const auto i)
        {
            const auto& gg = gridDiscretization(asImp_().problem(i));

            std::vector<int> colors(gg.gridView().size(0), -1);

            // pre-reserve some memory for helper arrays to avoid reallocation
            std::vector<int> neighborColors; neighborColors.reserve(200);
            std::vector<bool> colorUsed; colorUsed.reserve(200);

            auto& elementSets = std::get<i>(elementSets_);
            elementSets.clear();
            const auto& connectedElements = std::get<i>(std::tie(connectedElementsBulk, connectedElementsLowDim));

            for (const auto& element : elements(gg.gridView()))
            {
                // compute neighbor colors based on discretization-dependent stencil
                neighborColors.clear();
                Detail::addNeighborColors(gg, element, colors, connectedElements, neighborColors);

                if constexpr (i == bulkIdx)
                {
                    // make sure we also add all elements in the extended source stencil
                    // Detail::addNeighborColors only adds direct neighbors
                    // add everyone that also reads or modifies the extended source stencil dofs
                    const auto eIdx = bulkGG.elementMapper().index(element);
                    const auto addModifyingElementColors = [&](const auto dofIdx)
                    {
                        for (const auto nIdx : extendedSourceElements[dofIdx])
                            neighborColors.push_back(colors[nIdx]);
                    };
                    for (const auto dofIdx : asImp_().extendedSourceStencil(eIdx))
                    {
                        // the element of a cell-centered dof reads it too
                        if constexpr (!isCVFE<bulkIdx>())
                            neighborColors.push_back(colors[dofIdx]);
                        for (const auto nIdx : readingElementsBulk[dofIdx])
                            neighborColors.push_back(colors[nIdx]);
                        addModifyingElementColors(dofIdx);
                    }

                    // add everyone whose extended source stencil modifies the dofs this element reads
                    if constexpr (isCVFE<bulkIdx>())
                    {
                        for (const auto dofIdx : elementDofIndices(bulkIdx, eIdx))
                            addModifyingElementColors(dofIdx);
                    }
                    else
                    {
                        addModifyingElementColors(eIdx);
                        for (const auto nIdx : connectedElementsBulk[eIdx])
                            addModifyingElementColors(nIdx);
                    }
                }

                // add neighbor colors based on coupling stencil
                forEach(integralRange(Dune::index_constant<MDTraits::numSubDomains>{}), [&](const auto j)
                {
                    if constexpr (i != j)
                    {
                        // add all elements that manipulate the same coupled dofs
                        for (auto dofIdx : asImp_().couplingStencil(i, element, j))
                            for (auto eIdx : connectedElementsCoupling[j][dofIdx])
                                neighborColors.push_back(colors[eIdx]);
                    }
                });

                // find smallest color (positive integer) not in neighborColors
                const auto color = Detail::smallestAvailableColor(neighborColors, colorUsed);

                // assign color to element
                colors[gg.elementMapper().index(element)] = color;

                // add element to the set of elements with the same color
                if (color < elementSets.size())
                    elementSets[color].push_back(element.seed());
                else
                    elementSets.push_back(std::vector<ElementSeed<i>>{ element.seed() });
            }
        });

        std::cout << Fmt::format("Colored in {} seconds:\n", timer.elapsed());
        forEach(integralRange(Dune::index_constant<MDTraits::numSubDomains>{}), [&](const auto i)
        {
            std::cout << Fmt::format(
                "-- {} elements in subdomain {} with {} colors\n",
                gridDiscretization(asImp_().problem(i)).gridView().size(0),
                i(), std::get<i>(elementSets_).size()
            );
        });
    }

    //! The sets of elements of domain i that can be assembled in parallel, as seeds of the elements
    template<std::size_t i>
    const auto& elementSets(Dune::index_constant<i> domainId) const
    { return std::get<i>(elementSets_); }

    /*!
     * \brief Execute assembly kernel in parallel
     *
     * \param domainId the domain index of domain i
     * \param assembleElement kernel function to execute for one element
     */
    template<std::size_t i, class AssembleElementFunc>
    void assembleMultithreaded(Dune::index_constant<i> domainId, AssembleElementFunc&& assembleElement) const
    {
        if (std::get<i>(elementSets_).empty())
            DUNE_THROW(Dune::InvalidStateException, "Call computeColorsForAssembly before assembling in parallel!");

        // make this element loop run in parallel
        // for this we have to color the elements so that we don't get
        // race conditions when writing into the global matrix or modifying grid variable caches
        // each color can be assembled using multiple threads
        const auto& grid = gridDiscretization(this->problem(domainId)).gridView().grid();
        for (const auto& elements : std::get<i>(elementSets_))
        {
            Dumux::parallelFor(elements.size(), [&](const std::size_t n)
            {
                const auto element = grid.entity(elements[n]);
                assembleElement(element);
            });
        }
    }

    /*!
     * \brief Extended source stencil (for the bulk domain)
     * \note The default returns an empty stencil but some of the embedded
     *       coupling managers use this
     */
    const CouplingStencil<bulkIdx>& extendedSourceStencil(std::size_t eIdx) const
    { return asImp_().emptyStencil(bulkIdx); }

protected:

    //! computes the dof indices per element for the control-volume finite element methods
    template<std::size_t id>
    void precomputeElementDofIndices(Dune::index_constant<id> domainIdx)
    {
        if constexpr (isCVFE<domainIdx>())
        {
            auto& dofIndices = std::get<id>(elementDofIndices_);
            dofIndices.resize(gridView(domainIdx).size(0));
            auto elemDisc = localView(gridDiscretization(this->problem(domainIdx)));
            for (const auto& element : elements(gridView(domainIdx)))
            {
                elemDisc.bindElement(element);
                dofIndices[elemDisc.elementIndex()] = EmbeddedCoupling::elementDofIndices(elemDisc);
            }
        }
    }

    //! computes the vertex indices per element for the box method
    template<std::size_t id>
    [[deprecated("Will be removed after release 3.11. Use precomputeElementDofIndices.")]]
    void precomputeVertexIndices(Dune::index_constant<id> domainIdx)
    { precomputeElementDofIndices(domainIdx); }

    //! compute the shape function for a given point and geometry
    template<std::size_t i, class FVGG, class Geometry, class ShapeValues>
    void getShapeValues(Dune::index_constant<i> domainI, const FVGG& gridGeometry, const Geometry& geo, const GlobalPosition& globalPos, ShapeValues& shapeValues)
    {
        if constexpr (DiscretizationMethods::isCVFE<typename FVGG::DiscretizationMethod>)
            EmbeddedCoupling::basisFunctionValues(gridGeometry, geo, globalPos, shapeValues);
        else
            DUNE_THROW(Dune::InvalidStateException, "Shape values requested for a discretization without basis functions!");
    }

    //! Clear all internal data members
    void clear()
    {
        pointSources(bulkIdx).clear();
        pointSources(lowDimIdx).clear();
        couplingStencils(bulkIdx).clear();
        couplingStencils(lowDimIdx).clear();
        std::get<bulkIdx>(elementDofIndices_).clear();
        std::get<lowDimIdx>(elementDofIndices_).clear();
        pointSourceData_.clear();
        averageDistanceToBulkCell_.clear();
        std::get<bulkIdx>(additionalPointSourceContexts_).clear();
        std::get<lowDimIdx>(additionalPointSourceContexts_).clear();

        distributed_ = false;
        distributedSums_ = DistributedSums<Scalar>();
        for (auto& samples : samples_)
            samples.clear();
        for (auto& indices : sampleIndices_)
            indices.clear();
        sumSamples_.clear();
        idBulkSample_.clear();
        idLowDimSample_.clear();
        for (auto& v : remoteSumOffsets_) v.clear();
        for (auto& v : remoteSums_) v.clear();
        for (auto& v : remoteDerivativeOffsets_) v.clear();
        for (auto& v : remoteDerivatives_) v.clear();
        remoteColumns_.clear();
        pointSourceSampleKeys_.clear();
        importedPointSources_.clear();
        numOwnedSamples_ = 0;

        idCounter_ = 0;
    }

    //! compute the intersections between the two grids
    void glueGrids()
    {
        const auto& bulkGridGeometry = gridDiscretization(this->problem(bulkIdx));
        const auto& lowDimGridGeometry = gridDiscretization(this->problem(lowDimIdx));

        // intersect the bounding box trees
        glue_->build(bulkGridGeometry.boundingBoxTree(), lowDimGridGeometry.boundingBoxTree());
    }

    //! Return reference to point source data vector member
    std::vector<PointSourceData>& pointSourceData()
    { return pointSourceData_; }

    //! Return the contexts of point sources not in the point sources of domain i, e.g. distributed sources
    template<std::size_t i>
    std::vector<PointSourceContext<i>>& additionalPointSourceContexts(Dune::index_constant<i> dom)
    { return std::get<i>(additionalPointSourceContexts_); }

    //! Return reference to average distances to bulk cell
    std::vector<Scalar>& averageDistanceToBulkCell()
    { return averageDistanceToBulkCell_; }

    //! Return the point source if domain i
    template<std::size_t i>
    std::vector<PointSource<i>>& pointSources(Dune::index_constant<i> dom)
    { return std::get<i>(pointSources_); }

    //! Return reference to bulk coupling stencil member of domain i
    template<std::size_t i>
    CouplingStencils<i>& couplingStencils(Dune::index_constant<i> dom)
    { return std::get<i>(couplingStencils_); }

    //! Return a reference to the vertex indices
    template<std::size_t i>
    [[deprecated("Will be removed after release 3.11. Use elementDofIndices.")]]
    std::vector<GridIndex<i>>& vertexIndices(Dune::index_constant<i> dom, GridIndex<i> eIdx)
    { return std::get<i>(elementDofIndices_)[eIdx]; }

    //! Return a reference to the vertex indices container
    template<std::size_t i>
    [[deprecated("Will be removed after release 3.11. Use elementDofIndices.")]]
    std::vector<std::vector<GridIndex<i>>>& vertexIndices(Dune::index_constant<i> dom)
    { return std::get<i>(elementDofIndices_); }

    const GlueType& glue() const
    { return *glue_; }

    /*!
     * \name Point sources on distributed grids
     *
     * The owner of a lower-dimensional element creates the point sources of the element and numbers them,
     * so a point source is identified by the key (owner rank, number). A process holds the point sources
     * whose lower-dimensional element it owns and those with a bulk source in one of its interior bulk
     * elements, under process-local ids. The values of a point source are samples of the solution, sums of
     * terms on the processes holding the sampled elements; the processes holding a point source receive
     * the terms of the other processes whenever the solution is updated.
     */
    // \{

    //! the key of a point source or a sample: the rank of the owner and the number on the owner
    using DistributedKey = typename DistributedSums<Scalar>::Key;

    //! the two kinds of samples
    enum SampleKind : int { bulkSample = 0, lowDimSample = 1 };

    //! a term of a sample: the value at a point of an element of the given rank, times a weight
    struct SampleTerm
    {
        int rank;
        std::size_t element;
        GlobalPosition point;
        Scalar weight;
    };

    //! a sample as created by the owner of the lower-dimensional element
    struct SampleSpec
    {
        SampleKind kind;
        DistributedKey key;
        std::vector<int> contributors;
        std::vector<int> consumers;
        std::vector<SampleTerm> terms;
    };

    //! a part of a lower-dimensional element of this process inside one or more bulk elements of any process
    struct DistributedPiece
    {
        GridIndex<lowDimIdx> lowDimElementIdx;
        std::vector<GlobalPosition> corners;
        std::vector<std::pair<int, std::size_t>> bulkElements;
        std::vector<BulkGeometryCopy> bulkGeometries;

        PieceGeometry geometry() const
        { return PieceGeometry(Dune::GeometryTypes::simplex(lowDimDim), corners); }
    };

    //! whether the point sources are set up for distributed grids
    bool isDistributed() const
    { return distributed_; }

    //! the communication of the bulk grid
    decltype(auto) bulkComm() const
    { return gridView(bulkIdx).comm(); }

    /*!
     * \brief Start the setup on distributed grids: clear the data and compute the pieces
     * \return the pieces of the lower-dimensional elements of this process, see computeDistributedPieces
     * \note Collective.
     */
    std::vector<DistributedPiece> beginDistributedSetup()
    {
        if constexpr (lowDimDim != 1)
            DUNE_THROW(Dune::NotImplemented, "Embedded coupling on distributed grids for lower-dimensional domains of dimension " << lowDimDim);

        if (gridView(lowDimIdx).comm().size() != bulkComm().size())
            DUNE_THROW(Dune::InvalidStateException, "The bulk and the lower-dimensional grid have to be distributed over the same processes,"
                       << " but they are distributed over " << bulkComm().size() << " and " << gridView(lowDimIdx).comm().size() << " processes");

        clear();
        this->precomputeElementDofIndices(bulkIdx);
        this->precomputeElementDofIndices(lowDimIdx);
        distributed_ = true;
        distributedSums_ = DistributedSums<Scalar>(bulkComm());

        return computeDistributedPieces();
    }

    /*!
     * \brief The parts of the interior lower-dimensional elements of this process inside the bulk elements
     *
     * The intersections are computed on the processes of the bulk elements, as the sequential glue
     * computes them, and sent to the owner of the lower-dimensional element. Intersections of one
     * lower-dimensional element with identical corners, e.g. a segment on a face shared by two bulk
     * elements, are merged into one piece with several bulk elements, as in the sequential glue.
     * \note Collective.
     */
    std::vector<DistributedPiece> computeDistributedPieces() const
    {
        using namespace EmbeddedCoupling::Detail;
        const auto& comm = bulkComm();

        std::vector<std::vector<char>> messages(comm.size());
        {
            const auto glue = makeDistributedGlue(gridDiscretization(this->problem(bulkIdx)),
                                                  gridDiscretization(this->problem(lowDimIdx)));
            for (const auto& is : intersections(glue))
            {
                MessageWriter writer(messages[is.targetRank()]);
                writer.put(is.targetIndex());
                writer.put(is.domainIndex());
                const auto geometry = is.geometry();
                writer.put(geometry.corners());
                for (int c = 0; c < geometry.corners(); ++c)
                    writer.put(geometry.corner(c));
                const auto bulkGeometry = is.domainEntity().geometry();
                writer.put(bulkGeometry.type().id());
                writer.put(bulkGeometry.corners());
                for (int c = 0; c < bulkGeometry.corners(); ++c)
                    writer.put(bulkGeometry.corner(c));
            }
        }

        struct RawPiece
        {
            std::size_t lowDimElementIdx;
            int rank;
            std::size_t bulkElementIdx;
            std::vector<GlobalPosition> corners;
            BulkGeometryCopy bulkGeometry;
        };

        std::vector<RawPiece> rawPieces;
        for (const auto& [sender, message] : exchangeMessages(comm, messages))
        {
            MessageReader reader(message);
            while (!reader.done())
            {
                const auto lowDimElementIdx = reader.get<std::size_t>();
                const auto bulkElementIdx = reader.get<std::size_t>();
                std::vector<GlobalPosition> corners(reader.get<int>());
                for (auto& corner : corners)
                    corner = reader.get<GlobalPosition>();
                const auto topologyId = reader.get<unsigned int>();
                std::vector<GlobalPosition> bulkCorners(reader.get<int>());
                for (auto& corner : bulkCorners)
                    corner = reader.get<GlobalPosition>();
                rawPieces.push_back(RawPiece{
                    lowDimElementIdx, sender, bulkElementIdx, std::move(corners),
                    BulkGeometryCopy(Dune::GeometryType(topologyId, bulkDim), bulkCorners)
                });
            }
        }

        std::sort(rawPieces.begin(), rawPieces.end(), [](const auto& a, const auto& b)
        {
            return std::tie(a.lowDimElementIdx, a.rank, a.bulkElementIdx)
                 < std::tie(b.lowDimElementIdx, b.rank, b.bulkElementIdx);
        });

        std::vector<DistributedPiece> pieces;
        std::size_t firstPieceOfElement = 0;
        for (const auto& raw : rawPieces)
        {
            if (pieces.empty() || pieces.back().lowDimElementIdx != raw.lowDimElementIdx)
                firstPieceOfElement = pieces.size();

            auto piece = std::find_if(pieces.begin() + firstPieceOfElement, pieces.end(),
                                      [&](const auto& p){ return cornersMatch_(raw.corners, p.corners); });
            if (piece == pieces.end())
            {
                pieces.push_back(DistributedPiece{static_cast<GridIndex<lowDimIdx>>(raw.lowDimElementIdx), raw.corners, {}, {}});
                piece = std::prev(pieces.end());
            }

            piece->bulkElements.emplace_back(raw.rank, raw.bulkElementIdx);
            piece->bulkGeometries.push_back(raw.bulkGeometry);
        }

        return pieces;
    }

    //! A new point source owned by this process, with the lower-dimensional interpolation at the given point
    std::size_t addOwnedPointSource(GridIndex<lowDimIdx> lowDimElementIdx, const GlobalPosition& point)
    {
        const auto id = addPointSourceId_();
        addInterpolation(lowDimIdx, pointSourceData_[id], lowDimElementIdx, point);
        return id;
    }

    //! The key of a point source owned by this process
    DistributedKey ownedPointSourceKey(std::size_t id) const
    { return {bulkComm().rank(), id}; }

    //! The process-local id of the point source with the given key, added if this process does not hold it yet
    std::size_t pointSourceId(const DistributedKey& key)
    {
        if (key.first == bulkComm().rank())
            return key.second;

        const auto [it, inserted] = importedPointSources_.try_emplace(key, idCounter_);
        if (inserted)
            addPointSourceId_();
        return it->second;
    }

    //! The keys of the samples of the bulk and the lower-dimensional values of point source id
    void setPointSourceSamples(std::size_t id, const DistributedKey& bulkKey, const DistributedKey& lowDimKey)
    { pointSourceSampleKeys_[id] = {bulkKey, lowDimKey}; }

    //! Set the average distance of point source id to its bulk element
    void setAverageDistance(std::size_t id, Scalar distance)
    { averageDistanceToBulkCell_[id] = distance; }

    //! Add the interpolation of domain i at a point of an element to the point source data
    template<std::size_t i>
    void addInterpolation(Dune::index_constant<i> domainI, PointSourceData& data, std::size_t eIdx, const GlobalPosition& point)
    {
        if constexpr (isCVFE<i>())
        {
            const auto& gg = gridDiscretization(this->problem(domainI));
            std::vector<Dune::FieldVector<Scalar, 1>> shapeValues;
            this->getShapeValues(domainI, gg, gg.element(eIdx).geometry(), point, shapeValues);
            if constexpr (i == bulkIdx)
                data.addBulkInterpolation(shapeValues, this->elementDofIndices(domainI, eIdx), eIdx);
            else
                data.addLowDimInterpolation(shapeValues, this->elementDofIndices(domainI, eIdx), eIdx);
        }
        else
        {
            if constexpr (i == bulkIdx)
                data.addBulkInterpolation(eIdx);
            else
                data.addLowDimInterpolation(eIdx);
        }
    }

    //! Add the dofs of an element of domain j to the coupling stencil of an element of domain i
    template<std::size_t i, std::size_t j>
    void addToCouplingStencil(Dune::index_constant<i> domainI, std::size_t eIdxI,
                               Dune::index_constant<j> domainJ, std::size_t eIdxJ)
    {
        auto& stencil = this->couplingStencils(domainI)[eIdxI];
        if constexpr (isCVFE<j>())
        {
            const auto& dofs = this->elementDofIndices(domainJ, eIdxJ);
            stencil.insert(stencil.end(), dofs.begin(), dofs.end());
        }
        else
            stencil.push_back(eIdxJ);
    }

    //! A new key for a sample of the given kind created by this process
    DistributedKey newSampleKey(SampleKind kind)
    { return {bulkComm().rank(), 2*numOwnedSamples_++ + kind}; }

    //! Add a term to a sample; for a cell-centered domain the point does not matter, so terms of the same element are added up
    void addSampleTerm(SampleSpec& spec, const SampleTerm& term)
    {
        const bool cellCentered = spec.kind == bulkSample ? !isCVFE<bulkIdx>() : !isCVFE<lowDimIdx>();
        if (cellCentered)
        {
            const auto it = std::find_if(spec.terms.begin(), spec.terms.end(), [&](const auto& t){
                return t.rank == term.rank && t.element == term.element;
            });
            if (it != spec.terms.end())
            {
                it->weight += term.weight;
                return;
            }
        }
        spec.terms.push_back(term);
    }

    /*!
     * \brief Set up the samples on all processes that contribute terms to them or need their values
     * \param specs the samples created by this process
     * \note Collective.
     */
    void distributeSamples(const std::vector<SampleSpec>& specs)
    {
        using namespace EmbeddedCoupling::Detail;
        const auto& comm = bulkComm();

        std::vector<std::vector<char>> messages(comm.size());
        for (const auto& spec : specs)
        {
            auto ranks = spec.contributors;
            ranks.insert(ranks.end(), spec.consumers.begin(), spec.consumers.end());
            std::sort(ranks.begin(), ranks.end());
            ranks.erase(std::unique(ranks.begin(), ranks.end()), ranks.end());

            for (const int p : ranks)
            {
                MessageWriter writer(messages[p]);
                writer.put(static_cast<int>(spec.kind));
                writer.put(spec.key);
                writer.put(spec.contributors);
                writer.put(spec.consumers);
                writer.put(static_cast<std::size_t>(std::count_if(spec.terms.begin(), spec.terms.end(),
                                                                  [&](const auto& t){ return t.rank == p; })));
                for (const auto& term : spec.terms)
                {
                    if (term.rank != p)
                        continue;
                    writer.put(term.element);
                    writer.put(term.point);
                    writer.put(term.weight);
                }
            }
        }

        const int rank = comm.rank();
        for (const auto& [sender, message] : exchangeMessages(comm, messages))
        {
            MessageReader reader(message);
            while (!reader.done())
            {
                const auto kind = static_cast<SampleKind>(reader.get<int>());
                const auto key = reader.get<DistributedKey>();
                const auto contributors = reader.getVector<int>();
                const auto consumers = reader.getVector<int>();

                const auto s = sampleIndex_(kind, key);
                const auto numTerms = reader.get<std::size_t>();
                for (std::size_t t = 0; t < numTerms; ++t)
                {
                    const auto element = reader.get<std::size_t>();
                    const auto point = reader.get<GlobalPosition>();
                    const auto weight = reader.get<Scalar>();
                    addSampleTerm_(kind, samples_[kind][s], element, point, weight);
                }

                const auto contains = [](const auto& ranks, int r){ return std::find(ranks.begin(), ranks.end(), r) != ranks.end(); };
                const auto hasOther = [&](const auto& ranks){ return std::any_of(ranks.begin(), ranks.end(), [&](int r){ return r != rank; }); };
                const bool contributes = contains(contributors, rank);
                const bool consumes = contains(consumers, rank);
                if ((contributes && hasOther(consumers)) || (consumes && hasOther(contributors)))
                {
                    auto& sample = samples_[kind][s];
                    sample.sum = distributedSums_.insert(key, numSampleComponents_(kind));
                    if (contributes)
                        for (const int p : consumers)
                            distributedSums_.addDestination(sample.sum, p);
                    if (consumes)
                        for (const int p : contributors)
                            distributedSums_.addSource(sample.sum, p);
                }
            }
        }
    }

    /*!
     * \brief Finish the setup on distributed grids: connect the point sources to their samples, fix the
     *        communication pattern and exchange the values of the current solution
     * \note Collective.
     */
    void finishDistributedSetup()
    {
        idBulkSample_.resize(idCounter_);
        idLowDimSample_.resize(idCounter_);
        for (std::size_t id = 0; id < idCounter_; ++id)
        {
            const auto& [bulkKey, lowDimKey] = pointSourceSampleKeys_[id];
            idBulkSample_[id] = sampleIndices_[bulkSample].at(bulkKey);
            idLowDimSample_[id] = sampleIndices_[lowDimSample].at(lowDimKey);
        }

        sumSamples_.assign(distributedSums_.size(), {bulkSample, invalidIndex});
        for (const auto kind : {bulkSample, lowDimSample})
            for (std::size_t s = 0; s < samples_[kind].size(); ++s)
                if (samples_[kind][s].sum != invalidIndex)
                    sumSamples_[samples_[kind][s].sum] = {kind, s};

        distributedSums_.finalize();

        for (auto& indices : sampleIndices_)
            indices.clear();
        importedPointSources_.clear();
        pointSourceSampleKeys_.clear();

        using namespace Dune::Hybrid;
        forEach(integralRange(std::integral_constant<std::size_t, MDTraits::numSubDomains>{}), [&](const auto domainIdx)
        {
            for (auto&& stencil : this->couplingStencils(domainIdx))
            {
                std::sort(stencil.second.begin(), stencil.second.end());
                stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());
            }
        });

        setupRemoteDerivatives_();
        exchangeDistributedValues_();
    }

    /*!
     * \brief Call f(eIdx, id) for every point source id whose values enter the source term of the bulk element eIdx
     * \note Coupling managers whose bulk sources are not point sources of the bulk domain implement their own.
     */
    template<class Function>
    void forEachBulkElementSourceId(Function&& f) const
    {
        for (const auto& source : std::get<bulkIdx>(pointSources_))
            f(source.elementIndex(), source.id());
    }

    //! a bulk point source, sent by the owner of the lower-dimensional element to the process of the bulk element
    struct BulkPointSourceRecord
    {
        DistributedKey key;
        GlobalPosition position;
        Scalar quadratureWeight;
        Scalar integrationElement;
        int embeddings;
        std::size_t bulkElementIdx;
        DistributedKey bulkSampleKey;
        DistributedKey lowDimSampleKey;
        Scalar distance;
    };

    //! Append a bulk point source to a message
    static void writeBulkPointSource(std::vector<char>& message, const BulkPointSourceRecord& record)
    {
        EmbeddedCoupling::Detail::MessageWriter writer(message);
        writer.put(record.key);
        writer.put(record.position);
        writer.put(record.quadratureWeight);
        writer.put(record.integrationElement);
        writer.put(record.embeddings);
        writer.put(record.bulkElementIdx);
        writer.put(record.bulkSampleKey);
        writer.put(record.lowDimSampleKey);
        writer.put(record.distance);
    }

    /*!
     * \brief Create the bulk point sources sent by the owners of the lower-dimensional elements
     * \param messages per rank the bulk point sources written with writeBulkPointSource
     * \param f called as f(id, record) for every bulk point source, with the process-local id of the point source
     * \note Collective.
     */
    template<class Function>
    void receiveBulkPointSources(std::vector<std::vector<char>>& messages, Function&& f)
    {
        using namespace EmbeddedCoupling::Detail;
        const int rank = bulkComm().rank();
        for (const auto& [sender, message] : exchangeMessages(bulkComm(), messages))
        {
            MessageReader reader(message);
            while (!reader.done())
            {
                BulkPointSourceRecord record;
                record.key = reader.get<DistributedKey>();
                record.position = reader.get<GlobalPosition>();
                record.quadratureWeight = reader.get<Scalar>();
                record.integrationElement = reader.get<Scalar>();
                record.embeddings = reader.get<int>();
                record.bulkElementIdx = reader.get<std::size_t>();
                record.bulkSampleKey = reader.get<DistributedKey>();
                record.lowDimSampleKey = reader.get<DistributedKey>();
                record.distance = reader.get<Scalar>();

                const auto id = pointSourceId(record.key);
                pointSources(bulkIdx).emplace_back(record.position, id, record.quadratureWeight,
                                                   record.integrationElement, record.bulkElementIdx);
                pointSources(bulkIdx).back().setEmbeddings(record.embeddings);
                addInterpolation(bulkIdx, pointSourceData_[id], record.bulkElementIdx, record.position);
                setAverageDistance(id, record.distance);
                setPointSourceSamples(id, record.bulkSampleKey, record.lowDimSampleKey);

                if (sender == rank)
                    addToCouplingStencil(bulkIdx, record.bulkElementIdx, lowDimIdx, pointSourceData_[id].lowDimElementIdx());

                f(id, record);
            }
        }
    }

    //! \copydoc receiveBulkPointSources
    void receiveBulkPointSources(std::vector<std::vector<char>>& messages)
    { receiveBulkPointSources(messages, [](std::size_t, const BulkPointSourceRecord&){}); }

    /*!
     * \brief Couple the lower-dimensional elements of this process to the dofs of this process in the bulk
     *        samples of their point sources
     * \note Call after distributeSamples.
     */
    void addOwnedLowDimCouplingStencils()
    {
        for (std::size_t id = 0; id < idCounter_; ++id)
        {
            const auto lowDimElementIdx = pointSourceData_[id].lowDimElementIdx();
            if (lowDimElementIdx == std::numeric_limits<GridIndex<lowDimIdx>>::max())
                continue;

            const auto& dofs = sampleDofs(bulkSample, pointSourceSampleKeys_[id][bulkSample]);
            auto& stencil = this->couplingStencils(lowDimIdx)[lowDimElementIdx];
            stencil.insert(stencil.end(), dofs.begin(), dofs.end());
        }
    }

    //! The dofs of this process in the terms of a sample, during the setup
    const std::vector<std::size_t>& sampleDofs(SampleKind kind, const DistributedKey& key)
    { return samples_[kind][sampleIndex_(kind, key)].dofs; }

    //! The distributed bounding box tree of the bulk elements
    auto makeDistributedBulkTree() const
    {
        using EntitySet = GridViewGeometricEntitySet<GridView<bulkIdx>, 0, ElementMapper<bulkIdx>>;
        const auto& gridGeometry = gridDiscretization(this->problem(bulkIdx));
        return std::make_shared<DistributedBoundingBoxTree<EntitySet>>(
            std::make_shared<EntitySet>(gridGeometry.gridView(), gridGeometry.elementMapper())
        );
    }

    /*!
     * \brief Sort the extended source stencils and remove duplicates and the dofs of the element itself
     */
    template<class Stencils>
    void makeExtendedSourceStencilsUnique(Stencils& stencils) const
    {
        for (auto&& stencil : stencils)
        {
            std::sort(stencil.second.begin(), stencil.second.end());
            stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());

            if constexpr (isCVFE<bulkIdx>())
            {
                const auto& indices = this->elementDofIndices(bulkIdx, stencil.first);
                stencil.second.erase(std::remove_if(stencil.second.begin(), stencil.second.end(),
                                                    [&](auto i){ return std::find(indices.begin(), indices.end(), i) != indices.end(); }),
                                     stencil.second.end());
            }
            else
            {
                stencil.second.erase(std::remove_if(stencil.second.begin(), stencil.second.end(),
                                                    [&](auto i){ return i == stencil.first; }),
                                     stencil.second.end());
            }
        }
    }

    //! Print the number of point sources and the setup time
    void printDistributedSetup(const Dune::Timer& watch) const
    {
        const auto& comm = bulkComm();
        const auto numSources = comm.sum(pointSources(lowDimIdx).size());
        if (comm.rank() == 0)
            std::cout << "Initialized " << numSources << " point sources on " << comm.size()
                      << " processes in " << watch.elapsed() << " seconds." << std::endl;
    }

    /*!
     * \brief Add values of the lower-dimensional elements that the coupling terms need
     * \param numValues the number of values per element
     * \param f the function f(lowDimElementIdx, values) writing the values into the std::span<Scalar> values
     * \return the index of the first of these values, see lowDimElementValue
     */
    template<class Function>
    std::size_t addLowDimElementValues(std::size_t numValues, Function&& f)
    {
        if (distributed_)
            DUNE_THROW(Dune::InvalidStateException, "Values of the lower-dimensional elements have to be added before the point sources are computed");

        const auto first = numLowDimElementValues_;
        lowDimElementValueFunctions_.emplace_back(first, numValues, std::forward<Function>(f));
        numLowDimElementValues_ += numValues;
        return first;
    }

    /*!
     * \brief The value i of the lower-dimensional element of point source id, see addLowDimElementValues
     * \note On distributed grids the values are those of the last solution update on every process.
     */
    Scalar lowDimElementValue(std::size_t id, std::size_t i) const
    {
        if (distributed_)
        {
            const auto s = idLowDimSample_[id];
            const auto& sample = samples_[lowDimSample][s];
            if (sample.lowDimElement == invalidIndex)
                return distributedSums_.received(sample.sum)[PrimaryVariables<lowDimIdx>::size() + i];
            return lowDimElementValues_[s*numLowDimElementValues_ + i];
        }

        for (const auto& [first, numValues, f] : lowDimElementValueFunctions_)
        {
            if (i < first || i >= first + numValues)
                continue;
            std::vector<Scalar> values(numValues);
            f(pointSourceData_[id].lowDimElementIdx(), std::span<Scalar>(values));
            return values[i - first];
        }

        DUNE_THROW(Dune::RangeError, "No value " << i << " of the lower-dimensional elements");
    }

    // \}

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    //! id generator for point sources
    std::size_t idCounter_ = 0;

private:
    //! the terms of a sample on this process and its sum over the processes
    struct Sample
    {
        std::vector<std::size_t> dofs;
        std::vector<Scalar> weights;
        std::size_t sum = invalidIndex;
        std::size_t lowDimElement = invalidIndex;
    };

    //! a column of the Jacobian entries of a sum: the global index of a degree of freedom of another process and its weight
    struct RemoteColumn
    {
        std::int64_t index;
        Scalar weight;
    };

    //! a received sum of a coupling manager deflected for numeric differentiation on this thread
    struct SumDeflection
    {
        const void* manager = nullptr;
        std::size_t sum = invalidIndex;
        std::size_t component = 0;
        Scalar delta = 0.0;
    };

    //! The number of rows of the element residual of an element of domain i
    template<std::size_t i>
    std::size_t numElementResidualRows_(Dune::index_constant<i> domainI, std::size_t eIdx) const
    {
        if constexpr (isCVFE<i>())
            return this->elementDofIndices(domainI, eIdx).size();
        else
            return 1;
    }

    //! For every element, the received sums of the samples of its point sources, and storage for the derivatives
    void setupRemoteDerivatives_()
    {
        std::array<std::vector<std::vector<std::size_t>>, 2> elementSums;
        elementSums[bulkIdx].resize(gridView(bulkIdx).size(0));
        elementSums[lowDimIdx].resize(gridView(lowDimIdx).size(0));
        const auto addSums = [&](std::size_t domain, std::size_t eIdx, std::size_t id)
        {
            for (const auto* sample : {&samples_[bulkSample][idBulkSample_[id]], &samples_[lowDimSample][idLowDimSample_[id]]})
                if (sample->sum != invalidIndex && distributedSums_.receives(sample->sum))
                    elementSums[domain][eIdx].push_back(sample->sum);
        };
        asImp_().forEachBulkElementSourceId([&](std::size_t eIdx, std::size_t id){ addSums(bulkIdx, eIdx, id); });
        for (const auto& source : pointSources(lowDimIdx))
            addSums(lowDimIdx, source.elementIndex(), source.id());

        using namespace Dune::Hybrid;
        forEach(integralRange(std::integral_constant<std::size_t, MDTraits::numSubDomains>{}), [&](const auto domainIdx)
        {
            using GG = std::decay_t<decltype(gridDiscretization(this->problem(domainIdx)))>;
            if constexpr (Dumux::Detail::hasPeriodicDofMap<GG>())
                if (!gridDiscretization(this->problem(domainIdx)).periodicDofMap().empty())
                    DUNE_THROW(Dune::NotImplemented, "Embedded coupling on distributed grids with periodic boundaries");

            constexpr auto numEq = PrimaryVariables<domainIdx>::size();
            auto& offsets = remoteSumOffsets_[domainIdx];
            auto& sums = remoteSums_[domainIdx];
            auto& derivativeOffsets = remoteDerivativeOffsets_[domainIdx];
            offsets.assign(1, 0);
            sums.clear();
            derivativeOffsets.clear();
            std::size_t numDerivatives = 0;
            for (std::size_t eIdx = 0; eIdx < elementSums[domainIdx].size(); ++eIdx)
            {
                auto& elementSumsOfElement = elementSums[domainIdx][eIdx];
                std::sort(elementSumsOfElement.begin(), elementSumsOfElement.end());
                elementSumsOfElement.erase(std::unique(elementSumsOfElement.begin(), elementSumsOfElement.end()), elementSumsOfElement.end());
                for (const auto sum : elementSumsOfElement)
                {
                    sums.push_back(sum);
                    derivativeOffsets.push_back(numDerivatives);
                    const auto numComponents = sumSamples_[sum].first == bulkSample ? PrimaryVariables<bulkIdx>::size()
                                                                                     : PrimaryVariables<lowDimIdx>::size();
                    numDerivatives += numElementResidualRows_(domainIdx, eIdx)*numEq*numComponents;
                }
                offsets.push_back(sums.size());
            }
            remoteDerivatives_[domainIdx].assign(numDerivatives, 0.0);
        });
        remoteColumns_.clear();
    }

    /*!
     * \brief Set the derivatives with respect to received sums of the equations of the bound element to zero that
     *        the local assembler replaces by Dirichlet conditions or internal Dirichlet constraints
     */
    template<std::size_t i, class LocalAssemblerI>
    void zeroConstrainedDerivatives_(Dune::index_constant<i> domainI, const LocalAssemblerI& localAssemblerI)
    {
        constexpr auto numEq = PrimaryVariables<i>::size();
        const auto& elemDisc = [&]() -> decltype(auto) {
            if constexpr (requires { localAssemblerI.elemDisc(); })
                return localAssemblerI.elemDisc();
            else
                return localAssemblerI.fvGeometry();
        }();
        using ElementDiscretization = std::decay_t<decltype(elemDisc)>;
        [[maybe_unused]] const auto& problem = localAssemblerI.problem();
        const auto& element = localAssemblerI.element();
        const auto eIdx = gridDiscretization(this->problem(domainI)).elementMapper().index(element);

        [[maybe_unused]] const auto zeroEquation = [&](std::size_t localRow, std::size_t eqIdx)
        {
            for (std::size_t m = remoteSumOffsets_[i][eIdx]; m < remoteSumOffsets_[i][eIdx+1]; ++m)
            {
                const auto numComponents = sumSamples_[remoteSums_[i][m]].first == bulkSample ? PrimaryVariables<bulkIdx>::size()
                                                                                                : PrimaryVariables<lowDimIdx>::size();
                Scalar* derivatives = remoteDerivatives_[i].data() + remoteDerivativeOffsets_[i][m];
                std::fill_n(derivatives + (localRow*numEq + eqIdx)*numComponents, numComponents, 0.0);
            }
        };

        // element discretizations without control volumes only have constraints of the problem
        if constexpr (requires { scvs(elemDisc); })
        {
            for (const auto& scv : scvs(elemDisc))
            {
                const std::size_t localRow = [&]() -> std::size_t {
                    if constexpr (isCVFE<i>())
                        return scv.localDofIndex();
                    else
                        return 0;
                }();
                if constexpr (isCVFE<i>() && !Dumux::Detail::hasProblemBoundaryTypesForFaceFunction<std::decay_t<decltype(problem)>, ElementDiscretization>()
                              && requires { localAssemblerI.elemBcTypes(); })
                {
                    if (localAssemblerI.elemBcTypes().hasDirichlet())
                    {
                        const auto bcTypes = localAssemblerI.elemBcTypes().get(elemDisc, scv);
                        for (std::size_t eqIdx = 0; eqIdx < numEq; ++eqIdx)
                            if (bcTypes.isDirichlet(eqIdx) && !bcTypes.isCouplingDirichlet(eqIdx))
                                zeroEquation(localRow, eqIdx);
                    }
                }

                using ProblemI = std::decay_t<decltype(problem)>;
                if constexpr (requires { ProblemI::enableInternalDirichletConstraints(); })
                {
                    if constexpr (ProblemI::enableInternalDirichletConstraints())
                    {
                        const auto constraints = problem.hasInternalDirichletConstraint(element, scv);
                        for (std::size_t eqIdx = 0; eqIdx < numEq; ++eqIdx)
                            if (constraints[eqIdx])
                                zeroEquation(localRow, eqIdx);
                    }
                }
            }
        }
    }

    //! Call f(domainI, dofIdx, eqIdx, sum, pvIdx, derivative) for every derivative with respect to a received sum
    template<class Function>
    void forEachRemoteDerivative_(Function&& f) const
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(std::integral_constant<std::size_t, MDTraits::numSubDomains>{}), [&](const auto domainIdx)
        {
            constexpr auto numEq = PrimaryVariables<domainIdx>::size();

            // the equations the assembler replaces by the constraints of the problem
            std::set<std::pair<std::size_t, std::size_t>> constrained;
            if constexpr (Dumux::Detail::hasSubProblemGlobalConstraints<Problem<domainIdx>>())
                for (const auto& constraintData : this->problem(domainIdx).constraints())
                    for (std::size_t eqIdx = 0; eqIdx < numEq; ++eqIdx)
                        if (constraintData.constraintInfo().isConstraintEquation(eqIdx))
                            constrained.emplace(constraintData.dofIndex(), eqIdx);

            const auto& offsets = remoteSumOffsets_[domainIdx];
            for (std::size_t eIdx = 0; eIdx + 1 < offsets.size(); ++eIdx)
            {
                for (std::size_t m = offsets[eIdx]; m < offsets[eIdx+1]; ++m)
                {
                    const auto sum = remoteSums_[domainIdx][m];
                    const auto numComponents = sumSamples_[sum].first == bulkSample ? PrimaryVariables<bulkIdx>::size()
                                                                                     : PrimaryVariables<lowDimIdx>::size();
                    const Scalar* derivatives = remoteDerivatives_[domainIdx].data() + remoteDerivativeOffsets_[domainIdx][m];
                    for (std::size_t k = 0; k < numElementResidualRows_(domainIdx, eIdx); ++k)
                    {
                        const std::size_t dofIdx = [&]{
                            if constexpr (isCVFE<domainIdx>())
                                return static_cast<std::size_t>(this->elementDofIndices(domainIdx, eIdx)[k]);
                            else
                                return eIdx;
                        }();
                        for (std::size_t eqIdx = 0; eqIdx < numEq; ++eqIdx)
                            if (!constrained.contains({dofIdx, eqIdx}))
                                for (std::size_t c = 0; c < numComponents; ++c)
                                    f(domainIdx, dofIdx, eqIdx, sum, c, derivatives[(k*numEq + eqIdx)*numComponents + c]);
                    }
                }
            }
        });
    }

    //! the line mode on distributed grids: a point source per quadrature point of a piece and bulk element
    void computeDistributedPointSourceData_(std::size_t order)
    {
        Dune::Timer watch;
        const auto pieces = beginDistributedSetup();
        const auto& comm = bulkComm();
        const int rank = comm.rank();

        std::vector<SampleSpec> specs;
        std::vector<std::vector<char>> sourceMessages(comm.size());
        for (const auto& piece : pieces)
        {
            const auto geometry = piece.geometry();
            const auto& quad = Dune::QuadratureRules<Scalar, lowDimDim>::rule(geometry.type(), order);
            const int embeddings = piece.bulkElements.size();
            for (const auto& qp : quad)
            {
                const auto globalPos = geometry.global(qp.position());
                const auto qpweight = qp.weight();
                const auto ie = geometry.integrationElement(qp.position());

                const auto lowDimKey = newSampleKey(lowDimSample);
                const auto lowDimSpec = specs.size();
                specs.push_back(SampleSpec{lowDimSample, lowDimKey, {rank}, {rank}, {{rank, piece.lowDimElementIdx, globalPos, 1.0}}});

                for (int n = 0; n < embeddings; ++n)
                {
                    const auto [bulkRank, bulkElementIdx] = piece.bulkElements[n];
                    const auto id = addOwnedPointSource(piece.lowDimElementIdx, globalPos);
                    pointSources(lowDimIdx).emplace_back(globalPos, id, qpweight, ie, piece.lowDimElementIdx);
                    pointSources(lowDimIdx).back().setEmbeddings(embeddings);

                    const auto distance = averageDistancePointGeometry(globalPos, piece.bulkGeometries[n]);
                    setAverageDistance(id, distance);

                    const auto bulkKey = newSampleKey(bulkSample);
                    specs.push_back(SampleSpec{bulkSample, bulkKey, {bulkRank}, {rank, bulkRank}, {{bulkRank, bulkElementIdx, globalPos, 1.0}}});
                    specs[lowDimSpec].consumers.push_back(bulkRank);
                    setPointSourceSamples(id, bulkKey, lowDimKey);

                    writeBulkPointSource(sourceMessages[bulkRank], BulkPointSourceRecord{
                        ownedPointSourceKey(id), globalPos, qpweight, ie, embeddings, bulkElementIdx, bulkKey, lowDimKey, distance
                    });
                }
            }
        }

        distributeSamples(specs);
        receiveBulkPointSources(sourceMessages);
        addOwnedLowDimCouplingStencils();
        finishDistributedSetup();
        printDistributedSetup(watch);
    }

    //! A new point source id with default data
    std::size_t addPointSourceId_()
    {
        pointSourceData_.emplace_back();
        averageDistanceToBulkCell_.push_back(0.0);
        pointSourceSampleKeys_.emplace_back();
        return idCounter_++;
    }

    //! The local index of the sample with the given key, added if this process does not have it yet
    std::size_t sampleIndex_(SampleKind kind, const DistributedKey& key)
    {
        const auto [it, inserted] = sampleIndices_[kind].try_emplace(key, samples_[kind].size());
        if (inserted)
            samples_[kind].emplace_back();
        return it->second;
    }

    //! Add the value at a point of an element of this process times weight to a sample
    void addSampleTerm_(SampleKind kind, Sample& sample, std::size_t eIdx, const GlobalPosition& point, Scalar weight)
    {
        const auto addWeight = [&](std::size_t dof, Scalar w)
        {
            const auto it = std::find(sample.dofs.begin(), sample.dofs.end(), dof);
            if (it != sample.dofs.end())
                sample.weights[it - sample.dofs.begin()] += w;
            else
            {
                sample.dofs.push_back(dof);
                sample.weights.push_back(w);
            }
        };

        const auto addTerm = [&](auto domainI)
        {
            if constexpr (isCVFE<decltype(domainI)::value>())
            {
                const auto& gg = gridDiscretization(this->problem(domainI));
                std::vector<Dune::FieldVector<Scalar, 1>> shapeValues;
                this->getShapeValues(domainI, gg, gg.element(eIdx).geometry(), point, shapeValues);
                const auto& dofs = this->elementDofIndices(domainI, eIdx);
                for (std::size_t j = 0; j < dofs.size(); ++j)
                    addWeight(dofs[j], weight*shapeValues[j][0]);
            }
            else
                addWeight(eIdx, weight);
        };

        if (kind == bulkSample)
            addTerm(bulkIdx);
        else
        {
            addTerm(lowDimIdx);
            sample.lowDimElement = eIdx;
        }
    }

    //! The number of values of a sample: the primary variables, for lower-dimensional samples also the element values
    std::size_t numSampleComponents_(SampleKind kind) const
    {
        return kind == bulkSample ? PrimaryVariables<bulkIdx>::size()
                                  : PrimaryVariables<lowDimIdx>::size() + numLowDimElementValues_;
    }

    //! Add the terms of this process of a sample for the solution sol to values
    template<class Sol, class Values>
    static void addLocalTerms_(const Sample& sample, const Sol& sol, Values& values, std::size_t numComponents)
    {
        for (std::size_t j = 0; j < sample.dofs.size(); ++j)
            for (std::size_t c = 0; c < numComponents; ++c)
                values[c] += sol[sample.dofs[j]][c]*sample.weights[j];
    }

    //! The value of a sample: the terms of this process for the current solution and the received terms
    template<class Values, class Sol>
    Values sampleValue_(const Sample& sample, const Sol& sol) const
    {
        Values values(0.0);
        addLocalTerms_(sample, sol, values, Values::size());
        if (sample.sum != invalidIndex)
        {
            const auto received = distributedSums_.received(sample.sum);
            for (std::size_t c = 0; c < Values::size(); ++c)
                values[c] += received[c];
            if (sumDeflection_.manager == this && sumDeflection_.sum == sample.sum)
                values[sumDeflection_.component] += sumDeflection_.delta;
        }
        return values;
    }

    //! Evaluate the element values of the lower-dimensional samples of this process and send the terms of
    //! this process of the samples to the processes that need them
    void exchangeDistributedValues_()
    {
        const auto& lowDimSamples = samples_[lowDimSample];
        const auto n = numLowDimElementValues_;
        lowDimElementValues_.resize(lowDimSamples.size()*n);
        for (std::size_t s = 0; s < lowDimSamples.size(); ++s)
            if (lowDimSamples[s].lowDimElement != invalidIndex)
                for (const auto& [first, numValues, f] : lowDimElementValueFunctions_)
                    f(lowDimSamples[s].lowDimElement, std::span<Scalar>(lowDimElementValues_.data() + s*n + first, numValues));

        static constexpr std::size_t numLowDimPriVars = PrimaryVariables<lowDimIdx>::size();
        distributedSums_.exchange([&](std::size_t sum, std::span<Scalar> values)
        {
            const auto [kind, s] = sumSamples_[sum];
            const auto& sample = samples_[kind][s];
            if (kind == bulkSample)
                addLocalTerms_(sample, this->curSol(bulkIdx), values, PrimaryVariables<bulkIdx>::size());
            else
            {
                addLocalTerms_(sample, this->curSol(lowDimIdx), values, numLowDimPriVars);
                if (sample.lowDimElement != invalidIndex)
                    std::copy_n(lowDimElementValues_.begin() + s*n, n, values.begin() + numLowDimPriVars);
            }
        });
    }

    //! Whether two intersections are identical, with the tolerance of the sequential glue
    static bool cornersMatch_(const std::vector<GlobalPosition>& corners, const std::vector<GlobalPosition>& otherCorners)
    {
        if (otherCorners.size() != corners.size())
            return false;

        using std::max;
        Scalar eps2 = std::numeric_limits<Scalar>::min();
        for (std::size_t i = 1; i < corners.size(); ++i)
            eps2 = max(eps2, (corners[i] - corners[0]).two_norm2());
        eps2 *= 1.5e-7*1.5e-7;

        return std::all_of(corners.begin(), corners.end(), [&](const auto& corner)
        {
            return std::any_of(otherCorners.begin(), otherCorners.end(),
                               [&](const auto& other){ return (corner - other).two_norm2() < eps2; });
        });
    }

    //! the point source in both domains
    std::tuple<std::vector<PointSource<bulkIdx>>, std::vector<PointSource<lowDimIdx>>> pointSources_;
    std::tuple<std::vector<PointSourceContext<bulkIdx>>, std::vector<PointSourceContext<lowDimIdx>>> additionalPointSourceContexts_;
    std::vector<PointSourceData> pointSourceData_;
    std::vector<Scalar> averageDistanceToBulkCell_;

    //! Stencil data
    std::tuple<std::vector<std::vector<GridIndex<bulkIdx>>>,
               std::vector<std::vector<GridIndex<lowDimIdx>>>> elementDofIndices_;
    std::tuple<CouplingStencils<bulkIdx>, CouplingStencils<lowDimIdx>> couplingStencils_;
    std::tuple<CouplingStencil<bulkIdx>, CouplingStencil<lowDimIdx>> emptyStencil_;

    //! The glue object
    std::shared_ptr<GlueType> glue_;

    //! integration order for coupling source
    int integrationOrder_ = 1;

    //! coloring for multithreaded assembly
    std::tuple<
        std::deque<std::vector<ElementSeed<bulkIdx>>>,
        std::deque<std::vector<ElementSeed<lowDimIdx>>>
    > elementSets_;

    //! point sources on distributed grids
    bool distributed_ = false;
    DistributedSums<Scalar> distributedSums_;
    std::array<std::vector<Sample>, 2> samples_;
    std::array<std::map<DistributedKey, std::size_t>, 2> sampleIndices_;
    std::vector<std::pair<SampleKind, std::size_t>> sumSamples_;
    std::vector<std::size_t> idBulkSample_;
    std::vector<std::size_t> idLowDimSample_;

    //! per domain and element (CSR) the received sums its source term depends on, and the derivatives with respect to them
    std::array<std::vector<std::size_t>, 2> remoteSumOffsets_;
    std::array<std::vector<std::size_t>, 2> remoteSums_;
    std::array<std::vector<std::size_t>, 2> remoteDerivativeOffsets_;
    std::array<std::vector<Scalar>, 2> remoteDerivatives_;
    std::vector<std::vector<RemoteColumn>> remoteColumns_;
    std::vector<Scalar> remoteProductValues_;
    static inline thread_local SumDeflection sumDeflection_{};
    std::vector<std::array<DistributedKey, 2>> pointSourceSampleKeys_;
    std::map<DistributedKey, std::size_t> importedPointSources_;
    std::size_t numOwnedSamples_ = 0;

    //! values of the lower-dimensional elements the coupling terms need: first value, number of values, function
    std::vector<std::tuple<std::size_t, std::size_t, std::function<void(std::size_t, std::span<Scalar>)>>> lowDimElementValueFunctions_;
    std::size_t numLowDimElementValues_ = 0;
    //! the values of the lower-dimensional elements of the samples of this process at the last solution update
    std::vector<Scalar> lowDimElementValues_;
    std::size_t userElementValues_ = invalidIndex;
};

} // end namespace Dumux

#endif
