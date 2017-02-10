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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Facet
 * \brief Creates mappings of the entities from the bulk and the facet subdomain
 */

#ifndef DUMUX_MIXEDDIMENSION_CCMPFA_FACET_COUPLINGMAPPER_HH
#define DUMUX_MIXEDDIMENSION_CCMPFA_FACET_COUPLINGMAPPER_HH

#include <dune/common/version.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/gridfactory.hh>

#include <dumux/common/math.hh>
#include <dumux/common/propertysystem.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/valgrind.hh>

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(BulkProblemTypeTag);
NEW_PROP_TAG(LowDimProblemTypeTag);
}

/*!
 * \ingroup FacetCoupling
 * \brief Creates mappings of the entities from the bulk and the facet subdomain
 */
template<class TypeTag>
class CCMpfaFacetCouplingMapper
{
    using GridCreator = typename GET_PROP_TYPE(TypeTag, GridCreator);

    using BulkProblemTypeTag = typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag);
    using LowDimProblemTypeTag = typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag);

    using BulkProblem = typename GET_PROP_TYPE(BulkProblemTypeTag, Problem);
    using LowDimProblem = typename GET_PROP_TYPE(LowDimProblemTypeTag, Problem);

    using BulkGridView = typename GET_PROP_TYPE(BulkProblemTypeTag, GridView);
    using LowDimGridView = typename GET_PROP_TYPE(LowDimProblemTypeTag, GridView);

    using BulkIndexType = typename BulkGridView::IndexSet::IndexType;
    using LowDimIndexType = typename LowDimGridView::IndexSet::IndexType;

    using BulkElement = typename BulkGridView::template Codim<0>::Entity;
    using LowDimElement = typename LowDimGridView::template Codim<0>::Entity;

    using BulkSubControlVolumeFace = typename GET_PROP_TYPE(BulkProblemTypeTag, SubControlVolumeFace);
    using BulkMpfaHelper = typename GET_PROP_TYPE(BulkProblemTypeTag, MpfaHelper);

    static constexpr int lowDimDim = LowDimGridView::dimension;

public:

    struct BulkDataInLowDim
    {
        // indicates whether or not this element is coupled
        bool isCoupled;

        // The list of bulk element indices connected to a low dim element
        std::vector<BulkIndexType> couplingStencil;

        // We store for each directly connected element a list of
        // global scvf indices over which the two domains are coupled
        std::vector<std::pair<BulkIndexType, std::vector<BulkIndexType>>> elementScvfList;

        // The constructor
        BulkDataInLowDim() : isCoupled(false) {}
    };

    struct LowDimDataInBulk
    {
        // indicates whether or not this element is coupled
        bool isCoupled;

        // The list of low dim elements with influence on this bulk element
        std::vector<LowDimIndexType> couplingStencil;

        // to each scvf of the bulk element we store
        // the low dim element index that it is connected to
        std::vector<std::pair<BulkIndexType, LowDimIndexType>> scvfToLowDimData;

        // The constructor
        LowDimDataInBulk() : isCoupled(false) {}

        // For a given scvf of the element this method returns a pair of a bool and an index. The boolean
        // indicates if this scvf does couple to a low dim element, the index is the low dim element index
        std::pair<bool, LowDimIndexType> getScvfCouplingData(const BulkSubControlVolumeFace& scvf) const
        { return getScvfCouplingData(scvf.index()); }

        std::pair<bool, LowDimIndexType> getScvfCouplingData(BulkIndexType scvfIdx) const
        {
            for (const auto& pair : scvfToLowDimData)
                if (pair.first == scvfIdx)
                    return std::make_pair(true, pair.second);
            return std::make_pair(false, 0);
        }
    };

    //! The constructor
    CCMpfaFacetCouplingMapper(BulkProblem &bulkProblem, LowDimProblem &lowDimProblem)
    : isInitialized_(false),
      bulkProblemPtr_(&bulkProblem),
      lowDimProblemPtr_(&lowDimProblem)
    {}

    /*!
     * \brief initializes the maps
     *
     * \param bulkProblem The bulk sub problem
     * \param lowDimProblem The lower-dimensional sub problem
     */
    void init()
    {
        const auto& bulkGridView = bulkProblem_().gridView();
        const auto& lowDimGridView = lowDimProblem_().gridView();

        lowDimCouplingData_.resize(lowDimGridView.size(0));
        bulkCouplingData_.resize(bulkGridView.size(0));

        for (const auto& lowDimElement : elements(lowDimGridView))
        {
            const auto lowDimElementGlobalIdx = lowDimProblem_().elementMapper().index(lowDimElement);

            // The indices of bulk elements in which the actual low dim element is a facet
            const auto& bulkElementIndices = GridCreator::getCoupledBulkElementIndices(lowDimElement);
            if (bulkElementIndices.size() > 1)
            {
                lowDimCouplingData_[lowDimElementGlobalIdx].isCoupled = true;

                for (auto bulkElementIdx : bulkElementIndices)
                {
                    bulkCouplingData_[bulkElementIdx].isCoupled = true;

                    const auto bulkElement = bulkProblem_().model().globalFvGeometry().element(bulkElementIdx);
                    auto bulkFvGeometry = localView(bulkProblem_().model().globalFvGeometry());
                    bulkFvGeometry.bindElement(bulkElement);

                    // find the scvfs in the bulk element that "touch" the facet element
                    std::vector<BulkIndexType> scvfIndices;
                    for (const auto& scvf : scvfs(bulkFvGeometry))
                    {
                        if (scvf.boundary())
                            continue;

                        const auto anyOtherIdx = bulkElementIndices[0] == bulkElementIdx ?
                                                 bulkElementIndices[1] :
                                                 bulkElementIndices[0];

                        // check if any of the other bulk indices is in the outside indices of the scvf
                        if (BulkMpfaHelper::contains(scvf.outsideScvIndices(), anyOtherIdx))
                        {
                            scvfIndices.push_back(scvf.index());

                            // add data to the bulk data map
                            bulkCouplingData_[bulkElementIdx].scvfToLowDimData.emplace_back(std::make_pair(scvf.index(), lowDimElementGlobalIdx));

                            // also, insert scvf stencil to the coupling stencil of the low dim element
                            const auto scvfStencil = [&] ()
                            {
                                if (bulkProblem_().model().globalFvGeometry().isInBoundaryInteractionVolume(scvf))
                                    return bulkProblem_().model().globalFvGeometry().boundaryInteractionVolumeSeed(scvf).globalScvIndices();
                                else
                                    return bulkProblem_().model().globalFvGeometry().interactionVolumeSeed(scvf).globalScvIndices();
                            } ();

                            auto& lowDimCouplingStencil = lowDimCouplingData_[lowDimElementGlobalIdx].couplingStencil;
                            lowDimCouplingStencil.insert(lowDimCouplingStencil.end(), scvfStencil.begin(), scvfStencil.end());

                            // also, all the bulk elements in the scvf stencil will be coupled to the actual low dim element
                            for (auto bulkIdx : scvfStencil)
                            {
                                bulkCouplingData_[bulkIdx].isCoupled = true;
                                bulkCouplingData_[bulkIdx].couplingStencil.push_back(lowDimElementGlobalIdx);
                            }
                        }
                    }

                    // insert data into the low dim map
                    lowDimCouplingData_[lowDimElementGlobalIdx].elementScvfList.emplace_back(std::make_pair(bulkElementIdx, std::move(scvfIndices)));
                }
            }
            else if (bulkElementIndices.size() == 1)
                DUNE_THROW(Dune::NotImplemented, "Coupled facet elements on the bulk boundary");

            // make the coupling stencil of the low dim element unique
            auto& stencil = lowDimCouplingData_[lowDimElementGlobalIdx].couplingStencil;
            std::sort(stencil.begin(), stencil.end());
            stencil.erase(std::unique(stencil.begin(), stencil.end()), stencil.end());
        }

        // make the coupling stencils in the bulk map unique
        for (auto& data : bulkCouplingData_)
        {
            auto& stencil = data.couplingStencil;
            std::sort(stencil.begin(), stencil.end());
            stencil.erase(std::unique(stencil.begin(), stencil.end()), stencil.end());
        }

        // state the initialization status
        isInitialized_ = true;
    }

    const LowDimDataInBulk& getBulkCouplingData(const BulkElement& bulkElement) const
    { return getBulkCouplingData(bulkProblem_().elementMapper().index(bulkElement)); }

    const LowDimDataInBulk& getBulkCouplingData(BulkIndexType bulkElementIdx) const
    { return bulkCouplingData_[bulkElementIdx]; }

    const BulkDataInLowDim& getLowDimCouplingData(const LowDimElement& lowDimElement) const
    { return getLowDimCouplingData(lowDimProblem_().elementMapper().index(lowDimElement)); }

    const BulkDataInLowDim& getLowDimCouplingData(LowDimIndexType lowDimElementIdx) const
    { return lowDimCouplingData_[lowDimElementIdx]; }

    bool isInitialized() const
    { return isInitialized_; }

private:
    const BulkProblem& bulkProblem_() const
    { return *bulkProblemPtr_; }

    const LowDimProblem& lowDimProblem_() const
    { return *lowDimProblemPtr_; }

    bool isInitialized_;
    const BulkProblem* bulkProblemPtr_;
    const LowDimProblem* lowDimProblemPtr_;

    std::vector<LowDimDataInBulk> bulkCouplingData_;
    std::vector<BulkDataInLowDim> lowDimCouplingData_;
};

} // end namespace

#endif // DUMUX_FACETCOUPLINGMAPPER_HH
