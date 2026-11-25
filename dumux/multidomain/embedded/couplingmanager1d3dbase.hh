// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedCoupling
 * \brief Coupling manager for low-dimensional domains embedded in the bulk domain.
 */

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGER_1D3D_BASE_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGER_1D3D_BASE_HH

#include <limits>
#include <span>
#include <vector>

#include <dumux/common/tag.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/griddiscretization.hh>
#include <dumux/common/indextraits.hh>

#include <dumux/multidomain/embedded/couplingmanagerbase.hh>

namespace Dumux {

/*!
 * \ingroup EmbeddedCoupling
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *        Point sources on each integration point are computed by an AABB tree.
 */

template<class MDTraits, class Implementation, class PSTraits = DefaultPointSourceTraits<MDTraits>>
class Embedded1d3dCouplingManagerBase
: public EmbeddedCouplingManagerBase<MDTraits, Implementation, PSTraits>
{
    using ThisType = Embedded1d3dCouplingManagerBase<MDTraits, Implementation, PSTraits>;
    using ParentType = EmbeddedCouplingManagerBase<MDTraits, Implementation, PSTraits>;

    using Scalar = typename MDTraits::Scalar;
    using SolutionVector = typename MDTraits::SolutionVector;

    static constexpr auto bulkIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto lowDimIdx = typename MDTraits::template SubDomain<1>::Index();

    // the sub domain type aliases
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;
    template<std::size_t id> using Problem = GetPropType<SubDomainTypeTag<id>, Properties::Problem>;
    template<std::size_t id> using GridGeometry = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using GridIndex = typename IndexTraits<GridView<id>>::GridIndex;

public:

    using ParentType::ParentType;

    void init(std::shared_ptr<Problem<bulkIdx>> bulkProblem,
              std::shared_ptr<Problem<lowDimIdx>> lowDimProblem,
              const SolutionVector& curSol)
    {
        if (radiusValue_ == invalidValue_)
            radiusValue_ = this->addLowDimElementValues(1, [this](std::size_t eIdx, std::span<Scalar> values)
            { values[0] = this->problem(lowDimIdx).spatialParams().radius(eIdx); });

        ParentType::init(bulkProblem, lowDimProblem, curSol);
        computeLowDimVolumeFractions_();
    }

    //! Compute the low dim volume fraction in the bulk domain cells
    [[deprecated("Will be removed after release 3.11. Only needed for the deprecated lowDimVolume and lowDimVolumeFraction.")]]
    void computeLowDimVolumeFractions()
    { computeLowDimVolumeFractions_(); }

    /*!
     * \brief Methods to be accessed by the subproblems
     */
    // \{

    //! Returns the radius of the lower dimensional element
    Scalar radius(std::size_t id) const
    {
        if (this->isDistributed())
            return this->lowDimElementValue(id, radiusValue_);

        const auto& data = this->pointSourceData()[id];
        return this->problem(lowDimIdx).spatialParams().radius(data.lowDimElementIdx());
    }

    //! The volume the lower dimensional domain occupies in the bulk domain element
    // For one-dimensional low dim domain we assume radial tubes
    [[deprecated("Will be removed after release 3.11. This is not a general implementation: the volume of "
                 "a cylinder around a segment is counted in full in every bulk element containing the segment "
                 "(several for a segment on a face or an edge) and not in the elements the cylinder reaches "
                 "beyond the segment. Use EmbeddedCoupling::lowDimVolumePerBulkElement "
                 "from dumux/multidomain/embedded/lowdimvolume.hh.")]]
    Scalar lowDimVolume(const Element<bulkIdx>& element) const
    {
        if (this->isDistributed())
            DUNE_THROW(Dune::NotImplemented, "lowDimVolume on distributed grids");
        const auto eIdx = gridDiscretization(this->problem(bulkIdx)).elementMapper().index(element);
        return lowDimVolumeInBulkElement_[eIdx];
    }

    //! The volume fraction the lower dimensional domain occupies in the bulk domain element
    // For one-dimensional low dim domain we assume radial tubes
    [[deprecated("Will be removed after release 3.11. This is not a general implementation: see lowDimVolume. "
                 "Divide the result of EmbeddedCoupling::lowDimVolumePerBulkElement "
                 "from dumux/multidomain/embedded/lowdimvolume.hh by the element volume.")]]
    Scalar lowDimVolumeFraction(const Element<bulkIdx>& element) const
    {
        if (this->isDistributed())
            DUNE_THROW(Dune::NotImplemented, "lowDimVolumeFraction on distributed grids");
        const auto eIdx = gridDiscretization(this->problem(bulkIdx)).elementMapper().index(element);
        return lowDimVolumeInBulkElement_[eIdx] / element.geometry().volume();
    }

    // \}

private:
    void computeLowDimVolumeFractions_()
    {
        // resize the storage vector
        lowDimVolumeInBulkElement_.resize(this->gridView(bulkIdx).size(0));
        // get references to the grid geometries
        const auto& lowDimGridGeometry = gridDiscretization(this->problem(lowDimIdx));
        const auto& bulkGridGeometry = gridDiscretization(this->problem(bulkIdx));

        // compute the low dim volume fractions
        for (const auto& is : intersections(this->glue()))
        {
            // all inside elements are identical...
            const auto& inside = is.targetEntity(0);
            const auto intersectionGeometry = is.geometry();
            const auto lowDimElementIdx = lowDimGridGeometry.elementMapper().index(inside);

            // compute the volume the low-dim domain occupies in the bulk domain if it were full-dimensional
            const auto radius = this->problem(lowDimIdx).spatialParams().radius(lowDimElementIdx);
            for (int outsideIdx = 0; outsideIdx < is.numDomainNeighbors(); ++outsideIdx)
            {
                const auto& outside = is.domainEntity(outsideIdx);
                const auto bulkElementIdx = bulkGridGeometry.elementMapper().index(outside);
                lowDimVolumeInBulkElement_[bulkElementIdx] += intersectionGeometry.volume()*M_PI*radius*radius;
            }
        }
    }

    //! vector for the volume fraction of the lowdim domain in the bulk domain cells
    std::vector<Scalar> lowDimVolumeInBulkElement_;

    static constexpr std::size_t invalidValue_ = std::numeric_limits<std::size_t>::max();
    //! the index of the radius in the values of the lower-dimensional elements of the parent
    std::size_t radiusValue_ = invalidValue_;
};

} // end namespace Dumux

#endif
