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
 * \ingroup CCWMpfaDiscretization
 * \brief Data handle class for face data of wmpfa methods.
 */
#ifndef DUMUX_DISCRETIZATION_CC_WMPFA_FACEDATAHANDLE_HH
#define DUMUX_DISCRETIZATION_CC_WMPFA_FACEDATAHANDLE_HH

#include <cassert>
#include <vector>

#include <dumux/common/parameters.hh>
#include "dumux/common/vectordecomposition.hh"

namespace Dumux {

namespace WMpfaHelper{

template<class Coeff, class Stencil>
void eraseZeros(Coeff& c, Stencil& s)
{
    assert(c.size() == s.size());
    for(std::size_t i=0; i<c.size(); ++i)
    {
        if(std::abs(c[i]) < 1.0e-30)
        {
            c.erase(c.begin()+i);
            s.erase(s.begin()+i);
        }
    }
}

}


namespace WMpfaDataHandle
{
//! Empty data handle class
class EmptyDataHandle {};

//! Data handle for quantities related to advection
template<class IntOp, class PhysicsTraits, bool EnableAdvection>
class AdvectionDataHandle
{
    using Scalar = typename IntOp::Scalar;
    //ToDo get the index types
    using SubControlVolume = typename IntOp::GridGeometry::SubControlVolume;
    using GridIndexType = typename SubControlVolume::Traits::GridIndexType;
    using LocalIndexType = typename SubControlVolume::Traits::LocalIndexType;
    using Position = typename IntOp::Position;

    struct Entry{
        Scalar coefficient;
        GridIndexType index;
        Position position;
    };

public:
    using Entries = std::vector<Entry>;

    using Interpolator = typename IntOp::AdvectionInterpolator;

    static constexpr int numPhases = PhysicsTraits::numPhases;

    void clear()
    {
        entries_.clear();
    }

    void prepare()
    {
        clear();
    }

    template<class Problem, class TF, class EG, class SCVF>
    void decompose(const Problem& problem, const Interpolator& intOp, const TF& tensor, const EG& fvGeometry, const SCVF& scvf)
    {
        boundaryFace_ = scvf.boundary();
        const auto coNormal = mv(tensor(scvf.insideScvIdx()), scvf.unitOuterNormal());
        auto&& [indices, coeff, found] = VectorDecomposition::calculateVectorDecomposition(coNormal, intOp.getDistanceVectors(fvGeometry));

        if(!found)
            DUNE_THROW(Dune::InvalidStateException, "CoNormal decomposition not found");
        else
        {
            WMpfaHelper::eraseZeros(coeff, indices);
            updateEntries(problem, intOp, indices, coeff, fvGeometry, scvf);
        }
    }

    template<class Problem, class Indices, class Coeff, class EG, class SCVF>
    void updateEntries(const Problem& problem, const Interpolator& intOp, const Indices& indices, const Coeff& coeff, const EG& fvGeometry, const SCVF& scvf)
    {
        const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
        entries_.push_back({0.0, scv.dofIndex(), scv.center()});
        for(std::size_t i=0; i<indices.size(); ++i)
        {
            const auto& intData = intOp.getInterpolationData(indices[i]);
            bool skipEntries = false;
            if (fvGeometry.hasBoundaryScvf())
            {
                const auto& scvfJ = fvGeometry.scvf(intData.scvfIdx());
                if(scvfJ.boundary())
                {
                    const auto& bcTypes = problem.boundaryTypes(fvGeometry.gridGeometry().element(scvfJ.insideScvIdx()), scvfJ);
                    if(bcTypes.hasNeumann())
                        skipEntries = true;
                }
            }

            if(!skipEntries)
            {
                for(const auto& e : intData.entries())
                {
                    if(e.dofIndex() != scvf.insideScvIdx())
                    {
                        auto c = coeff[i]*e.weight();
                        entries_[0].coefficient += c;
                        //ToDo pos is not needed for each entry!!!
                        entries_.push_back({-1.0*c, e.dofIndex(), intData.position()});
                    }
                }
            }
        }
        scvfIndices_ = scvf.index();
        scvIndices_ = scvf.insideScvIdx();
    }

    const auto& subFluxData() const
    {
        return entries_;
    }

    bool valid()
    {
        return entries_.size() > 0;
    }

private:
    GridIndexType scvfIndices_;
    GridIndexType scvIndices_;
    Entries entries_;
    bool boundaryFace_ = {false};
};

//! Data handle for quantities related to diffusion
template<class IntOp, class PhysicsTraits, bool EnableDiffusion>
class DiffusionDataHandle
{
    static constexpr int numPhases = PhysicsTraits::numPhases;
    static constexpr int numComponents = PhysicsTraits::numComponents;

    using Interpolator = typename IntOp::DiffusionInterpolator;

public:

    template<class TF, class FVElementGeometry, class SubControlVolumeFace>
    void decompose(const Interpolator& intOp, const TF& tensor, const FVElementGeometry& fvGeometry, const SubControlVolumeFace& scvf)
    {}
private:
    using Scalar = typename IntOp::Scalar;
    using OneSidedCoefficients = std::vector<Scalar>;
    using Coefficients = std::pair<OneSidedCoefficients,OneSidedCoefficients>;
    std::array<Coefficients, numPhases * numComponents> coNormalDecomposition_;
};

//! Data handle for quantities related to heat conduction
template<class IntOp, class PhysicsTraits, bool EnableConduction>
class HeatConductionDataHandle
{
    using Interpolator = typename IntOp::HeatConductionInterpolator;

public:

    template<class TF, class FVElementGeometry, class SubControlVolumeFace>
    void decompose(const Interpolator& intOp, const TF& tensor, const FVElementGeometry& fvGeometry, const SubControlVolumeFace& scvf)
    {}

private:
    using Scalar = typename IntOp::Scalar;
    using OneSidedStencil = std::vector<std::size_t>;
    using Stencils = std::pair<OneSidedStencil,OneSidedStencil>;
    using OneSidedCoefficients = std::vector<Scalar>;
    using Coefficients = std::pair<OneSidedCoefficients,OneSidedCoefficients>;
    std::array<Coefficients, 1 > coefficients_;
    std::array<Stencils, 1 > stencils;
};


//! Process-dependent data handles when related process is disabled
template<class IntOp, class PhysicsTraits>
class AdvectionDataHandle<IntOp, PhysicsTraits, false> : public EmptyDataHandle {};
template<class IntOp, class PhysicsTraits>
class DiffusionDataHandle<IntOp, PhysicsTraits, false> : public EmptyDataHandle {};
template<class IntOp, class PhysicsTraits>
class HeatConductionDataHandle<IntOp, PhysicsTraits, false> : public EmptyDataHandle {};

}

template<class IntOp, class PT>
class FaceDataHandle
{

public:
    //! export the underlying process-specific handle types
    using AdvectionHandle = WMpfaDataHandle::AdvectionDataHandle<IntOp, PT, PT::enableAdvection>;
    using DiffusionHandle = WMpfaDataHandle::DiffusionDataHandle<IntOp, PT, PT::enableMolecularDiffusion>;
    using HeatConductionHandle = WMpfaDataHandle::HeatConductionDataHandle<IntOp, PT, PT::enableHeatConduction>;

    //! return references to the handle containing data related to advection
    const AdvectionHandle& advectionHandle() const { return advectionHandle_; }
    AdvectionHandle& advectionHandle() { return advectionHandle_; }

    //! return references to the handle containing data related to diffusion
    const DiffusionHandle& diffusionHandle() const { return diffusionHandle_; }
    DiffusionHandle& diffusionHandle() { return diffusionHandle_; }

    //! return references to the handle containing data related to heat conduction
    const HeatConductionHandle& heatConductionHandle() const { return heatConductionHandle_; }
    HeatConductionHandle& heatConductionHandle() { return heatConductionHandle_; }

private:
    AdvectionHandle advectionHandle_;
    DiffusionHandle diffusionHandle_;
    HeatConductionHandle heatConductionHandle_;
};

} // end namespace Dumux

#endif
