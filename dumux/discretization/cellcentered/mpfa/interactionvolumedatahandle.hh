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
 * \brief Data handle class for interaction volumes of mpfa methods.
 *        This class is passed to interaction volumes to store the necessary data in it.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_INTERACTIONVOLUMEDATAHANDLE_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_INTERACTIONVOLUMEDATAHANDLE_HH

#include <vector>

#include <dune/common/dynvector.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/matrixvectorhelper.hh>

namespace Dumux
{

//! Empty data handle class
class EmptyDataHandle
{
public:
    template< class InteractionVolume >
    void resize(const InteractionVolume& iv) {}
};

//! Data handle for quantities related to advection
template<class MatVecTraits, class PhysicsTraits, bool EnableAdvection>
class AdvectionDataHandle
{
    // obtain matrix & vector types from interaction volume
    using AMatrix = typename MatVecTraits::AMatrix;
    using CMatrix = typename MatVecTraits::CMatrix;
    using TMatrix = typename MatVecTraits::TMatrix;
    using CellVector = typename MatVecTraits::CellVector;
    using FaceVector = typename MatVecTraits::FaceVector;
    using FaceScalar = typename FaceVector::value_type;
    using OutsideGravityStorage = std::vector< Dune::DynamicVector<FaceScalar> >;

    static constexpr int numPhases = PhysicsTraits::numPhases;

public:

    //! prepare data handle for subsequent fill
    template< class IV >
    void resize(const IV& iv)
    {
        using LocalIndexType = typename IV::Traits::IndexSet::LocalIndexType;

        // resize matrices
        resizeMatrix(T_, iv.numFaces(), iv.numKnowns());
        resizeMatrix(A_, iv.numUnknowns(), iv.numUnknowns());
        resizeMatrix(CA_, iv.numFaces(), iv.numUnknowns());

        // lambdas to resize a vector
        auto resizeToNumFaces = [&iv] (auto& v) { resizeVector(v, iv.numFaces()); };
        auto resizeToNumKnowns = [&iv] (auto& v) { resizeVector(v, iv.numKnowns()); };

        // resize pressure/gravity vector
        std::for_each(p_.begin(), p_.end(), resizeToNumKnowns);
        std::for_each(g_.begin(), g_.end(), resizeToNumFaces);

        // on surface grids, resize additional containers
        static constexpr int dim = IV::Traits::GridView::dimension;
        static constexpr int dimWorld = IV::Traits::GridView::dimensionworld;
        if (dim < dimWorld)
        {
            outsideT_.resize(iv.numFaces());
            std::for_each(outsideG_.begin(), outsideG_.end(), resizeToNumFaces);

            // resize the entries for each face
            for (LocalIndexType i = 0; i < iv.numFaces(); ++i)
            {
                const auto& localScvf = iv.localScvf(i);
                const auto numNeighbors = localScvf.neighboringLocalScvIndices().size()-1;
                outsideT_[i].resize(numNeighbors);

                std::for_each(outsideT_[i].begin(), outsideT_[i].end(), resizeToNumKnowns);
                std::for_each(outsideG_.begin(), outsideG_.end(), [&] (auto& v) { resizeVector(v[i], numNeighbors); });
            }
        }
    }

    //! Access to the iv-wide pressure of one phase
    const CellVector& pressures(unsigned int pIdx) const { return p_[pIdx]; }
    CellVector& pressures(unsigned int pIdx) { return p_[pIdx]; }

    //! The gravitational flux contributions for a phase on all faces
    const FaceVector& gravity(unsigned int pIdx) const { return g_[pIdx]; }
    FaceVector& gravity(unsigned int pIdx) { return g_[pIdx]; }

    //! Access to the gravitational flux contributions for all phases
    const std::array< FaceVector, numPhases >& gravity() const { return g_; }
    std::array< FaceVector, numPhases >& gravity() { return g_; }


    //! Projection matrix for Neumann/gravity contribution computation
    const CMatrix& advectionCA() const { return CA_; }
    CMatrix& advectionCA() { return CA_; }

    //! Inverse of the iv-local system matrix
    const AMatrix& advectionA() const { return A_; }
    AMatrix& advectionA() { return A_; }

    //! The transmissibility matrix (i.e. C*(A^-1)*B + D)
    const TMatrix& advectionT() const { return T_; }
    TMatrix& advectionT() { return T_; }

    //! The transmissibilities for "outside" faces (used on surface grids)
    const std::vector< std::vector<CellVector> >& advectionTout() const { return outsideT_; }
    std::vector< std::vector<CellVector> >& advectionTout() { return outsideT_; }

    //! The gravitational acceleration for "outside" faces (used on surface grids)
    const std::array< OutsideGravityStorage, numPhases >& gravityOutside() const { return outsideG_; }
    std::array< OutsideGravityStorage, numPhases >& gravityOutside() { return outsideG_; }

    //! The gravitational acceleration for one phase on "outside" faces (used on surface grids)
    const OutsideGravityStorage& gravityOutside(unsigned int pIdx) const { return outsideG_[pIdx]; }
    OutsideGravityStorage& gravityOutside(unsigned int pIdx) { return outsideG_[pIdx]; }

private:
    TMatrix T_;  //!< The transmissibilities such that f_i = T_ij*p_j (... + Neumann/gravity contributions)
    AMatrix A_;  //!< Inverse of the iv-local system matrix (needed e.g. for face pressure reconstruction)
    CMatrix CA_; //!< A_ right multiplied to C (needed e.g. for Neumann/gravity contribution computation)
    std::array< CellVector, numPhases > p_;           //!< The interaction volume-wide phase pressures
    std::array< FaceVector, numPhases > g_;           //!< The gravitational acceleration at each scvf (only for enabled gravity)
    std::vector< std::vector<CellVector> > outsideT_; //!< The transmissibilities for "outside" faces (only on surface grids)
    std::array< OutsideGravityStorage, numPhases > outsideG_;  //!< The gravitational acceleration on "outside" faces (only on surface grids)
};

//! Data handle for quantities related to diffusion
template<class MatVecTraits, class PhysicsTraits, bool EnableDiffusion>
class DiffusionDataHandle
{
    using TMatrix = typename MatVecTraits::TMatrix;
    using CellVector = typename MatVecTraits::CellVector;
    using OutsideTContainer = std::vector< std::vector<CellVector> >;

    static constexpr int numPhases = PhysicsTraits::numPhases;
    static constexpr int numComponents = PhysicsTraits::numComponents;

public:
    //! diffusion caches need to set phase and component index
    void setPhaseIndex(unsigned int phaseIdx) { phaseIdx_ = phaseIdx; }
    void setComponentIndex(unsigned int compIdx) { compIdx_ = compIdx; }

    //! prepare data handle for subsequent fill
    template< class InteractionVolume >
    void resize(const InteractionVolume& iv)
    {
        for (unsigned int pIdx = 0; pIdx < numPhases; ++pIdx)
        {
            for (unsigned int cIdx = 0; cIdx < numComponents; ++cIdx)
            {
                // resize transmissibility matrix & mole fraction vector
                resizeMatrix(T_[pIdx][cIdx], iv.numFaces(), iv.numKnowns());
                resizeVector(xj_[pIdx][cIdx], iv.numKnowns());

                // resize outsideTij on surface grids
                using GridView = typename InteractionVolume::Traits::GridView;
                if (int(GridView::dimension) < int(GridView::dimensionworld))
                {
                    outsideT_[pIdx][cIdx].resize(iv.numFaces());

                    using LocalIndexType = typename InteractionVolume::Traits::IndexSet::LocalIndexType;
                    for (LocalIndexType i = 0; i < iv.numFaces(); ++i)
                    {
                        const auto numNeighbors = iv.localScvf(i).neighboringLocalScvIndices().size() - 1;
                        outsideT_[pIdx][cIdx][i].resize(numNeighbors);
                        for (LocalIndexType j = 0; j < numNeighbors; ++j)
                            resizeVector(outsideT_[pIdx][cIdx][i][j], iv.numKnowns());
                    }
                }
            }
        }
    }

    //! Access to the iv-wide mole fractions of a component in one phase
    const CellVector& moleFractions(unsigned int pIdx, unsigned int compIdx) const { return xj_[pIdx][compIdx]; }
    CellVector& moleFractions(unsigned int pIdx, unsigned int compIdx) { return xj_[pIdx][compIdx]; }

    //! The transmissibilities associated with diffusive fluxes
    const TMatrix& diffusionT() const { return T_[phaseIdx_][compIdx_]; }
    TMatrix& diffusionT() { return T_[phaseIdx_][compIdx_]; }

    //! The transmissibilities for "outside" faces (used on surface grids)
    const OutsideTContainer& diffusionTout() const { return outsideT_[phaseIdx_][compIdx_]; }
    OutsideTContainer& diffusionTout() { return outsideT_[phaseIdx_][compIdx_]; }

private:
    //! diffusion-related variables
    unsigned int phaseIdx_;                                             //!< The phase index set for the context
    unsigned int compIdx_;                                              //!< The component index set for the context
    std::array< std::array<TMatrix, numComponents>, numPhases > T_;     //!< The transmissibilities such that f_i = T_ij*x_j
    std::array< std::array<CellVector, numComponents>, numPhases > xj_; //!< The interaction volume-wide mole fractions
    std::array< std::array<OutsideTContainer, numComponents>, numPhases> outsideT_;
};

//! Data handle for quantities related to heat conduction
template<class MatVecTraits, class PhysicsTraits, bool enableHeatConduction>
class HeatConductionDataHandle
{
    using TMatrix = typename MatVecTraits::TMatrix;
    using CellVector = typename MatVecTraits::CellVector;

public:
    //! prepare data handle for subsequent fill
    template< class InteractionVolume >
    void resize(const InteractionVolume& iv)
    {
        //! resize transmissibility matrix & temperature vector
        resizeMatrix(T_, iv.numFaces(), iv.numKnowns());
        resizeVector(Tj_, iv.numKnowns());

        //! resize outsideTij on surface grids
        using GridView = typename InteractionVolume::Traits::GridView;
        if (int(GridView::dimension) < int(GridView::dimensionworld))
        {
            outsideT_.resize(iv.numFaces());

            using LocalIndexType = typename InteractionVolume::Traits::IndexSet::LocalIndexType;
            for (LocalIndexType i = 0; i < iv.numFaces(); ++i)
            {
                const auto numNeighbors = iv.localScvf(i).neighboringLocalScvIndices().size() - 1;
                outsideT_[i].resize(numNeighbors);
                for (LocalIndexType j = 0; j < numNeighbors; ++j)
                    resizeVector(outsideT_[i][j], iv.numKnowns());
            }
        }
    }

    //! Access to the iv-wide temperatures
    const CellVector& temperatures() const { return Tj_; }
    CellVector& temperatures() { return Tj_; }

    //! The transmissibilities associated with conductive fluxes
    const TMatrix& heatConductionT() const { return T_; }
    TMatrix& heatConductionT() { return T_; }

    //! The transmissibilities for "outside" faces (used on surface grids)
    const std::vector< std::vector<CellVector> >& heatConductionTout() const { return outsideT_; }
    std::vector< std::vector<CellVector> >& heatConductionTout() { return outsideT_; }

private:
    // heat conduction-related variables
    TMatrix T_;                                       //!< The transmissibilities such that f_i = T_ij*T_j
    CellVector Tj_;                                   //!< The interaction volume-wide temperatures
    std::vector< std::vector<CellVector> > outsideT_; //!< The transmissibilities for "outside" faces (only necessary on surface grids)
};

//! Process-dependent data handles when related process is disabled
template<class MatVecTraits, class PhysicsTraits>
class AdvectionDataHandle<MatVecTraits, PhysicsTraits, false> : public EmptyDataHandle {};
template<class MatVecTraits, class PhysicsTraits>
class DiffusionDataHandle<MatVecTraits, PhysicsTraits, false> : public EmptyDataHandle {};
template<class MatVecTraits, class PhysicsTraits>
class HeatConductionDataHandle<MatVecTraits, PhysicsTraits, false> : public EmptyDataHandle {};

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Class for the interaction volume data handle.
 *
 * \tparam MVT The matrix/vector traits collecting type information used by the iv
 * \tparam PT The physics traits collecting data on the physical processes to be considered
 */
template<class MVT, class PT>
class InteractionVolumeDataHandle
: public AdvectionDataHandle<MVT, PT, PT::enableAdvection>
, public DiffusionDataHandle<MVT, PT, PT::enableMolecularDiffusion>
, public HeatConductionDataHandle<MVT, PT, PT::enableHeatConduction>
{
    using AdvectionHandle = AdvectionDataHandle<MVT, PT, PT::enableAdvection>;
    using DiffusionHandle = DiffusionDataHandle<MVT, PT, PT::enableMolecularDiffusion>;
    using HeatConductionHandle = HeatConductionDataHandle<MVT, PT, PT::enableHeatConduction>;

public:
    //! prepare data handles for subsequent fills
    template< class InteractionVolume >
    void resize(const InteractionVolume& iv)
    {
        AdvectionHandle::resize(iv);
        DiffusionHandle::resize(iv);
        HeatConductionHandle::resize(iv);
    }
};

} // end namespace Dumux

#endif
