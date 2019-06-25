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
 * \ingroup CCMpfaDiscretization
 * \brief Data handle class for interaction volumes of mpfa methods.
 *        This class is passed to interaction volumes to store the necessary data in it.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_INTERACTIONVOLUMEDATAHANDLE_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_INTERACTIONVOLUMEDATAHANDLE_HH

#include <cassert>
#include <vector>

#include <dune/common/dynvector.hh>

#include <dumux/common/parameters.hh>

namespace Dumux {
namespace CCMpfaDataHandleBases {

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Common base class to all handles. Stores arrays of the
 *        matrices involved in the interaction volume-local systems
 *        of equations. Apart from the transmissibility matrix we
 *        store those matrices that are needed e.g. for later face
 *        pressure reconstruction.
 *        The fluxes as well as the local systems of equations can
 *        be expressed as functions of the intermediate unknown face
 *        face values \f$\bar{\mathbf{u}}\f$ and the known cell/Dirichlet
 *        values \f$\mathbf{u}\f$ using the matrices \f$\mathbf{A}\f$,
 *        \f$\mathbf{B}\f$, \f$\mathbf{C}\f$, \f$\mathbf{D}\f$ and the
 *        vector of Neumann fluxes \f$\mathbf{N}\f$ as follows:
 *
 *        Fluxes: \f$\mathbf{f} = \mathbf{C}\bar{\mathbf{u}} + \mathbf{D}\mathbf{u}\f$
 *        Local eq system: \f$\mathbf{A}\bar{\mathbf{u}} = \mathbf{B}\mathbf{u} + \mathbf{N}\f$
 *
 * \tparam MVT The matrix/vector traits collecting type information used by the iv
 * \tparam size1 first size specifier for the arrays
 * \tparam size2 second size specifier for the arrays
 */
template<class MVT, int size1, int size2>
class SystemMatricesHandle
{
    using AMatrix = typename MVT::AMatrix;
    using BMatrix = typename MVT::BMatrix;
    using CMatrix = typename MVT::CMatrix;
    using TMatrix = typename MVT::TMatrix;
    using CellVector = typename MVT::CellVector;
    using OutsideTij = std::vector< std::vector<CellVector> >;
    using OmegaStorage = typename MVT::OmegaStorage;

public:
    //! Access functions to context-dependent data
    const CMatrix& CA() const { return CA_[contextIdx1_][contextIdx2_]; }
    CMatrix& CA() { return CA_[contextIdx1_][contextIdx2_]; }

    const AMatrix& A() const { return A_[contextIdx1_][contextIdx2_]; }
    AMatrix& A() { return A_[contextIdx1_][contextIdx2_]; }

    const BMatrix& AB() const { return AB_[contextIdx1_][contextIdx2_]; }
    BMatrix& AB() { return AB_[contextIdx1_][contextIdx2_]; }

    const TMatrix& T() const { return T_[contextIdx1_][contextIdx2_]; }
    TMatrix& T() { return T_[contextIdx1_][contextIdx2_]; }

    const OutsideTij& tijOutside() const { return tijOutside_[contextIdx1_][contextIdx2_]; }
    OutsideTij& tijOutside() { return tijOutside_[contextIdx1_][contextIdx2_]; }

    const OmegaStorage& omegas() const { return wijk_[contextIdx1_][contextIdx2_]; }
    OmegaStorage& omegas() { return wijk_[contextIdx1_][contextIdx2_]; }

    //! functionality to set the context indices
    void setContextIndex1(unsigned int idx) const { assert(idx < size1); contextIdx1_ = idx; }
    void setContextIndex2(unsigned int idx) const { assert(idx < size2); contextIdx2_ = idx; }

protected:
    //! indices to be set before accessing data
    mutable unsigned int contextIdx1_{0};
    mutable unsigned int contextIdx2_{0};

    std::array< std::array<OmegaStorage, size2>, size1 > wijk_;     //!< The omega factors that form the entries of the matrices

    std::array< std::array<TMatrix, size2>, size1 > T_;             //!< The transmissibility matrix
    std::array< std::array<AMatrix, size2>, size1 > A_;             //!< Inverse of the iv-local system matrix
    std::array< std::array<BMatrix, size2>, size1 > AB_;            //!< A_ left multiplied to B
    std::array< std::array<CMatrix, size2>, size1 > CA_;            //!< A_ right multiplied to C
    std::array< std::array<OutsideTij, size2>, size1 > tijOutside_; //!< The transmissibilities for "outside" faces (on surface grids)
};

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Common base class to all handles. Stores arrays of the vectors
 *        involved in the interaction volume-local systems of equations.
 *
 * \tparam MVT The matrix/vector traits collecting type information used by the iv
 * \tparam size1 first size specifier for the arrays
 * \tparam size2 second size specifier for the arrays
 */
template<class MVT, int size1, int size2>
class SystemVectorsHandle
{
    using CellVector = typename MVT::CellVector;

public:
    //! Access to the iv-wide known cell/Dirichlet values
    const CellVector& uj() const { return u_[contextIdx1_][contextIdx2_]; }
    CellVector& uj() { return u_[contextIdx1_][contextIdx2_]; }

protected:
    //! functionality to set the context indices
    void setContextIndex1(unsigned int idx) const { assert(idx < size1); contextIdx1_ = idx; }
    void setContextIndex2(unsigned int idx) const { assert(idx < size2); contextIdx2_ = idx; }

    //! indices to be set before accessing data
    mutable unsigned int contextIdx1_{0};
    mutable unsigned int contextIdx2_{0};

    //! The interaction volume-local known values
    std::array< std::array<CellVector, size2>, size1 > u_;
};

} // end namespace CCMpfaDataHandleBases

//! Empty data handle class
class EmptyDataHandle {};

//! Data handle for quantities related to advection
template<class MatVecTraits, class PhysicsTraits, bool EnableAdvection>
class AdvectionDataHandle
: public CCMpfaDataHandleBases::SystemMatricesHandle<MatVecTraits, 1, 1>
, public CCMpfaDataHandleBases::SystemVectorsHandle<MatVecTraits, PhysicsTraits::numPhases, 1>
{
    // we only have one local system for all phases since we
    // solve them w.r.t. permeability tensor (unique for all phases)
    using Base1 = CCMpfaDataHandleBases::SystemMatricesHandle<MatVecTraits, 1, 1>;

    // we do have cell/Dirichlet values for all phases though!
    static constexpr int numPhases = PhysicsTraits::numPhases;
    using Base2 = CCMpfaDataHandleBases::SystemVectorsHandle<MatVecTraits, numPhases, 1>;

    using UnknownVector = typename MatVecTraits::AMatrix::row_type;
    using FaceVector = typename MatVecTraits::FaceVector;
    using FaceScalar = typename FaceVector::value_type;
    using OutsideGravityStorage = std::vector< std::vector<FaceScalar> >;

public:
    //! Set the phase index of the context
    void setPhaseIndex(unsigned int phaseIdx) const { Base2::setContextIndex1(phaseIdx); }

    //! The gravitational flux contributions for a phase on all faces
    const FaceVector& g() const { return g_[Base2::contextIdx1_]; }
    FaceVector& g() { return g_[Base2::contextIdx1_]; }

    //! The deltaG vector for gravity within the iv-local eq-system
    const UnknownVector& deltaG() const { return deltaG_[Base2::contextIdx1_]; }
    UnknownVector& deltaG() { return deltaG_[Base2::contextIdx1_]; }

    //! The gravitational acceleration for one phase on "outside" faces (used on surface grids)
    const OutsideGravityStorage& gOutside() const { return outsideG_[Base2::contextIdx1_]; }
    OutsideGravityStorage& gOutside() { return outsideG_[Base2::contextIdx1_]; }

private:
    std::array< FaceVector, numPhases > g_; //!< The gravitational acceleration at each scvf (only for enabled gravity)
    std::array< UnknownVector, numPhases > deltaG_; //!< The gravity coefficients forming part of iv-local eq-system
    std::array< OutsideGravityStorage, numPhases > outsideG_;  //!< The gravitational acceleration on "outside" faces (only on surface grids)
};

//! Data handle for quantities related to diffusion
template<class MatVecTraits, class PhysicsTraits, bool EnableDiffusion>
class DiffusionDataHandle
: public CCMpfaDataHandleBases::SystemMatricesHandle<MatVecTraits, PhysicsTraits::numPhases, PhysicsTraits::numComponents>
, public CCMpfaDataHandleBases::SystemVectorsHandle<MatVecTraits, PhysicsTraits::numPhases, PhysicsTraits::numComponents>
{
    static constexpr int numPhases = PhysicsTraits::numPhases;
    static constexpr int numComponents = PhysicsTraits::numComponents;
    using Base1 = CCMpfaDataHandleBases::SystemMatricesHandle<MatVecTraits, numPhases, numComponents>;
    using Base2 = CCMpfaDataHandleBases::SystemVectorsHandle<MatVecTraits, numPhases, numComponents>;

public:
    //! diffusion caches need to set phase and component index
    void setPhaseIndex(unsigned int phaseIdx) const
    { Base1::setContextIndex1(phaseIdx); Base2::setContextIndex1(phaseIdx); }
    void setComponentIndex(unsigned int compIdx) const
    { Base1::setContextIndex2(compIdx); Base2::setContextIndex2(compIdx); }
};

//! Data handle for quantities related to heat conduction
template<class MatVecTraits, class PhysicsTraits, bool enableHeatConduction>
class HeatConductionDataHandle
: public CCMpfaDataHandleBases::SystemMatricesHandle<MatVecTraits, 1, 1>
, public CCMpfaDataHandleBases::SystemVectorsHandle<MatVecTraits, 1, 1>
{};

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
{

public:
    //! export the underlying process-specific handle types
    using AdvectionHandle = AdvectionDataHandle<MVT, PT, PT::enableAdvection>;
    using DiffusionHandle = DiffusionDataHandle<MVT, PT, PT::enableMolecularDiffusion>;
    using HeatConductionHandle = HeatConductionDataHandle<MVT, PT, PT::enableHeatConduction>;

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
