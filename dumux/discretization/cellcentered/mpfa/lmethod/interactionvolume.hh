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
 * \brief Base classes for interaction volume of mpfa methods.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_L_INTERACTIONVOLUME_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_L_INTERACTIONVOLUME_HH

#include <dumux/common/math.hh>
#include <dumux/implicit/cellcentered/mpfa/properties.hh>
#include <dumux/discretization/cellcentered/mpfa/interactionvolumebase.hh>
#include <dumux/discretization/cellcentered/mpfa/omethod/interactionvolume.hh>
#include <dumux/discretization/cellcentered/mpfa/facetypes.hh>
#include <dumux/discretization/cellcentered/mpfa/methods.hh>

#include "interactionvolumeseed.hh"
#include "interactionregions.hh"

namespace Dumux
{
//! Specialization of the interaction volume traits class for the mpfa-o method
template<class TypeTag>
class CCMpfaLInteractionVolumeTraits : public CCMpfaInteractionVolumeTraitsBase<TypeTag>
{
    using BaseTraits = CCMpfaInteractionVolumeTraitsBase<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    using BoundaryInteractionVolume = CCMpfaInteractionVolumeImplementation<TypeTag, MpfaMethods::oMethod>;

    // in the interaction region there will be dim faces and dim+1 cell pressures
    using PositionVector = std::vector<GlobalPosition>;
    using Matrix = Dune::FieldMatrix<Scalar, dim, dim+1>;
    using Vector = Dune::FieldVector<Scalar, dim+1>;

    using typename BaseTraits::LocalIndexSet;
    using typename BaseTraits::GlobalIndexSet;
    using Seed = CCMpfaLInteractionVolumeSeed<GlobalIndexSet, LocalIndexSet>;
};

/*!
 * \ingroup Mpfa
 * \brief Base class for the interaction volumes of the mpfa-l method.
 */
template<class TypeTag>
class CCMpfaInteractionVolumeImplementation<TypeTag, MpfaMethods::lMethod>
       : public CCMpfaInteractionVolumeBase<TypeTag, CCMpfaLInteractionVolumeTraits<TypeTag>>
{
    // The interaction volume implementation has to be friend
    friend typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using Traits = CCMpfaLInteractionVolumeTraits<TypeTag>;
    using ParentType = CCMpfaInteractionVolumeBase<TypeTag, Traits>;

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using MpfaHelper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using InteriorBoundaryData = typename GET_PROP_TYPE(TypeTag, InteriorBoundaryData);

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using LocalBasis = std::array<GlobalPosition, dim>;
    using InteractionRegion = Dumux::InteractionRegion<TypeTag>;

    using CoefficientVector = typename Traits::Vector;
    using Tensor = typename Traits::Tensor;

public:
    using Matrix = typename Traits::Matrix;
    using typename ParentType::GlobalLocalFaceDataPair;
    using typename ParentType::LocalIndexType;
    using typename ParentType::LocalIndexSet;
    using typename ParentType::LocalFaceData;
    using typename ParentType::GlobalIndexType;
    using typename ParentType::GlobalIndexSet;
    using typename ParentType::PositionVector;
    using typename ParentType::Seed;

    using ScvSeedType = typename Seed::LocalScvSeed;
    using OuterScvSeedType = typename Seed::LocalOuterScvSeed;

    CCMpfaInteractionVolumeImplementation(const Seed& seed,
                                          const Problem& problem,
                                          const FVElementGeometry& fvGeometry,
                                          const ElementVolumeVariables& elemVolVars)
    : seedPtr_(&seed),
      problemPtr_(&problem),
      fvGeometryPtr_(&fvGeometry),
      elemVolVarsPtr_(&elemVolVars),
      regionUnique_(seed.isUnique()),
      systemSolved_(false)
    {
        // set up the possible interaction regions
        setupInteractionRegions_(seed, fvGeometry);
    }

    template<typename GetTensorFunction>
    void solveLocalSystem(const GetTensorFunction& getTensor)
    {
        if (regionUnique_)
        {
            const auto& ir = interactionRegions_[0];
            T_ = assembleMatrix_(getTensor, ir);
            interactionRegionId_ = 0;
            updateGlobalLocalFaceData();
            systemSolved_ = true;
        }
        else
        {
            std::array<Matrix, 2> M;
            M[0] = assembleMatrix_(getTensor, interactionRegions_[0]);
            M[1] = assembleMatrix_(getTensor, interactionRegions_[1]);

            interactionRegionId_ = MpfaHelper::selectionCriterion(interactionRegions_[0], interactionRegions_[1], M[0], M[1]);
            updateGlobalLocalFaceData();
            systemSolved_ = true;
            T_ = M[interactionRegionId_];
        }
    }

    //! returns the local index pair. This holds the scvf's local index and a boolean whether or not the flux has to be inverted
    const LocalFaceData getLocalFaceData(const SubControlVolumeFace& scvf) const
    {
        assert(systemSolved_ && "Scvf Indices not set yet. You have to call solveLocalSystem() beforehand.");
        assert((scvf.index() == globalLocalScvfPairedData_[0].first->index()
                || scvf.index() == globalLocalScvfPairedData_[1].first->index())
                   && "The provided scvf is not the flux face of the interaction volume.");

        if (globalLocalScvfPairedData_[0].first->index() == scvf.index())
            return globalLocalScvfPairedData_[0].second;
        else
            return globalLocalScvfPairedData_[1].second;
    }

    //! returns the transmissibilities corresponding to the bound scvf
    CoefficientVector getTransmissibilities(const LocalFaceData& localFaceData) const
    {
        assert(systemSolved_ && "Transmissibilities not calculated yet. You have to call solveLocalSystem() beforehand.");

        auto tij = T_[localFaceData.localScvfIndex];
        if (localFaceData.isOutside)
            tij *= -1.0;
        return tij;
    }

    const GlobalIndexSet& volVarsStencil() const
    {
        assert(systemSolved_ && "stencil not set yet. You have to call solveLocalSystem() beforehand.");
        return interactionRegion().scvIndices;
    }

    const PositionVector& volVarsPositions() const
    {
        assert(systemSolved_ && "stencil not set yet. You have to call solveLocalSystem() beforehand.");
        return interactionRegion().scvCenters;
    }

    const std::vector<GlobalLocalFaceDataPair>& globalLocalScvfPairedData() const
    { return globalLocalScvfPairedData_; }

    //! Boundaries will be treated by a different mpfa method (e.g. o method). Thus, on
    //! faces in l-method interaction volumes there will never be a Neumann flux contribution.
    Scalar getNeumannFlux(const LocalFaceData& localFaceData, unsigned int eqIdx) const
    { return 0.0; }

    //! See comment of getNeumannFlux()
    CoefficientVector getNeumannFluxTransformationCoefficients(const LocalFaceData& localFaceData) const
    { return CoefficientVector(0.0); }

    const Matrix& matrix() const
    { return T_; }

    const InteractionRegion& interactionRegion() const
    { return interactionRegions_[interactionRegionId_]; }

    // The l-method can't treat interior boundaries
    const std::vector<InteriorBoundaryData> interiorBoundaryData() const
    { return std::vector<InteriorBoundaryData>(); }

private:
    //! Assembles and solves the local equation system
    //! Specialization for dim = 2
    template<typename GetTensorFunction, int d = dim>
    typename std::enable_if<d == 2, Matrix>::type
    assembleMatrix_(const GetTensorFunction& getTensor, const InteractionRegion& ir)
    {
        // the elements the scvs live in
        const auto& e1 = ir.elements[0];
        const auto& e2 = ir.elements[1];
        const auto& e3 = ir.elements[2];

        // the corresponding scvs
        const auto scv1 = fvGeometry_().scv(ir.scvIndices[0]);
        const auto scv2 = fvGeometry_().scv(ir.scvIndices[1]);
        const auto scv3 = fvGeometry_().scv(ir.scvIndices[2]);

        // Get diffusion tensors in the three scvs
        const auto T1 = getTensor(problem_(), e1, elemVolVars_()[scv1], fvGeometry_(), scv1);
        const auto T2 = getTensor(problem_(), e2, elemVolVars_()[scv2], fvGeometry_(), scv2);
        const auto T3 = getTensor(problem_(), e3, elemVolVars_()[scv3], fvGeometry_(), scv3);

        // required omega factors
        Scalar w111 = calculateOmega_(ir.normal[0], ir.nu[0], ir.detX[0], T1);
        Scalar w112 = calculateOmega_(ir.normal[0], ir.nu[1], ir.detX[0], T1);
        Scalar w123 = calculateOmega_(ir.normal[0], ir.nu[2], ir.detX[1], T2);
        Scalar w124 = calculateOmega_(ir.normal[0], ir.nu[3], ir.detX[1], T2);

        Scalar w211 = calculateOmega_(ir.normal[1], ir.nu[0], ir.detX[0], T1);
        Scalar w212 = calculateOmega_(ir.normal[1], ir.nu[1], ir.detX[0], T1);
        Scalar w235 = calculateOmega_(ir.normal[1], ir.nu[4], ir.detX[2], T3);
        Scalar w236 = calculateOmega_(ir.normal[1], ir.nu[5], ir.detX[2], T3);

        Scalar xi711 = calculateXi_(ir.nu[6], ir.nu[0], ir.detX[0]);
        Scalar xi712 = calculateXi_(ir.nu[6], ir.nu[1], ir.detX[0]);

        Dune::FieldMatrix<Scalar, dim, dim> C, A;
        Dune::FieldMatrix<Scalar, dim, dim+1> B;

        C[0][0] = -w111;
        C[0][1] = -w112;
        C[1][0] = -w211;
        C[1][1] = -w212;

        A[0][0] = w111 - w124 - w123*xi711;
        A[0][1] = w112 - w123*xi712;
        A[1][0] = w211 - w236*xi711;
        A[1][1] = w212 - w235 - w236*xi712;

        B[0][0] = w111 + w112 + w123*(1 - xi711 -xi712);
        B[0][1] = -w123 - w124;
        B[0][2] = 0.0;
        B[1][0] = w211 + w212 + w236*(1 - xi711 - xi712);
        B[1][1] = 0.0;
        B[1][2] = -w235 - w236;

        // T = C*A^-1*B + D
        A.invert();
        auto T = (A.leftmultiply(C)).rightmultiplyany(B);
        T[0][0] += w111 + w112;
        T[1][0] += w211 + w212;

        return T;
    }

    //! sets up the interaction regions for later transmissibility matrix calculation
    //! Specialization for dim = 2
    template<int d = dim>
    typename std::enable_if<d == 2>::type
    setupInteractionRegions_(const Seed& seed, const FVElementGeometry& fvGeometry)
    {
        const auto& fvGridGeometry = problem_().model().fvGridGeometry();

        interactionRegions_.reserve(2);
        auto&& scvSeed1 = seed.scvSeed(0);
        auto&& scvSeed2 = seed.scvSeed(1);
        auto&& outerScvSeed1 = seed.outerScvSeed(0);
        auto&& outerScvSeed2 = seed.outerScvSeed(1);
        auto e1 = fvGridGeometry.element(scvSeed1.globalIndex());
        auto e2 = fvGridGeometry.element(scvSeed2.globalIndex());
        auto e3 = fvGridGeometry.element(outerScvSeed1.globalIndex());
        auto e4 = fvGridGeometry.element(outerScvSeed2.globalIndex());

        // scvSeed1 is the one the seed construction began at
        if (scvSeed1.contiFaceLocalIdx() == 0)
        {
            interactionRegions_.emplace_back(problem_(), fvGeometry, scvSeed1, OuterScvSeedType(scvSeed2), outerScvSeed1, e1, e2, e3);
            interactionRegions_.emplace_back(problem_(), fvGeometry, scvSeed2, outerScvSeed2, OuterScvSeedType(scvSeed1), e2, e4, e1);
        }
        else
        {
            interactionRegions_.emplace_back(problem_(), fvGeometry, scvSeed1, outerScvSeed1, OuterScvSeedType(scvSeed2), e1, e3, e2);
            interactionRegions_.emplace_back(problem_(), fvGeometry, scvSeed2, OuterScvSeedType(scvSeed1), outerScvSeed2, e2, e1, e4);
        }

        // as an initial guess we set the first interaction region as the one used
        // this is updated after solution of the local system
        interactionRegionId_ = 0;
        updateGlobalLocalFaceData();
    }

    //! set the global&local face data according to the current choice for the interaction region
    void updateGlobalLocalFaceData()
    {
        const auto numPairs = globalLocalScvfPairedData_.size();

        // if no data has been stored yet, initialize
        if (numPairs == 0)
        {
            globalLocalScvfPairedData_.clear();
            globalLocalScvfPairedData_.reserve(2);
            globalLocalScvfPairedData_.emplace_back(&fvGeometry_().scvf(seed_().globalScvfIndices()[0]),
                                                    LocalFaceData(interactionRegion().contiFaceLocalIdx, /*dummy*/0, false));
            globalLocalScvfPairedData_.emplace_back(&fvGeometry_().scvf(seed_().globalScvfIndices()[1]),
                                                    LocalFaceData(interactionRegion().contiFaceLocalIdx, /*dummy*/0, true));
        }
        // if order was swapped, change bools indicating which face is the "outside" version of the local face
        else if (globalLocalScvfPairedData_[1].first->index() == interactionRegion().globalScvfIndices[0])
        {
            globalLocalScvfPairedData_[0].second.localScvfIndex = interactionRegion().contiFaceLocalIdx;
            globalLocalScvfPairedData_[1].second.localScvfIndex = interactionRegion().contiFaceLocalIdx;
            globalLocalScvfPairedData_[0].second.isOutside = true;
            globalLocalScvfPairedData_[1].second.isOutside = false;
        }
    }

    Scalar calculateOmega_(const GlobalPosition& normal,
                           const GlobalPosition& nu,
                           const Scalar detX,
                           const Tensor& T) const
    {
        GlobalPosition tmp;
        T.mv(nu, tmp);
        return (tmp*normal)/detX;
    }

    // calculates the omega factors entering the local matrix
    Scalar calculateOmega_(const GlobalPosition& normal,
                           const GlobalPosition& nu,
                           const Scalar detX,
                           const Scalar t) const
    {
        // make sure we have positive diffusion coefficients
        assert(t > 0.0 && "non-positive diffusion coefficients cannot be handled by mpfa methods");
        return (normal*nu)*t/detX;
    }

    //! Specialization for dim = 2
    template<int d = dim>
    typename std::enable_if<d == 2, Scalar>::type
    calculateXi_(const GlobalPosition& nu1,
                 const GlobalPosition& nu2,
                 const Scalar detX) const
    {
        LocalBasis basis({nu1, nu2});
        return MpfaHelper::calculateDetX(basis)/detX;
    }

    const Seed& seed_() const
    { return *seedPtr_; }

    const Problem& problem_() const
    { return *problemPtr_; }

    const FVElementGeometry& fvGeometry_() const
    { return *fvGeometryPtr_; }

    const ElementVolumeVariables& elemVolVars_() const
    { return *elemVolVarsPtr_; }

    const Seed* seedPtr_;
    const Problem* problemPtr_;
    const FVElementGeometry* fvGeometryPtr_;
    const ElementVolumeVariables* elemVolVarsPtr_;

    bool regionUnique_;
    bool systemSolved_;

    LocalIndexType interactionRegionId_;
    std::vector<InteractionRegion> interactionRegions_;
    std::vector<GlobalLocalFaceDataPair> globalLocalScvfPairedData_;

    Matrix T_;
};

} // end namespace

#endif
