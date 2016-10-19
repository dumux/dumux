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
 * \brief Base class for the flux variables
 */
#ifndef DUMUX_POROUSMEDIUM_IMPLICIT_FLUXVARIABLESCACHE_HH
#define DUMUX_POROUSMEDIUM_IMPLICIT_FLUXVARIABLESCACHE_HH

#include <dumux/implicit/properties.hh>
#include <dumux/implicit/cellcentered/mpfa/properties.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{
// forward declaration
template<class TypeTag, DiscretizationMethods Method>
class PorousMediumFluxVariablesCacheImplementation
{};

/*!
 * \ingroup ImplicitModel
 * \brief The flux variables cache classes for porous media.
 *        Store flux stencils and data required for flux calculation
 */
template<class TypeTag>
using PorousMediumFluxVariablesCache = PorousMediumFluxVariablesCacheImplementation<TypeTag, GET_PROP_VALUE(TypeTag, DiscretizationMethod)>;

// specialization for the Box Method
template<class TypeTag>
class PorousMediumFluxVariablesCacheImplementation<TypeTag, DiscretizationMethods::Box>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::vector<IndexType>;
    using TransmissibilityVector = std::vector<IndexType>;

    using CoordScalar = typename GridView::ctype;
    static const int dim = GridView::dimension;

    using FeCache = Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1>;
    using FeLocalBasis = typename FeCache::FiniteElementType::Traits::LocalBasisType;
    using ShapeJacobian = typename FeLocalBasis::Traits::JacobianType;
    using ShapeValue = typename Dune::FieldVector<Scalar, 1>;
    using JacobianInverseTransposed = typename Element::Geometry::JacobianInverseTransposed;

public:

    void update(const Problem& problem,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const SubControlVolumeFace &scvf)
    {
        const auto geometry = element.geometry();
        const auto& localBasis = fvGeometry.feLocalBasis();

        // evaluate shape functions and gradients at the integration point
        const auto ipLocal = geometry.local(scvf.center());
        jacInvT_ = geometry.jacobianInverseTransposed(ipLocal);
        localBasis.evaluateJacobian(ipLocal, shapeJacobian_);
        //localBasis.evaluateFunction(ipLocal, shapeValue_); // do we need the shapeValues for the flux?

        // The stencil info is obsolete for the box method.
        // It is here for compatibility with cc methods
        stencil_ = Stencil(0);
    }

    const std::vector<ShapeJacobian>& shapeJacobian() const
    { return shapeJacobian_; }

   /* const std::vector<ShapeValue>& shapeValue() const
    { return shapeValue_; }*/

    const JacobianInverseTransposed& jacInvT() const
    { return jacInvT_; }

    const Stencil& stencil() const
    {
        return stencil_;
    }

private:
    std::vector<ShapeJacobian> shapeJacobian_;
    //std::vector<ShapeValue> shapeValue_;
    JacobianInverseTransposed jacInvT_;

    Stencil stencil_;
};

// specialization for the cell centered tpfa method
template<class TypeTag>
class PorousMediumFluxVariablesCacheImplementation<TypeTag, DiscretizationMethods::CCTpfa>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::vector<IndexType>;

public:
    void update(const Problem& problem,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars,
                const SubControlVolumeFace &scvf)
    {
        FluxVariables fluxVars;
        stencil_ = fluxVars.computeStencil(problem, element, fvGeometry, scvf);
        tij_ = AdvectionType::calculateTransmissibilities(problem, element, fvGeometry, elemVolVars, scvf);
    }

    const Stencil& stencil() const
    { return stencil_; }

    const Scalar& tij() const
    { return tij_; }

private:
    Stencil stencil_;
    Scalar tij_;
};

// forward declaration of the base class of the mpfa flux variables cache
template<class TypeTag, bool EnableAdvection, bool EnableMolecularDiffusion, bool EnableEnergyBalance>
class PorousMediumMpfaFluxVariablesCache
{};

// specialization for cell centered mpfa methods
template<class TypeTag>
class PorousMediumFluxVariablesCacheImplementation<TypeTag, DiscretizationMethods::CCMpfa>
       : public PorousMediumMpfaFluxVariablesCache<TypeTag,
                                                   GET_PROP_VALUE(TypeTag, EnableAdvection),
                                                   GET_PROP_VALUE(TypeTag, EnableMolecularDiffusion),
                                                   GET_PROP_VALUE(TypeTag, EnableEnergyBalance)>
{};

// specialization for the case of pure advection
template<class TypeTag>
class PorousMediumMpfaFluxVariablesCache<TypeTag, true, false, false>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;

    static const int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
    static const int dim = GridView::dimension;

    using GlobalPosition = Dune::FieldVector<Scalar, dim>;

    // We always use the dynamic types here to be compatible on the boundary
    using Stencil = typename BoundaryInteractionVolume::GlobalIndexSet;
    using TransmissibilityVector = typename BoundaryInteractionVolume::Vector;
    using PositionVector = typename BoundaryInteractionVolume::PositionVector;

public:
    // the constructor
    PorousMediumMpfaFluxVariablesCache() : isUpdated_(false)
    {
        for (auto& nFlux : phaseNeumannFluxes_)
            nFlux = 0.0;
    }

    // update cached objects
    void updateAdvection(const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const SubControlVolumeFace &scvf,
                         const InteractionVolume& interactionVolume)
    {
        const auto& localIndexPair = interactionVolume.getLocalIndexPair(scvf);
        const auto& volVarsStencil = interactionVolume.volVarsStencil();
        const auto& volVarsPositions = interactionVolume.volVarsPositions();

        // the types coming from the inner interaction volumes might differ (thus, = assignment is not possible)
        const auto numVolVars = volVarsStencil.size();
        volVarsStencil_.clear();
        volVarsStencil_.reserve(numVolVars);
        volVarsPositions_.clear();
        volVarsPositions_.reserve(numVolVars);
        volVarsStencil_.insert(volVarsStencil_.begin(), volVarsStencil.begin(), volVarsStencil.end());
        volVarsPositions_.insert(volVarsPositions_.begin(), volVarsPositions.begin(), volVarsPositions.end());
        tij_ = interactionVolume.getTransmissibilities(localIndexPair);
    }

    void updatePhaseNeumannFlux(const Problem& problem,
                                const Element& element,
                                const FVElementGeometry& fvGeometry,
                                const ElementVolumeVariables& elemVolVars,
                                const SubControlVolumeFace &scvf,
                                const InteractionVolume& interactionVolume,
                                const unsigned int phaseIdx)
    {
        const auto& localIndexPair = interactionVolume.getLocalIndexPair(scvf);
        phaseNeumannFluxes_[phaseIdx] = interactionVolume.getNeumannFlux(localIndexPair);
    }

    const Stencil& advectionVolVarsStencil(const unsigned int phaseIdx) const
    { return volVarsStencil_; }

    const PositionVector& advectionVolVarsPositions(const unsigned int phaseIdx) const
    { return volVarsPositions_; }

    const TransmissibilityVector& advectionTij(const unsigned int phaseIdx) const
    { return tij_; }

    Scalar advectionNeumannFlux(const unsigned int phaseIdx) const
    { return phaseNeumannFluxes_[phaseIdx]; }

    bool isUpdated() const
    { return isUpdated_; }

    void setUpdated()
    {
        isUpdated_ = true;
    }

private:
    bool isUpdated_;
    Stencil volVarsStencil_;
    PositionVector volVarsPositions_;
    TransmissibilityVector tij_;
    std::array<Scalar, numPhases> phaseNeumannFluxes_;
};

} // end namespace

#endif
