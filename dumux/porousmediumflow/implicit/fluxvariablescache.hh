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
#include <dune/localfunctions/lagrange/pqkfactory.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{
// forward declaration
template<class TypeTag, DiscretizationMethods Method>
class PorousMediumFluxVariablesCacheImplementation;

namespace Properties
{
// forward declaration
NEW_PROP_TAG(NumPhases);
NEW_PROP_TAG(InteriorBoundaryData);
NEW_PROP_TAG(EnableInteriorBoundaries);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! The cache is dependent on the active physical processes (advection, diffusion, heat conduction)
//! For each type of process there is a base cache storing the data required to compute the respective fluxes
//! Specializations of the overall cache are provided for combinations of processes
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*!
 * \ingroup ImplicitModel
 * \brief The flux variables cache classes for porous media.
 *        Store data required for flux calculation. For each type of physical process (advection, diffusion, heat conduction)
 *        there is a base cache storing the data required to compute the respective fluxes. Specializations of the overall
 *        cache class are provided for different combinations of processes.
 */
template<class TypeTag>
using PorousMediumFluxVariablesCache = PorousMediumFluxVariablesCacheImplementation<TypeTag, GET_PROP_VALUE(TypeTag, DiscretizationMethod)>;

//! We only store discretization-related quantities for the box method.
//! Thus, we need no physics-dependent specialization.
template<class TypeTag>
class PorousMediumFluxVariablesCacheImplementation<TypeTag, DiscretizationMethods::Box>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
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
                const ElementVolumeVariables& elemVolVars,
                const SubControlVolumeFace &scvf)
    {
        const auto geometry = element.geometry();
        const auto& localBasis = fvGeometry.feLocalBasis();

        // evaluate shape functions and gradients at the integration point
        const auto ipLocal = geometry.local(scvf.center());
        jacInvT_ = geometry.jacobianInverseTransposed(ipLocal);
        localBasis.evaluateJacobian(ipLocal, shapeJacobian_);
        localBasis.evaluateFunction(ipLocal, shapeValues_); // do we need the shapeValues for the flux?
    }

    const std::vector<ShapeJacobian>& shapeJacobian() const
    { return shapeJacobian_; }

    const std::vector<ShapeValue>& shapeValues() const
    { return shapeValues_; }

    const JacobianInverseTransposed& jacInvT() const
    { return jacInvT_; }

private:
    std::vector<ShapeJacobian> shapeJacobian_;
    std::vector<ShapeValue> shapeValues_;
    JacobianInverseTransposed jacInvT_;
};

// forward declaration of the base class of the tpfa flux variables cache
template<class TypeTag, bool EnableAdvection, bool EnableMolecularDiffusion, bool EnableEnergyBalance>
class CCTpfaPorousMediumFluxVariablesCache;

// specialization for the cell centered tpfa method
template<class TypeTag>
class PorousMediumFluxVariablesCacheImplementation<TypeTag, DiscretizationMethods::CCTpfa>
      : public CCTpfaPorousMediumFluxVariablesCache<TypeTag, GET_PROP_VALUE(TypeTag, EnableAdvection),
                                                             GET_PROP_VALUE(TypeTag, EnableMolecularDiffusion),
                                                             GET_PROP_VALUE(TypeTag, EnableEnergyBalance)> {};

// specialization for the case of pure advection
template<class TypeTag>
class CCTpfaPorousMediumFluxVariablesCache<TypeTag, true, false, false> : public GET_PROP_TYPE(TypeTag, AdvectionType)::Cache {};

// specialization for the case of advection & diffusion
template<class TypeTag>
class CCTpfaPorousMediumFluxVariablesCache<TypeTag, true, true, false> : public GET_PROP_TYPE(TypeTag, AdvectionType)::Cache,
                                                                         public GET_PROP_TYPE(TypeTag, MolecularDiffusionType)::Cache {};

// specialization for the case of advection & heat conduction
template<class TypeTag>
class CCTpfaPorousMediumFluxVariablesCache<TypeTag, true, false, true> : public GET_PROP_TYPE(TypeTag, AdvectionType)::Cache,
                                                                         public GET_PROP_TYPE(TypeTag, HeatConductionType)::Cache {};

// specialization for the case of advection, diffusion & heat conduction
template<class TypeTag>
class CCTpfaPorousMediumFluxVariablesCache<TypeTag, true, true, true> : public GET_PROP_TYPE(TypeTag, AdvectionType)::Cache,
                                                                        public GET_PROP_TYPE(TypeTag, MolecularDiffusionType)::Cache,
                                                                        public GET_PROP_TYPE(TypeTag, HeatConductionType)::Cache {};

// TODO further specializations

// forward declaration of the base class of the mpfa flux variables cache
template<class TypeTag, bool EnableAdvection, bool EnableMolecularDiffusion, bool EnableEnergyBalance>
class CCMpfaPorousMediumFluxVariablesCache;

// specialization of the flux variables cache for the cell centered finite volume mpfa scheme
template<class TypeTag>
class PorousMediumFluxVariablesCacheImplementation<TypeTag, DiscretizationMethods::CCMpfa>
      : public CCMpfaPorousMediumFluxVariablesCache<TypeTag, GET_PROP_VALUE(TypeTag, EnableAdvection),
                                                             GET_PROP_VALUE(TypeTag, EnableMolecularDiffusion),
                                                             GET_PROP_VALUE(TypeTag, EnableEnergyBalance)>
{
    using InteriorBoundaryData = typename GET_PROP_TYPE(TypeTag, InteriorBoundaryData);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ParentType = CCMpfaPorousMediumFluxVariablesCache<TypeTag, GET_PROP_VALUE(TypeTag, EnableAdvection),
                                                                     GET_PROP_VALUE(TypeTag, EnableMolecularDiffusion),
                                                                     GET_PROP_VALUE(TypeTag, EnableEnergyBalance)>;

    static constexpr bool enableInteriorBoundaries = GET_PROP_VALUE(TypeTag, EnableInteriorBoundaries);

public:
    //! the constructor
    PorousMediumFluxVariablesCacheImplementation()
    : ParentType(),
      isUpdated_(false),
      interiorBoundaryInfoSelf_(false, -1)
    {}

    //! Returns whether or not this cache has been updated
    bool isUpdated() const
    { return isUpdated_; }

    //! Sets the update status from outside. This is used to only update
    //! the cache once after solving the local system. When visiting an scvf
    //! of the same interaction region again, the update is skipped.
    void setUpdateStatus(const bool status)
    {
        isUpdated_ = status;
    }

    //! maybe update data on interior Dirichlet boundaries
    template<class InteractionVolume>
    void updateInteriorBoundaryData(const InteractionVolume& iv, const SubControlVolumeFace &scvf)
    {
        if (enableInteriorBoundaries)
        {
            interiorBoundaryData_ = iv.interiorBoundaryData();

            // check if the actual scvf is on an interior Dirichlet boundary
            const auto scvfIdx = scvf.index();
            unsigned int indexInData = 0;
            for (auto&& data : interiorBoundaryData_)
            {
                if (data.scvfIndex() == scvfIdx)
                {
                    interiorBoundaryInfoSelf_ = std::make_pair(true, indexInData);
                    break;
                }

                indexInData++;
            }
        }
    }

    bool isInteriorBoundary() const
    { return interiorBoundaryInfoSelf_.first; }

    const std::vector<InteriorBoundaryData>& interiorBoundaryData() const
    { return interiorBoundaryData_; }

    const InteriorBoundaryData& interiorBoundaryDataSelf() const
    {
        assert(interiorBoundaryInfoSelf_.first && "Trying to obtain interior boundary data on a face that is not marked as such");
        assert(interiorBoundaryInfoSelf_.second != -1 && "The index to the interior boundary data of this face has not been set");
        return interiorBoundaryData_[interiorBoundaryInfoSelf_.second];
    }

private:
    // indicates whether or not this cache has been fully updated
    bool isUpdated_;

    // if this face is an interior Dirichlet boundary itself, store additional data
    std::pair<bool, unsigned int> interiorBoundaryInfoSelf_;

    // contains all the interior Dirichlet boundary data within the stencil of this face
    std::vector<InteriorBoundaryData> interiorBoundaryData_;
};

// specialization for the case of pure advection
template<class TypeTag>
class CCMpfaPorousMediumFluxVariablesCache<TypeTag, true, false, false> : public GET_PROP_TYPE(TypeTag, AdvectionType)::Cache {};

// specialization for the case of advection & diffusion
template<class TypeTag>
class CCMpfaPorousMediumFluxVariablesCache<TypeTag, true, true, false> : public GET_PROP_TYPE(TypeTag, AdvectionType)::Cache,
                                                                         public GET_PROP_TYPE(TypeTag, MolecularDiffusionType)::Cache {};

// specialization for the case of advection & heat conduction
template<class TypeTag>
class CCMpfaPorousMediumFluxVariablesCache<TypeTag, true, false, true> : public GET_PROP_TYPE(TypeTag, AdvectionType)::Cache,
                                                                         public GET_PROP_TYPE(TypeTag, HeatConductionType)::Cache {};

// specialization for the case of advection, diffusion & heat conduction
template<class TypeTag>
class CCMpfaPorousMediumFluxVariablesCache<TypeTag, true, true, true> : public GET_PROP_TYPE(TypeTag, AdvectionType)::Cache,
                                                                        public GET_PROP_TYPE(TypeTag, MolecularDiffusionType)::Cache,
                                                                        public GET_PROP_TYPE(TypeTag, HeatConductionType)::Cache {};

// forward declaration of the base class of the tpfa flux variables cache
template<class TypeTag, bool EnableAdvection, bool EnableMolecularDiffusion, bool EnableEnergyBalance>
class MimeticPorousMediumFluxVariablesCache;

// specialization for the cell centered tpfa method
template<class TypeTag>
class PorousMediumFluxVariablesCacheImplementation<TypeTag, DiscretizationMethods::Mimetic>
      : public MimeticPorousMediumFluxVariablesCache<TypeTag, GET_PROP_VALUE(TypeTag, EnableAdvection),
                                                             GET_PROP_VALUE(TypeTag, EnableMolecularDiffusion),
                                                             GET_PROP_VALUE(TypeTag, EnableEnergyBalance)> {};

// specialization for the case of pure advection
template<class TypeTag>
class MimeticPorousMediumFluxVariablesCache<TypeTag, true, false, false> : public GET_PROP_TYPE(TypeTag, AdvectionType)::Cache {};

// specialization for the case of advection & diffusion
template<class TypeTag>
class MimeticPorousMediumFluxVariablesCache<TypeTag, true, true, false> : public GET_PROP_TYPE(TypeTag, AdvectionType)::Cache,
                                                                         public GET_PROP_TYPE(TypeTag, MolecularDiffusionType)::Cache {};

// specialization for the case of advection & heat conduction
template<class TypeTag>
class MimeticPorousMediumFluxVariablesCache<TypeTag, true, false, true> : public GET_PROP_TYPE(TypeTag, AdvectionType)::Cache,
                                                                         public GET_PROP_TYPE(TypeTag, HeatConductionType)::Cache {};

// specialization for the case of advection, diffusion & heat conduction
template<class TypeTag>
class MimeticPorousMediumFluxVariablesCache<TypeTag, true, true, true> : public GET_PROP_TYPE(TypeTag, AdvectionType)::Cache,
                                                                        public GET_PROP_TYPE(TypeTag, MolecularDiffusionType)::Cache,
                                                                        public GET_PROP_TYPE(TypeTag, HeatConductionType)::Cache {};


// TODO further specializations

} // end namespace

#endif
