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
    static const int dimWorld = GridView::dimensionworld;

    using FeCache = Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1>;
    using FeLocalBasis = typename FeCache::FiniteElementType::Traits::LocalBasisType;
    using ShapeJacobian = typename FeLocalBasis::Traits::JacobianType;
    using ShapeValue = typename Dune::FieldVector<Scalar, 1>;
    using JacobianInverseTransposed = typename Element::Geometry::JacobianInverseTransposed;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

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
        localBasis.evaluateFunction(ipLocal, shapeValues_); // shape values for rho

        // compute the gradN at for every scv/dof
        gradN_.resize(fvGeometry.numScv());
        for (const auto& scv: scvs(fvGeometry))
            jacInvT_.mv(shapeJacobian_[scv.indexInElement()][0], gradN_[scv.indexInElement()]);

    }

    const std::vector<ShapeJacobian>& shapeJacobian() const
    { return shapeJacobian_; }

    const std::vector<ShapeValue>& shapeValues() const
    { return shapeValues_; }

    const JacobianInverseTransposed& jacInvT() const
    { return jacInvT_; }

    const GlobalPosition& gradN(unsigned int scvIdxInElement) const
    { return gradN_[scvIdxInElement]; }

private:
    std::vector<GlobalPosition> gradN_;
    std::vector<ShapeJacobian> shapeJacobian_;
    std::vector<ShapeValue> shapeValues_;
    JacobianInverseTransposed jacInvT_;
};

// forward declaration of the base class of the tpfa flux variables cache
template<class TypeTag, bool EnableAdvection, bool EnableMolecularDiffusion, bool EnableEnergyBalance>
class CCTpfaPorousMediumFluxVariablesCache : public GET_PROP_TYPE(TypeTag, AdvectionType)::Cache,
                                             public GET_PROP_TYPE(TypeTag, MolecularDiffusionType)::Cache,
                                             public GET_PROP_TYPE(TypeTag, HeatConductionType)::Cache {};

// specialization for the cell centered tpfa method
template<class TypeTag>
class PorousMediumFluxVariablesCacheImplementation<TypeTag, DiscretizationMethods::CCTpfa>
      : public CCTpfaPorousMediumFluxVariablesCache<TypeTag, GET_PROP_VALUE(TypeTag, EnableAdvection),
                                                             GET_PROP_VALUE(TypeTag, EnableMolecularDiffusion),
                                                             GET_PROP_VALUE(TypeTag, EnableEnergyBalance)> {};

// // specialization for the case of pure advection
// TODO ALL THESE SHOULDNOT BE NECESSARY AND ALWAYS DERIVE FROM ALL CACHES!
template<class TypeTag>
class CCTpfaPorousMediumFluxVariablesCache<TypeTag, true, false, false> : public GET_PROP_TYPE(TypeTag, AdvectionType)::Cache,
                                                                          public GET_PROP_TYPE(TypeTag, MolecularDiffusionType)::Cache {};

// specialization for the case of advection & diffusion
template<class TypeTag>
class CCTpfaPorousMediumFluxVariablesCache<TypeTag, true, true, false> : public GET_PROP_TYPE(TypeTag, AdvectionType)::Cache,
                                                                         public GET_PROP_TYPE(TypeTag, MolecularDiffusionType)::Cache {};

// specialization for the case of advection & heat conduction
template<class TypeTag>
class CCTpfaPorousMediumFluxVariablesCache<TypeTag, true, false, true> : public GET_PROP_TYPE(TypeTag, AdvectionType)::Cache,
                                                                         public GET_PROP_TYPE(TypeTag, MolecularDiffusionType)::Cache,
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

//! specialization of the flux variables cache for the cell centered finite volume mpfa scheme
//! stores data which is commonly used by all the different types of processes
template<class TypeTag>
class PorousMediumFluxVariablesCacheImplementation<TypeTag, DiscretizationMethods::CCMpfa>
      : public CCMpfaPorousMediumFluxVariablesCache<TypeTag, GET_PROP_VALUE(TypeTag, EnableAdvection),
                                                             GET_PROP_VALUE(TypeTag, EnableMolecularDiffusion),
                                                             GET_PROP_VALUE(TypeTag, EnableEnergyBalance)>
{
    using IndexType = typename GET_PROP_TYPE(TypeTag, GridView)::IndexSet::IndexType;
    using ParentType = CCMpfaPorousMediumFluxVariablesCache<TypeTag, GET_PROP_VALUE(TypeTag, EnableAdvection),
                                                                     GET_PROP_VALUE(TypeTag, EnableMolecularDiffusion),
                                                                     GET_PROP_VALUE(TypeTag, EnableEnergyBalance)>;
public:
    //! the constructor
    PorousMediumFluxVariablesCacheImplementation()
    : ParentType(),
      isUpdated_(false)
    {}

    //! Returns whether or not this cache has been updated
    bool isUpdated() const { return isUpdated_; }

    //! Sets the update status. When set to true, consecutive updates will be skipped
    void setUpdateStatus(bool status) { isUpdated_ = status; }

    //! Sets the index of the iv (this scvf is embedded in) in its container
    void setIvIndexInContainer(IndexType ivIndex) { ivIndexInContainer_ = ivIndex; }

    //! Returns the index of the iv (this scvf is embedded in) in its container
    IndexType ivIndexInContainer() const { return ivIndexInContainer_; }

private:
    //! indicates if cache has been fully updated
    bool isUpdated_;

    //! the index of the iv (this scvf is embedded in) in its container
    IndexType ivIndexInContainer_;
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

// TODO remaining specializations

} // end namespace

#endif
