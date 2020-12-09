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
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesFluxVariablesImpl
 */
#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_DIAMOND_FLUXVARIABLES_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_DIAMOND_FLUXVARIABLES_HH

#include <array>
#include <optional>

#include <dumux/common/math.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief The flux variables class for the Navier-Stokes model using the staggered grid discretization.
 */
template<class TypeTag>
class NavierStokesMomentumFluxVariables
{
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    using GridFluxVariablesCache = typename GridVariables::GridFluxVariablesCache;
    using ElementFluxVariablesCache = typename GridFluxVariablesCache::LocalView;
    using FluxVariablesCache = typename GridFluxVariablesCache::FluxVariablesCache;

    using GridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    using VelocityGradients = StaggeredVelocityGradients;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Extrusion = Extrusion_t<GridGeometry>;

    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using FeLocalBasis = typename GridGeometry::FeCache::FiniteElementType::Traits::LocalBasisType;
    using ShapeJacobian = typename FeLocalBasis::Traits::JacobianType;
    using ShapeValue = typename Dune::FieldVector<Scalar, 1>;

    using Tensor = Dune::FieldMatrix<Scalar, dim, dimWorld>;
    static_assert(NumEqVector::dimension == dimWorld, "Wrong dimension of velocity vector");

public:

    NavierStokesMomentumFluxVariables(const Problem& problem,
                                      const Element& element,
                                      const FVElementGeometry& fvGeometry,
                                      const SubControlVolumeFace& scvFace,
                                      const ElementVolumeVariables& elemVolVars,
                                      const ElementFluxVariablesCache& elemFluxVarsCache,
                                      const ElementBoundaryTypes& elemBcTypes)
    : problemPtr_(&problem)
    , elementPtr_(&element)
    , fvGeometryPtr_(&fvGeometry)
    , scvFacePtr_(&scvFace)
    , elemVolVarsPtr_(&elemVolVars)
    , elemFluxVarsCachePtr_(&elemFluxVarsCache)
    , elemBcTypesPtr_(&elemBcTypes)
    {}

    const Problem& problem() const
    { return *problemPtr_; }

    const Element& element() const
    { return *elementPtr_; }

    const SubControlVolumeFace& scvFace() const
    { return *scvFacePtr_; }

    const FVElementGeometry& fvGeometry() const
    { return *fvGeometryPtr_; }

    const ElementVolumeVariables& elemVolVars() const
    { return *elemVolVarsPtr_; }

    const ElementFluxVariablesCache& elemFluxVarsCache() const
    { return *elemFluxVarsCachePtr_; }

    const ElementBoundaryTypes& elemBcTypes() const
    { return *elemBcTypesPtr_; }

    /*!
     * \brief Returns the diffusive momentum flux due to viscous forces
     */
    NumEqVector advectiveMomentumFlux() const
    {
        if (!this->problem().enableInertiaTerms())
            return NumEqVector(0.0);

        const auto& fvGeometry = this->fvGeometry();
        const auto& elemVolVars = this->elemVolVars();
        const auto& geometry = this->element().geometry();
        const auto& scvf = this->scvFace();
        const auto ipLocal = geometry.local(scvf.ipGlobal());

        // TODO: Cache the shapeValues
        const auto& localBasis = fvGeometry.feLocalBasis();
        std::vector<ShapeValue> shapeValues;
        localBasis.evaluateFunction(ipLocal, shapeValues);

        // evaluate gradP - rho*g at integration point
        NumEqVector v(0.0);
        for (auto&& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];

            // interpolate velocity at scvf
            v.axpy(shapeValues[scv.indexInElement()][0], volVars.velocity());
        }

        //TODO interpolate rho
        const Scalar density = this->problem().density(this->element(), this->fvGeometry(), scvf);

        Scalar upwindFactor = 0.0;
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        if(v*scvf.unitOuterNormal() >= 0)
            upwindFactor = density*(insideVolVars.velocity()*scvf.unitOuterNormal());
        else
            upwindFactor = density*(outsideVolVars.velocity()*scvf.unitOuterNormal());

        return upwindFactor*v;
    }

    /*!
     * \brief Returns the diffusive momentum flux due to viscous forces
     */
    NumEqVector diffusiveMomentumFlux() const
    {
        NumEqVector result(0.0);

        //TODO Boundary handling
        const auto& fvGeometry = this->fvGeometry();
        const auto& elemVolVars = this->elemVolVars();
        const auto& geometry = this->element().geometry();
        const auto& scvf = this->scvFace();
        const auto ipLocal = geometry.local(scvf.ipGlobal());

        const auto& localBasis = fvGeometry.feLocalBasis();

        // TODO cache these values
        const auto jacInvT = geometry.jacobianInverseTransposed(ipLocal);
        std::vector<ShapeJacobian> shapeJacobian;
        localBasis.evaluateJacobian(ipLocal, shapeJacobian);

        // compute the gradN at for every scv/dof
        std::vector<GlobalPosition> gradN(fvGeometry.numScv(),0.0);
        for (const auto& scv: scvs(fvGeometry))
            jacInvT.mv(shapeJacobian[scv.localDofIndex()][0], gradN[scv.indexInElement()]);

        Tensor gradV(0.0);
        for (int dir = 0; dir < dim; ++dir)
            for (const auto& scv : scvs(fvGeometry))
                gradV[dir].axpy(elemVolVars[scv.indexInElement()].velocity(dir), gradN[scv.indexInElement()]);

        static const bool enableUnsymmetrizedVelocityGradient
            = getParamFromGroup<bool>(this->problem().paramGroup(), "FreeFlow.EnableUnsymmetrizedVelocityGradient", false);

        result = enableUnsymmetrizedVelocityGradient ? gradV*scvf.unitOuterNormal() : (gradV + getTransposed(gradV))*scvf.unitOuterNormal();

        const auto mu = this->problem().effectiveViscosity(this->element(), this->fvGeometry(), this->scvFace());
        result *= -mu * Extrusion::area(scvf) * extrusionFactor_(elemVolVars, scvf);

        static const bool enableDilatationTerm = getParamFromGroup<bool>(this->problem().paramGroup(), "FreeFlow.EnableDilatationTerm", false);
        if (enableDilatationTerm)
        {
            Scalar divergence = 0.0;
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto frontalScvf  = *(scvfs(fvGeometry, scv).begin());
                assert(frontalScvf.isFrontal() && !frontalScvf.boundary());
                divergence += VelocityGradients::velocityGradII(fvGeometry, frontalScvf, elemVolVars);
            }
            result += 2.0/3.0 * mu * trace(gradV) * scvf.unitOuterNormal() * Extrusion::area(scvf) * extrusionFactor_(elemVolVars, scvf);;
        }

        return result;
    }

    NumEqVector pressureContribution() const
    {
        const auto& scvf = this->scvFace();
        NumEqVector result(scvf.unitOuterNormal());

        // The pressure force needs to take the extruded scvf area into account.
        const auto pressure = this->problem().pressure(this->element(), this->fvGeometry(), scvf);

        // The pressure contribution calculated above might have a much larger numerical value compared to the viscous or inertial forces.
        // This may lead to numerical inaccuracies due to loss of significance (cancellantion) for the final residual value.
        // In the end, we are only interested in a pressure difference between the two relevant faces so we can
        // substract a reference value from the actual pressure contribution. Assuming an axisparallel cartesian grid,
        // scvf.area() will have the same value at both opposing faces such that the reference pressure contribution
        // cancels out in the final residual which combines the pressure contribution of two adjacent elements
        // We explicitly do extrude the area here because that might yield different results in both elements.
        // The multiplication by scvf.area() aims at having a reference value of the same order of magnitude as the actual pressure contribution.
        const auto referencePressure = this->problem().referencePressure(this->element(), this->fvGeometry(), scvf);
        // TODO why not using the same area and extrusion for both pressures?
        result *= (pressure*Extrusion::area(scvf)* extrusionFactor_(this->elemVolVars(), scvf);
                   - referencePressure*scvf.area());

        return result;
    }

private:

    template<class ElementVolumeVariables, class SubControlVolumeFace>
    Scalar extrusionFactor_(const ElementVolumeVariables& elemVolVars, const SubControlVolumeFace& scvf) const
    {
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        return harmonicMean(insideVolVars.extrusionFactor(), outsideVolVars.extrusionFactor());
    }


    const Problem* problemPtr_;                             //!< Pointer to the problem
    const Element* elementPtr_;                             //!< Pointer to the element at hand
    const FVElementGeometry* fvGeometryPtr_;                //!< Pointer to the current FVElementGeometry
    const SubControlVolumeFace* scvFacePtr_;                //!< Pointer to the sub control volume face for which the flux variables are created
    const ElementVolumeVariables* elemVolVarsPtr_;          //!< Pointer to the current element volume variables
    const ElementFluxVariablesCache* elemFluxVarsCachePtr_; //!< Pointer to the current element flux variables cache
    const ElementBoundaryTypes* elemBcTypesPtr_; //!< Pointer to element boundary types
};

} // end namespace Dumux

#endif // DUMUX_NAVIERSTOKES_DIAMOND_FLUXVARIABLES_HH
