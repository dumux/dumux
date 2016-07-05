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
 * \brief This file contains the data which is required to calculate
 *        diffusive mass fluxes due to molecular diffusion with Fick's law.
 */
#ifndef DUMUX_POROUSMEDIUMFLOW_FICKS_LAW_HH
#define DUMUX_POROUSMEDIUMFLOW_FICKS_LAW_HH

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>

#include <dumux/implicit/properties.hh>


namespace Dumux
{

namespace Properties
{
// forward declaration of properties
NEW_PROP_TAG(NumPhases);
NEW_PROP_TAG(FluidSystem);
NEW_PROP_TAG(EffectiveDiffusivityModel);
}

/*!
 * \ingroup CCTpfaFicksLaw
 * \brief Evaluates the diffusive mass flux according to Fick's law
 */
template <class TypeTag, typename DiscretizationMethod = void>
class FicksLaw
{};

// Specialization for the CC-Tpfa Method
template <class TypeTag>
class FicksLaw<TypeTag, typename std::enable_if<GET_PROP_VALUE(TypeTag, DiscretizationMethod) == GET_PROP(TypeTag, DiscretizationMethods)::CCTpfa>::type >
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, EffectiveDiffusivityModel) EffDiffModel;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace) SubControlVolumeFace;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::IndexSet::IndexType IndexType;
    typedef typename std::vector<IndexType> Stencil;

    using Element = typename GridView::template Codim<0>::Entity;

    enum { dim = GridView::dimension} ;
    enum { dimWorld = GridView::dimensionworld} ;
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases)} ;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimWorldMatrix;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:

    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const SubControlVolumeFace& scvFace,
                       const int phaseIdx,
                       const int compIdx)
    {
        // diffusion tensors are always solution dependent
        Scalar tij = calculateTransmissibility_(problem, scvFace, phaseIdx, compIdx);

        // Get the inside volume variables
        const auto insideScvIdx = scvFace.insideScvIdx();
        const auto& insideScv = problem.model().fvGeometries().subControlVolume(insideScvIdx);
        const auto& insideVolVars = problem.model().curVolVars(insideScv);

        // and the outside volume variables
        const auto& outsideVolVars = problem.model().curVolVars(scvFace.outsideScvIdx());

        // compute the diffusive flux
        const auto xInside = insideVolVars.moleFraction(phaseIdx, compIdx);
        const auto xOutside = outsideVolVars.moleFraction(phaseIdx, compIdx);
        const auto rho = 0.5*(insideVolVars.molarDensity(phaseIdx) + outsideVolVars.molarDensity(phaseIdx));

        return rho*tij*(xInside - xOutside);
    }

    static Stencil stencil(const Problem& problem, const SubControlVolumeFace& scvFace)
    {
        std::vector<IndexType> stencil;
        stencil.clear();
        if (!scvFace.boundary())
        {
            stencil.push_back(scvFace.insideScvIdx());
            stencil.push_back(scvFace.outsideScvIdx());
        }
        else
            stencil.push_back(scvFace.insideScvIdx());

        return stencil;
    }

private:


    static Scalar calculateTransmissibility_(const Problem& problem, const SubControlVolumeFace& scvFace, const int phaseIdx, const int compIdx)
    {
        Scalar tij;

        const auto insideScvIdx = scvFace.insideScvIdx();
        const auto& insideScv = problem.model().fvGeometries().subControlVolume(insideScvIdx);
        const auto& insideVolVars = problem.model().curVolVars(insideScvIdx);

        auto insideD = insideVolVars.diffusionCoefficient(phaseIdx, compIdx);
        insideD = EffDiffModel::effectiveDiffusivity(insideVolVars.porosity(), insideVolVars.saturation(phaseIdx), insideD);
        Scalar ti = calculateOmega_(problem, scvFace, insideD, insideScv);

        if (!scvFace.boundary())
        {
            const auto outsideScvIdx = scvFace.outsideScvIdx();
            const auto& outsideScv = problem.model().fvGeometries().subControlVolume(outsideScvIdx);
            const auto& outsideVolVars = problem.model().curVolVars(outsideScvIdx);

            auto outsideD = outsideVolVars.diffusionCoefficient(phaseIdx, compIdx);
            outsideD = EffDiffModel::effectiveDiffusivity(outsideVolVars.porosity(), outsideVolVars.saturation(phaseIdx), outsideD);
            Scalar tj = -1.0*calculateOmega_(problem, scvFace, outsideD, outsideScv);

            tij = scvFace.area()*(ti * tj)/(ti + tj);
        }
        else
        {
            tij = scvFace.area()*ti;
        }

        return tij;
    }

    static Scalar calculateOmega_(const Problem& problem, const SubControlVolumeFace& scvFace, const DimWorldMatrix &D, const SubControlVolume &scv)
    {
        GlobalPosition Dnormal;
        D.mv(scvFace.unitOuterNormal(), Dnormal);

        auto distanceVector = scvFace.center();
        distanceVector -= scv.center();
        distanceVector /= distanceVector.two_norm2();

        Scalar omega = Dnormal * distanceVector;
        omega *= problem.model().curVolVars(scv).extrusionFactor();

        return omega;
    }

    static Scalar calculateOmega_(const Problem& problem, const SubControlVolumeFace& scvFace, Scalar D, const SubControlVolume &scv)
    {
        auto distanceVector = scvFace.center();
        distanceVector -= scv.center();
        distanceVector /= distanceVector.two_norm2();

        Scalar omega = D * (distanceVector * scvFace.unitOuterNormal());
        omega *= problem.model().curVolVars(scv).extrusionFactor();

        return omega;
    }
};

// Specialization for the Box Method
template <class TypeTag>
class FicksLaw<TypeTag, typename std::enable_if<GET_PROP_VALUE(TypeTag, DiscretizationMethod) == GET_PROP(TypeTag, DiscretizationMethods)::Box>::type >
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using EffDiffModel = typename GET_PROP_TYPE(TypeTag, EffectiveDiffusivityModel);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;

    enum { dim = GridView::dimension} ;
    enum { dimWorld = GridView::dimensionworld} ;
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases)} ;

    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:

    void update(const Problem &problem,
                const Element& element,
                const SubControlVolumeFace &scvFace,
                int phaseIdx, int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "Fick's law for the Box method is not yet implemented!");

        problemPtr_ = &problem;
        scvFacePtr_ = &scvFace;

        phaseIdx_ = phaseIdx;
        compIdx_ = compIdx;

        updateStencil_();

        // TODO for non solution dependent diffusion tensors...
    }

    void update(const Problem &problem, const SubControlVolumeFace &scvFace,
                int phaseIdx, int compIdx,
                VolumeVariables* boundaryVolVars)
    {
        boundaryVolVars_ = boundaryVolVars;
        update(problem, scvFace, phaseIdx, compIdx);
    }

    void beginFluxComputation(bool boundaryVolVarsUpdated = false)
    {
        // diffusion tensors are always solution dependent
        updateTransmissibilities_();

        // Get the inside volume variables
        const auto insideScvIdx = scvFace_().insideScvIdx();
        const auto& insideScv = problem_().model().fvGeometries().subControlVolume(insideScvIdx);
        const auto* insideVolVars = &problem_().model().curVolVars(insideScv);

        // and the outside volume variables
        const VolumeVariables* outsideVolVars;
        if (!scvFace_().boundary())
            outsideVolVars = &problem_().model().curVolVars(scvFace_().outsideScvIdx());
        else
        {
            outsideVolVars = boundaryVolVars_;
            if (!boundaryVolVarsUpdated)
            {
                // update the boudary volvars for Dirichlet boundaries
                const auto element = problem_().model().fvGeometries().element(insideScv);
                const auto dirichletPriVars = problem_().dirichlet(element, scvFace_());
                boundaryVolVars_->update(dirichletPriVars, problem_(), element, insideScv);
            }
        }

        // compute the diffusive flux
        const auto xInside = insideVolVars->moleFraction(phaseIdx_, compIdx_);
        const auto xOutside = outsideVolVars->moleFraction(phaseIdx_, compIdx_);
        const auto rho = 0.5*(insideVolVars->molarDensity(phaseIdx_) + outsideVolVars->molarDensity(phaseIdx_));

        rhoDGradXNormal_ = rho*tij_*(xInside - xOutside);
    }



    /*!
     * \brief A function to calculate the mass flux over a sub control volume face
     *
     * \param phaseIdx The index of the phase of which the flux is to be calculated
     * \param compIdx The index of the transported component
     */
    Scalar flux() const
    {
        return rhoDGradXNormal_;
    }

    std::set<IndexType> stencil() const
    {
        return stencil_;
    }

protected:


    void updateTransmissibilities_()
    {
        const auto insideScvIdx = scvFace_().insideScvIdx();
        const auto& insideScv = problem_().model().fvGeometries().subControlVolume(insideScvIdx);
        const auto& insideVolVars = problem_().model().curVolVars(insideScvIdx);

        auto insideD = insideVolVars.diffusionCoefficient(phaseIdx_, compIdx_);
        insideD = EffDiffModel::effectiveDiffusivity(insideVolVars.porosity(), insideVolVars.saturation(phaseIdx_), insideD);
        Scalar ti = calculateOmega_(insideD, insideScv);

        if (!scvFace_().boundary())
        {
            const auto outsideScvIdx = scvFace_().outsideScvIdx();
            const auto& outsideScv = problem_().model().fvGeometries().subControlVolume(outsideScvIdx);
            const auto& outsideVolVars = problem_().model().curVolVars(outsideScvIdx);

            auto outsideD = outsideVolVars.diffusionCoefficient(phaseIdx_, compIdx_);
            outsideD = EffDiffModel::effectiveDiffusivity(outsideVolVars.porosity(), outsideVolVars.saturation(phaseIdx_), outsideD);
            Scalar tj = -1.0*calculateOmega_(outsideD, outsideScv);

            tij_ = scvFace_().area()*(ti * tj)/(ti + tj);
        }
        else
        {
            tij_ = scvFace_().area()*ti;
        }
    }

    void updateStencil_()
    {
        // fill the stencil
        if (!scvFace_().boundary())
            stencil_= {scvFace_().insideScvIdx(), scvFace_().outsideScvIdx()};
        else
            // fill the stencil
            stencil_ = {scvFace_().insideScvIdx()};
    }

    Scalar calculateOmega_(const DimWorldMatrix &D, const SubControlVolume &scv) const
    {
        GlobalPosition Dnormal;
        D.mv(scvFace_().unitOuterNormal(), Dnormal);

        auto distanceVector = scvFace_().center();
        distanceVector -= scv.center();
        distanceVector /= distanceVector.two_norm2();

        Scalar omega = Dnormal * distanceVector;
        omega *= problem_().model().curVolVars(scv).extrusionFactor();

        return omega;
    }

    Scalar calculateOmega_(Scalar D, const SubControlVolume &scv) const
    {
        auto distanceVector = scvFace_().center();
        distanceVector -= scv.center();
        distanceVector /= distanceVector.two_norm2();

        Scalar omega = D * (distanceVector * scvFace_().unitOuterNormal());
        omega *= problem_().model().curVolVars(scv).extrusionFactor();

        return omega;
    }

    const Problem &problem_() const
    {
        return *problemPtr_;
    }

    const SubControlVolumeFace& scvFace_() const
    {
        return *scvFacePtr_;
    }

    const Problem *problemPtr_;
    const SubControlVolumeFace *scvFacePtr_; //!< Pointer to the sub control volume face for which the flux variables are created
    std::set<IndexType> stencil_;         //!< Indices of the cells of which the pressure is needed for the flux calculation

    //! Boundary volume variables (they only get updated on Dirichlet boundaries)
    VolumeVariables* boundaryVolVars_;

    IndexType phaseIdx_;
    IndexType compIdx_;
    //! Precomputed values
    Scalar tij_; //!< transmissibility for the flux calculation
    Scalar rhoDGradXNormal_; //! rho*D(grad(x))*n
};

} // end namespace

#endif
