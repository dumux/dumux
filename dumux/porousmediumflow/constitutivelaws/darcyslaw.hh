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
 *        volume and mass fluxes of fluid phases over a face of a finite volume by means
 *        of the Darcy approximation. Specializations are provided for the different discretization methods.
 */
#ifndef DUMUX_POROUSMEDIUMFLOW_DARCYS_LAW_HH
#define DUMUX_POROUSMEDIUMFLOW_DARCYS_LAW_HH

#include <memory>

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
NEW_PROP_TAG(ProblemEnableGravity);
}

/*!
 * \ingroup DarcysLaw
 * \brief Evaluates the normal component of the Darcy velocity
 * on a (sub)control volume face. Specializations are provided
 * for the different discretization methods.
 */
template <class TypeTag, typename DiscretizationMethod = void>
class DarcysLaw
{};

// Specialization for the CC-Tpfa method
template <class TypeTag>
class DarcysLaw<TypeTag, typename std::enable_if<GET_PROP_VALUE(TypeTag, DiscretizationMethod) == GET_PROP(TypeTag, DiscretizationMethods)::CCTpfa>::type >
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace) SubControlVolumeFace;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GridView::IndexSet::IndexType IndexType;
    typedef std::vector<IndexType> Stencil;

    using Element = typename GridView::template Codim<0>::Entity;
    enum { dim = GridView::dimension} ;
    enum { dimWorld = GridView::dimensionworld} ;
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases)} ;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimWorldMatrix;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:

    static Scalar flux(const Problem& problem,
                       const SubControlVolumeFace& scvFace,
                       const IndexType phaseIdx)
    {
        const auto& tij = getTransmissibilities(problem, scvFace);

        // Get the inside volume variables
        const auto insideScvIdx = scvFace.insideScvIdx();
        const auto& insideScv = problem.model().fvGeometries().subControlVolume(insideScvIdx);
        const auto& insideVolVars = problem.model().curVolVars(insideScv);

        // and the outside volume variables
        const auto& outsideVolVars = problem.model().curVolVars(scvFace.outsideScvIdx());

        auto hInside = insideVolVars.pressure(phaseIdx);
        auto hOutside = outsideVolVars.pressure(phaseIdx);

        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
        {
            // do averaging for the density
            const auto rhoInside = insideVolVars.density(phaseIdx);
            const auto rhoOutide = outsideVolVars.density(phaseIdx);
            const auto rho = (rhoInside + rhoOutide)*0.5;


            // ask for the gravitational acceleration in the inside neighbor
            const auto xInside = insideScv.center();
            const auto gInside = problem.gravityAtPos(xInside);

            hInside -= rho*(gInside*xInside);

            // and the outside neighbor
            if (scvFace.boundary())
            {
                const auto xOutside = scvFace.center();
                const auto gOutside = problem.gravityAtPos(xOutside);
                hOutside -= rho*(gOutside*xOutside);
            }
            else
            {
                const auto outsideScvIdx = scvFace.outsideScvIdx();
                const auto& outsideScv = problem.model().fvGeometries().subControlVolume(outsideScvIdx);
                const auto xOutside = outsideScv.center();
                const auto gOutside = problem.gravityAtPos(xOutside);
                hOutside -= rho*(gOutside*xOutside);
            }
        }

        return tij*(hInside - hOutside);
    }

    static Stencil stencil(const Problem& problem, const SubControlVolumeFace& scvFace)
    {
        std::vector<IndexType> stencil;
        if (!scvFace.boundary())
        {
            stencil.push_back(scvFace.insideScvIdx());
            stencil.push_back(scvFace.outsideScvIdx());
        }
        else
            stencil.push_back(scvFace.insideScvIdx());

        return stencil;
    }

    template <typename T = TypeTag>
    static const typename std::enable_if<GET_PROP_VALUE(T, EnableFluxVariablesCache), Scalar>::type& getTransmissibilities(const Problem& problem, const SubControlVolumeFace& scvFace)
    {
        return problem.model().fluxVarsCache(scvFace).tij();
    }

    template <typename T = TypeTag>
    static const typename std::enable_if<!GET_PROP_VALUE(T, EnableFluxVariablesCache), Scalar>::type getTransmissibilities(const Problem& problem, const SubControlVolumeFace& scvFace)
    {
        return calculateTransmissibilities(problem, scvFace);
    }

    static Scalar calculateTransmissibilities(const Problem& problem, const SubControlVolumeFace& scvFace)
    {
        Scalar tij;

        const auto insideScvIdx = scvFace.insideScvIdx();
        const auto& insideScv = problem.model().fvGeometries().subControlVolume(insideScvIdx);
        const auto insideK = problem.spatialParams().intrinsicPermeability(insideScv);
        Scalar ti = calculateOmega_(problem, scvFace, insideK, insideScv);

        if (!scvFace.boundary())
        {
            const auto outsideScvIdx = scvFace.outsideScvIdx();
            const auto& outsideScv = problem.model().fvGeometries().subControlVolume(outsideScvIdx);
            const auto outsideK = problem.spatialParams().intrinsicPermeability(outsideScv);
            Scalar tj = -1.0*calculateOmega_(problem, scvFace, outsideK, outsideScv);

            tij = scvFace.area()*(ti * tj)/(ti + tj);
        }
        else
        {
            tij = scvFace.area()*ti;
        }

        return tij;
    }

private:

    static Scalar calculateOmega_(const Problem& problem, const SubControlVolumeFace& scvFace, const DimWorldMatrix &K, const SubControlVolume &scv)
    {
        GlobalPosition Knormal;
        K.mv(scvFace.unitOuterNormal(), Knormal);

        auto distanceVector = scvFace.center();
        distanceVector -= scv.center();
        distanceVector /= distanceVector.two_norm2();

        Scalar omega = Knormal * distanceVector;
        omega *= problem.model().curVolVars(scv).extrusionFactor();

        return omega;
    }

    static Scalar calculateOmega_(const Problem& problem, const SubControlVolumeFace& scvFace, Scalar K, const SubControlVolume &scv)
    {
        auto distanceVector = scvFace.center();
        distanceVector -= scv.center();
        distanceVector /= distanceVector.two_norm2();

        Scalar omega = K * (distanceVector * scvFace.unitOuterNormal());
        omega *= problem.model().curVolVars(scv).extrusionFactor();

        return omega;
    }
};

// Specialization for the Box Method
template <class TypeTag>
class DarcysLaw<TypeTag, typename std::enable_if<GET_PROP_VALUE(TypeTag, DiscretizationMethod) == GET_PROP(TypeTag, DiscretizationMethods)::Box>::type >
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
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

    void update(const Problem& problem,
                const Element& element,
                const SubControlVolumeFace& scvFace)
    {
        DUNE_THROW(Dune::NotImplemented, "Darcy's law for the Box method is not yet implemented!");

        problemPtr_ = &problem;
        scvFacePtr_ = &scvFace;
        elementPtr_ = &element;
        enableGravity_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity);

        updateTransmissibilities_();
    }

    void update(const Problem& problem,
                const Element& element,
                const SubControlVolumeFace &scvFace,
                VolumeVariables* boundaryVolVars)
    {
        update(problem, element, scvFace);
    }

    void updateTransmissibilities(const Problem &problem, const SubControlVolumeFace &scvFace)
    {
        updateTransmissibilities_();
    }

    void beginFluxComputation(bool boundaryVolVarsUpdated = false)
    {
        // Get the inside volume variables
        const auto insideScvIdx = scvFace_().insideScvIdx();
        const auto& insideScv = problem_().model().fvGeometries().subControlVolume(insideScvIdx);
        const auto* insideVolVars = &problem_().model().curVolVars(insideScv);

        // and the outside volume variables
        const VolumeVariables* outsideVolVars;
        outsideVolVars = &problem_().model().curVolVars(scvFace_().outsideScvIdx());

        // loop over all phases to compute the volume flux
        for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
        {
            auto hInside = insideVolVars->pressure(phaseIdx);
            auto hOutside = outsideVolVars->pressure(phaseIdx);

            if (enableGravity_)
            {
                // do averaging for the density
                const auto rhoInside = insideVolVars->density(phaseIdx);
                const auto rhoOutide = outsideVolVars->density(phaseIdx);
                const auto rho = (rhoInside + rhoOutide)*0.5;


                // ask for the gravitational acceleration in the inside neighbor
                const auto xInside = insideScv.center();
                const auto gInside = problem_().gravityAtPos(xInside);

                hInside -= rho*(gInside*xInside);

                const auto outsideScvIdx = scvFace_().outsideScvIdx();
                const auto& outsideScv = problem_().model().fvGeometries().subControlVolume(outsideScvIdx);
                const auto xOutside = outsideScv.center();
                const auto gOutside = problem_().gravityAtPos(xOutside);
                hOutside -= rho*(gOutside*xOutside);
            }

            kGradPNormal_[phaseIdx] = tij_*(hInside - hOutside);

            if (std::signbit(kGradPNormal_[phaseIdx]))
            {
                upWindIndices_[phaseIdx] = std::make_pair(scvFace_().outsideScvIdx(), scvFace_().insideScvIdx());
            }
            else
            {
                upWindIndices_[phaseIdx] = std::make_pair(scvFace_().insideScvIdx(), scvFace_().outsideScvIdx());
            }
        }
    }

    /*!
     * \brief A function to calculate the mass flux over a sub control volume face
     *
     * \param phaseIdx The index of the phase of which the flux is to be calculated
     * \param upwindFunction A function which does the upwinding
     */
    template<typename FunctionType>
    Scalar flux(IndexType phaseIdx, FunctionType upwindFunction) const
    {
        return kGradPNormal_[phaseIdx]*upwindFunction(upVolVars(phaseIdx), dnVolVars(phaseIdx));
    }

    // for compatibility with cell-centered models
    const std::set<IndexType>& stencil() const
    {
        return std::set<IndexType>();
    }

    const VolumeVariables& upVolVars(IndexType phaseIdx) const
    {
        return problem_().model().curVolVars(upWindIndices_[phaseIdx].first);
    }

    const VolumeVariables& dnVolVars(IndexType phaseIdx) const
    {
        return problem_().model().curVolVars(upWindIndices_[phaseIdx].second);
    }

private:

    void updateTransmissibilities_()
    {
        const auto insideScvIdx = scvFace_().insideScvIdx();
        const auto& insideScv = problem_().model().fvGeometries().subControlVolume(insideScvIdx);
        const auto insideK = problem_().spatialParams().intrinsicPermeability(insideScv);
        Scalar ti = calculateOmega_(insideK, insideScv);

        const auto outsideScvIdx = scvFace_().outsideScvIdx();
        const auto& outsideScv = problem_().model().fvGeometries().subControlVolume(outsideScvIdx);
        const auto outsideK = problem_().spatialParams().intrinsicPermeability(outsideScv);
        Scalar tj = -1.0*calculateOmega_(outsideK, outsideScv);

        tij_ = scvFace_().area()*(ti * tj)/(ti + tj);
    }

    const Problem &problem_() const
    {
        return *problemPtr_;
    }

    const SubControlVolumeFace& scvFace_() const
    {
        return *scvFacePtr_;
    }

    const Element& element_() const
    {
        return *elementPtr_;
    }

    const Problem *problemPtr_; //! Pointer to the problem
    const SubControlVolumeFace *scvFacePtr_; //! Pointer to the sub control volume face for which the flux variables are created
    const Element *elementPtr_; //! Point to the element
    bool enableGravity_; //! If we have a problem considering gravitational effects

    //! The upstream (first) and downstream (second) volume variable indices
    std::array<std::pair<IndexType, IndexType>, numPhases> upWindIndices_;
    Scalar tij_ = 0;

    //! Precomputed values
    std::array<Scalar, numPhases> kGradPNormal_; //! K(grad(p) - rho*g)*n
    GlobalPosition normalK_; //! (K^T)n
};

} // end namespace

#endif
