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
 *        of the Darcy approximation.
 */
#ifndef DUMUX_CC_TPFA_DARCYS_LAW_HH
#define DUMUX_CC_TPFA_DARCYS_LAW_HH

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
 * \ingroup CCTpfaDarcysLaw
 * \brief Evaluates the normal component of the Darcy velocity
 * on a (sub)control volume face.
 */
template <class TypeTag>
class CCTpfaDarcysLaw
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace) SubControlVolumeFace;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::IndexSet::IndexType IndexType;
    using Element = typename GridView::template Codim<0>::Entity;

    enum { dim = GridView::dimension} ;
    enum { dimWorld = GridView::dimensionworld} ;
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases)} ;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimWorldMatrix;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:

    void update(const Problem &problem,
                const Element& element,
                const SubControlVolumeFace &scvFace)
    {
        problemPtr_ = &problem;
        scvFacePtr_ = &scvFace;
        enableGravity_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity);

        updateTransmissibilities_();
        updateStencil_();
    }

    void update(const Problem& problem,
                const Element& element,
                const SubControlVolumeFace &scvFace,
                VolumeVariables* boundaryVolVars)
    {
        boundaryVolVars_ = boundaryVolVars;
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

                // and the outside neighbor
                if (scvFace_().boundary())
                {
                    const auto xOutside = scvFace_().center();
                    const auto gOutside = problem_().gravityAtPos(xOutside);
                    hOutside -= rho*(gOutside*xOutside);
                }
                else
                {
                    const auto outsideScvIdx = scvFace_().outsideScvIdx();
                    const auto& outsideScv = problem_().model().fvGeometries().subControlVolume(outsideScvIdx);
                    const auto xOutside = outsideScv.center();
                    const auto gOutside = problem_().gravityAtPos(xOutside);
                    hOutside -= rho*(gOutside*xOutside);
                }
            }

            kGradPNormal_[phaseIdx] = tij_*(hInside - hOutside);

            if (std::signbit(kGradPNormal_[phaseIdx]))
            {
                if (scvFace_().boundary())
                    upWindIndices_[phaseIdx] = std::make_pair(-1, scvFace_().insideScvIdx());
                else
                    upWindIndices_[phaseIdx] = std::make_pair(scvFace_().outsideScvIdx(), scvFace_().insideScvIdx());
            }
            else
            {
                if (scvFace_().boundary())
                    upWindIndices_[phaseIdx] = std::make_pair(scvFace_().insideScvIdx(), -1);
                else
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

    const std::set<IndexType>& stencil() const
    {
        return stencil_;
    }

    const VolumeVariables& upVolVars(IndexType phaseIdx) const
    {
        if(upWindIndices_[phaseIdx].first != -1)
            return problem_().model().curVolVars(upWindIndices_[phaseIdx].first);
        else
            return *boundaryVolVars_;
    }

    const VolumeVariables& dnVolVars(IndexType phaseIdx) const
    {
        if(upWindIndices_[phaseIdx].second != -1)
            return problem_().model().curVolVars(upWindIndices_[phaseIdx].second);
        else
            return *boundaryVolVars_;
    }

private:

    void updateTransmissibilities_()
    {
        if (!scvFace_().boundary())
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
        else
        {
            const auto insideScvIdx = scvFace_().insideScvIdx();
            const auto& insideScv = problem_().model().fvGeometries().subControlVolume(insideScvIdx);
            const auto insideK = problem_().spatialParams().intrinsicPermeability(insideScv);
            Scalar ti = calculateOmega_(insideK, insideScv);

            tij_ = scvFace_().area()*ti;
        }
    }

    void updateStencil_()
    {
        // fill the stencil
        if (!scvFace_().boundary())
            stencil_.insert({scvFace_().insideScvIdx(), scvFace_().outsideScvIdx()});
        else
            // fill the stencil
            stencil_.insert(scvFace_().insideScvIdx());
    }

    Scalar calculateOmega_(const DimWorldMatrix &K, const SubControlVolume &scv) const
    {
        GlobalPosition Knormal;
        K.mv(scvFace_().unitOuterNormal(), Knormal);

        auto distanceVector = scvFace_().center();
        distanceVector -= scv.center();
        distanceVector /= distanceVector.two_norm2();

        Scalar omega = Knormal * distanceVector;
        omega *= problem_().model().curVolVars(scv).extrusionFactor();

        return omega;
    }

    Scalar calculateOmega_(Scalar K, const SubControlVolume &scv) const
    {
        auto distanceVector = scvFace_().center();
        distanceVector -= scv.center();
        distanceVector /= distanceVector.two_norm2();

        Scalar omega = K * (distanceVector * scvFace_().unitOuterNormal());
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

    const Problem *problemPtr_;              //! Pointer to the problem
    const SubControlVolumeFace *scvFacePtr_; //! Pointer to the sub control volume face for which the flux variables are created
    bool enableGravity_;                     //! If we have a problem considering gravitational effects
    std::set<IndexType> stencil_;         //! Indices of the cells of which the pressure is needed for the flux calculation

    //! Boundary volume variables (they only get updated on Dirichlet boundaries)
    VolumeVariables* boundaryVolVars_;

    //! The upstream (first) and downstream (second) volume variable indices
    std::array<std::pair<IndexType, IndexType>, numPhases> upWindIndices_;

    //! Precomputed values
    Scalar tij_; //! transmissibility for the flux calculation tij(ui - uj)
    std::array<Scalar, numPhases> kGradPNormal_; //! K(grad(p) - rho*g)*n
};

} // end namespace

#endif
