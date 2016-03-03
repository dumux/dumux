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
#ifndef DUMUX_CC_TPFA_DARCY_FLUX_VARIABLES_HH
#define DUMUX_CC_TPFA_DARCY_FLUX_VARIABLES_HH

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
 * \ingroup ImplicitFluxVariables
 * \brief Evaluates the normal component of the Darcy velocity
 * on a (sub)control volume face.
 */
template <class TypeTag>
class CCTpfaImplicitDarcyFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace) SubControlVolumeFace;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::IndexSet::IndexType IndexType;

    enum { dim = GridView::dimension} ;
    enum { dimWorld = GridView::dimensionworld} ;
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases)} ;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimWorldMatrix;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:

    void update(const Problem &problem, const SubControlVolumeFace &scvFace)
    {
        problemPtr_ = &problem;
        scvFacePtr_ = &scvFace;

        updateTransmissibilities_();
        updateStencil_();
    }

    void updateTransmissibilities(const Problem &problem, const SubControlVolumeFace &scvFace)
    {
        updateTransmissibilities_();
    }

    /*!
     * \brief A function to calculate the mass flux over a sub control volume face
     *
     * \param phaseIdx The index of the phase of which the flux is to be calculated
     * \param upwindFunction A function which does the upwinding
     */
    template<typename FunctionType>
    Scalar computeFlux(IndexType phaseIdx, FunctionType upwindFunction) const
    {
        const auto insideScvIdx = scvFace_().insideScvIdx();
        const auto& insideScv = problem_().model().fvGeometries().subControlVolume(insideScvIdx);

        // Set the inside-/outside volume variables provisionally
        const auto &insideVolVars = problem_().model().curVolVars(insideScv);
        VolumeVariables tmp;
        VolumeVariables& outsideVolVars = tmp;

        // calculate the piezometric heights in the two cells
        // and set the downstream volume variables
        Scalar hInside = insideVolVars.pressure(phaseIdx);
        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
        {
            // ask for the gravitational acceleration in the inside cell
            GlobalPosition gInside(problem_().gravityAtPos(insideScv.center()));
            Scalar rhoInside = insideVolVars.density(phaseIdx);

            hInside -= rhoInside*(gInside*insideScv.center());
        }

        Scalar hOutside;
        if (!scvFace_().boundary())
        {
            outsideVolVars = problem_().model().curVolVars(scvFace_().outsideScvIdx());
            hOutside = outsideVolVars.pressure(phaseIdx);

            // if switched on, ask for the gravitational acceleration in the outside cell
            if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
            {
                const auto outsideScvIdx = scvFace_().outsideScvIdx();
                const auto& outsideScv = problem_().model().fvGeometries().subControlVolume(outsideScvIdx);
                GlobalPosition gOutside(problem_().gravityAtPos(outsideScv.center()));
                Scalar rhoOutside = outsideVolVars.density(phaseIdx);

                hOutside -= rhoOutside*(gOutside*outsideScv.center());
            }
        }
        else
        {
            // if (complexBCTreatment)
            // outsideVolVars = problem_().model().constructBoundaryVolumeVariables(scvFace_());
            // else
            // {
                auto element = problem_().model().fvGeometries().element(insideScv);
                auto dirichletPriVars = problem_().dirichlet(element, scvFace_());
                outsideVolVars.update(dirichletPriVars, problem_(), element, insideScv);
            // }
            hOutside = outsideVolVars.pressure(phaseIdx);

            // if switched on, ask for the gravitational acceleration in the outside cell
            if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
            {
                GlobalPosition gOutside(problem_().gravityAtPos(scvFace_().center()));
                Scalar rhoOutside = outsideVolVars.density(phaseIdx);

                hOutside -= rhoOutside*(gOutside*scvFace_().center());
            }
        }

        auto volumeFlux = tij_*(hInside - hOutside);

        if (std::signbit(volumeFlux))
            return volumeFlux*upwindFunction(outsideVolVars, insideVolVars);
        else
            return volumeFlux*upwindFunction(insideVolVars, outsideVolVars);
    }

    std::set<IndexType> stencil() const
    {
        return std::set<IndexType>(stencil_.begin(), stencil_.end());
    }

protected:


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
            stencil_= {scvFace_().insideVolVarsIdx(), scvFace_().outsideVolVarsIdx()};
        else
            // fill the stencil
            stencil_ = {scvFace_().insideVolVarsIdx()};
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

    const Problem *problemPtr_;
    const SubControlVolumeFace *scvFacePtr_; //!< Pointer to the sub control volume face for which the flux variables are created
    Scalar tij_;                             //!< transmissibility for the flux calculation
    std::vector<IndexType> stencil_;         //!< Indices of the cells of which the pressure is needed for the flux calculation
};

} // end namespace

#endif // DUMUX_CC_TPFA_DARCY_FLUX_VARIABLES_HH
