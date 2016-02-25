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

    enum { dim = GridView::dimension} ;
    enum { dimWorld = GridView::dimensionworld} ;
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases)} ;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimWorldMatrix;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:

    void update(const Problem &problem, const SubControlVolumeFace &scvFace)
    {
        problemPtr_ = &problem;
        scvFace_ = &scvFace;

        updateTransmissivitesAndStencil_();
    }

    /*!
     * \brief A function to calculate the mass flux over a sub control volume face
     *
     * \param phaseIdx The index of the phase of which the flux is to be calculated
     * \param upwindFunction A function which does the upwinding
     */
    template<FunctionType>
    Scalar calculateFlux(const int phaseIdx, FunctionType upwindFunction)
    {
        // Set the inside-/outside volume variables provisionally
        const VolumeVariables &insideVolVars = problem_().model().curVolVars(stencil_[0]);
        const VolumeVariables &outsideVolVars = insideVolVars;

        // calculate the piezometric heights in the two cells
        // and set the downstream volume variables
        Scalar hInside = insideVolVars.pressure(phaseIdx);
        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
        {
            // ask for the gravitational acceleration in the inside cell
            const auto insideScvIdx = scvFace_.insideScvIdx();
            const auto& insideScv = problem_().model().fvGeometries().subControlVolume(insideScvIdx);
            GlobalPosition gInside(problem_().gravityAtPos(insideScv.center()));
            Scalar rhoInside = insideVolVars.density(phaseIdx);

            hInside -= rhoInside*gInside*insideScv.center();
        }

        Scalar hOutside;
        if (!scvFace_.boundary())
        {
            outsideVolVars = problem_().model().curVolVars(stencil_[1]);
            hOutside = outsideVolVars.pressure(phaseIdx);

            // if switched on, ask for the gravitational acceleration in the outside cell
            if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
            {
                const auto outsideScvIdx = scvFace_.outsideScvIdx();
                const auto& outsideScv = problem_().model().fvGeometries().subControlVolume(outsideScvIdx);
                GlobalPosition gOutside(problem_().gravityAtPos(outsideScv.center()));
                Scalar rhoOutside = outsideVolVars.density(phaseIdx);

                hOutside -= rhoOutside*gOutside*outsideScv.center();
            }
        }
        else
        {
            // if (complexBCTreatment)
            // outsideVolVars = problem_().model().constructBoundaryVolumeVariables(scvFace_);
            // else
            // {
                PrimaryVariables dirichletPV = problem_().dirichlet(scvFace_);
                const auto insideScvIdx = scvFace_.insideScvIdx();
                const auto& insideScv = problem_().model().fvGeometries().subControlVolume(insideScvIdx);
                outsideVolVars.update(dirichletPV, problem_(), outsideScv);
            // }
            h2 = outsideVolVars.pressure(phaseIdx);

            // if switched on, ask for the gravitational acceleration in the outside cell
            if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
            {
                GlobalPosition gOutside(problem_().gravityAtPos(scvFace_.center()));
                Scalar rhoOutside = outsideVolVars.density(phaseIdx);

                hOutside += rhoOutside*gOutside*scvFace_.center();
            }
        }

        volumeFlux_ = t_*(hOutside - hInside);

        if (volumeFlux > 0)
            return volumeFlux_*upwindFunction(insideVolVars, outsideVolVars);
        else
            return volumeFlux_*upwindFunction(outsideVolVars, insideVolVars);
    }

    std::set<unsigned int> stencil() const
    {
        return std::set<unsigned int>(stencil_.begin(), stencil.end());
    }

protected:


    void updateTransmissivitesAndStencil_()
    {
        if (!scvFace_.boundary())
        {
            const auto insideScvIdx = scvFace_.insideScvIdx();
            const auto& insideScv = problem_().model().fvGeometries().subControlVolume(insideScvIdx);
            const auto insideK = problem_().spatialParams().intrinsicPermeability(insideScv);
            Scalar omega1 = calculateOmega_(insideK, insideScv);

            const auto outsideScvIdx = scvFace_.outsideScvIdx();
            const auto& outsideScv = problem_().model().fvGeometries().subControlVolume(outsideScvIdx);
            const auto outsideK = problem_().spatialParams().intrinsicPermeability(outsideScv);
            Scalar omega2 = calculateOmega_(outsideK, outsideScv);

            t_ = omega1*omega2;
            t_ /= omega1 - omega2;
            t_ *= -1;

            stencil_[0] = scvFace_.insideVolVarsIdx();
            stencil_[1] = scvFace_.outsideVolVarsIdx();
        }
        else
        {
            const auto insideScvIdx = scvFace_.insideScvIdx();
            const auto& insideScv = problem_().model().fvGeometries().subControlVolume(insideScvIdx);
            const auto insideK = problem_().spatialParams().intrinsicPermeability(insideScv);
            Scalar omega1 = calculateOmega_(insideK, insideScv);

            t_ = omega1;

            stencil_[0] = scvFace_.insideVolVarsIdx();
        }
    }

    Scalar calculateOmega_(const DimWorldMatrix &K, const SubControlVolume &scv) const
    {
        GlobalPosition connectVec = scvFace_.center();
        connectVec -= scv.geometry().center();
        connectVec /= connectVec.two_norm2();

        GlobalPosition Kconnect(0);
        K.mv(connectVec, Kconnect);

        Scalar omega = Kconnect * scvFace_.normal();
        omega *= scvFace_.area()*problem_().model().curVolVars(scv).extrusionFactor();
        omega *= -1;

        return omega;
    }

    Scalar calculateOmega_(Scalar k, const SubControlVolume &scv) const
    {
        GlobalPosition connectVec = scvFace_.center();
        connectVec -= scv.geometry().center();
        connectVec /= connectVec.two_norm2();

        Scalar omega = connectVec * scvFace_.normal();
        omega *= scvFace_.area()*problem_().model().curVolVars(scv).extrusionFactor();
        omega *= k;
        omega *= -1;

        return omega;
    }

    const Problem &problem_()
    {
        return *problemPtr_;
    }

    const Problem *problem_;
    const SubControlVolumeFace *scvFace_;       //!< Pointer to the sub control volume face for which the flux variables are created
    Scalar t_;                                            //!< transmissivities for the flux calculation
    std::array<Scalar, 2> stencil_;                       //!< Indices of the cells of which the pressure is needed for the flux calculation
};

} // end namespace

#endif // DUMUX_CC_TPFA_DARCY_FLUX_VARIABLES_HH
