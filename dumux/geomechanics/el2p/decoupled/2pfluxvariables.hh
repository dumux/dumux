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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 * \brief This file contains the data which is required to calculate
 *        all fluxes of fluid phases over a face of a finite volume in the
 *        two phase discrete fracture-matrix model.
 */
#ifndef DUMUX_MODELS_2P_FLUX_VARIABLES_HH
#define DUMUX_MODELS_2P_FLUX_VARIABLES_HH

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/porousmediumflow/implicit/darcyfluxvariables.hh>

#include <dumux/porousmediumflow/2p/implicit/properties.hh>

namespace Dumux
{

/*!
 * \ingroup TwoPModel
 * \ingroup ImplicitFluxVariables
 * \brief Contains the data which is required to calculate the fluxes of
 *        the fluid phases over a face of a finite volume for the two-phase
 *        discrete fracture-matrix model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the intergration point, etc.
 */
template <class TypeTag>
class TwoPFluxVariables : public ImplicitDarcyFluxVariables<TypeTag>
{
    friend class ImplicitDarcyFluxVariables<TypeTag>; // be friends with parent
    typedef ImplicitDarcyFluxVariables<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum { dim = GridView::dimension} ;
    enum { dimWorld = GridView::dimensionworld} ;
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases)} ;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimWorldMatrix;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

    typedef Dune::FieldVector<Scalar, dim> DimVector;

public:
    /*!
     * \brief The old constructor
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the fully implicit scheme
     * \param fIdx The local index of the SCV (sub-control-volume) face
     * \param elemVolVars The volume variables of the current element
     * \param onBoundary A boolean variable to specify whether the flux variables
     * are calculated for interior SCV faces or boundary faces, default=false
     */
    DUNE_DEPRECATED_MSG("FluxVariables now have to be default constructed and updated.")
    TwoPFluxVariables(const Problem &problem,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 int fIdx,
                 const ElementVolumeVariables &elemVolVars,
                 const bool onBoundary = false)
        : ImplicitDarcyFluxVariables<TypeTag>(problem, element, fvGeometry, fIdx, elemVolVars, onBoundary) {}

    /*!
     * \brief Default constructor
     * \note This can be removed when the deprecated constructor is removed.
     */
    TwoPFluxVariables() = default;

        /*!
     * \brief Calculation of the time derivative of solid displacement
     *        \param problem The considered problem file
     *        \param element The considered element of the grid
     *        \param elemVolVars The parameters stored in the considered element
     */
    void calculateDDt_(const Problem &problem,
                    const Element &element,
                    const ElementVolumeVariables &elemVolVars)
    {
        Scalar dt = problem.timeManager().timeStepSize();

        DimVector tmp(0.0);
        // calculate time derivative of solid displacement vector
        for (int coordDir = 0; coordDir < dim; ++coordDir)
                tmp[coordDir] = problem.getdU(element, this->fvGeometry_(), this->face().i)[coordDir] / dt;

        // multiply time derivative of solid displacement vector with
        // normal vector of current scv-face
        timeDerivUNormal_ = tmp * this->face().normal;
    }

    /*!
     * \brief Actual calculation of the normal Darcy velocities.
     *
     * \param problem The problem
     * \param element The finite element
     * \param elemVolVars The volume variables of the current element
     */
    void calculateNormalVelocity_(const Problem &problem,
                                  const Element &element,
                                  const ElementVolumeVariables &elemVolVars)
    {
        // calculate the mean intrinsic permeability
        const SpatialParams &spatialParams = problem.spatialParams();
        DimWorldMatrix K;
        if (GET_PROP_VALUE(TypeTag, ImplicitIsBox))
        {
            spatialParams.meanK(K,
                                problem.getEffPermeability(element,
                                                                    this->fvGeometry_(),
                                                                    this->face().i),
                                problem.getEffPermeability(element,
                                                                    this->fvGeometry_(),
                                                                    this->face().j));
        }
        else
        {
            const Element& elementI = this->fvGeometry_().neighbors[this->face().i];
            FVElementGeometry fvGeometryI;
            fvGeometryI.subContVol[0].global = elementI.geometry().center();

            const Element& elementJ = this->fvGeometry_().neighbors[this->face().j];
            FVElementGeometry fvGeometryJ;
            fvGeometryJ.subContVol[0].global = elementJ.geometry().center();

            spatialParams.meanK(K,
                                problem.getEffPermeability(elementI, fvGeometryI, 0),
                                problem.getEffPermeability(elementJ, fvGeometryJ, 0));
        }

        // loop over all phases
        for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
        {
            // calculate the flux in the normal direction of the
            // current sub control volume face:
            //
            // v = - (K_f grad phi) * n
            // with K_f = rho g / mu K
            //
            // Mind, that the normal has the length of it's area.
            // This means that we are actually calculating
            //  Q = - (K grad phi) dot n /|n| * A


            K.mv(this->potentialGrad_[phaseIdx], this->kGradP_[phaseIdx]);
            this->kGradPNormal_[phaseIdx] = this->kGradP_[phaseIdx]*this->face().normal;

            // determine the upwind direction
            if (this->kGradPNormal_[phaseIdx] < 0)
            {
                this->upstreamIdx_[phaseIdx] = this->face().i;
                this->downstreamIdx_[phaseIdx] = this->face().j;
            }
            else
            {
                this->upstreamIdx_[phaseIdx] = this->face().j;
                this->downstreamIdx_[phaseIdx] = this->face().i;
            }

            // obtain the upwind volume variables
            const VolumeVariables& upVolVars = elemVolVars[ this->upstreamIdx(phaseIdx) ];
            const VolumeVariables& downVolVars = elemVolVars[ this->downstreamIdx(phaseIdx) ];

            // the minus comes from the Darcy relation which states that
            // the flux is from high to low potentials.
            // set the velocity
            this->velocity_[phaseIdx] = this->kGradP_[phaseIdx];
            this->velocity_[phaseIdx] *= - ( this->mobilityUpwindWeight_*upVolVars.mobility(phaseIdx)
                    + (1.0 - this->mobilityUpwindWeight_)*downVolVars.mobility(phaseIdx)) ;

            // set the volume flux
            this->volumeFlux_[phaseIdx] = this->velocity_[phaseIdx] * this->face().normal;
        }// loop all phases
    }

    Scalar timeDerivUNormal() const
    {
        return timeDerivUNormal_;
    }

    Scalar timeDerivUNormal_;
    //! change of solid displacement with time at integration point
    GlobalPosition dU_;

private:
    const FVElementGeometry* fvGeometryPtr_; //!< Information about the geometry of discretization

};

} // end namespace

#endif // DUMUX_MODELS_2P_FLUX_VARIABLES_HH
