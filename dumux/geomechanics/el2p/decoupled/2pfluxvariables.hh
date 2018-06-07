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
#ifndef DUMUX_DECOUPLED_TWOP_FLUX_VARIABLES_HH
#define DUMUX_DECOUPLED_TWOP_FLUX_VARIABLES_HH

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
class DecoupledTwoPFluxVariables : public ImplicitDarcyFluxVariables<TypeTag>
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
    DecoupledTwoPFluxVariables(const Problem &problem,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 int fIdx,
                 const ElementVolumeVariables &elemVolVars,
                 const bool onBoundary = false)
        : ParentType(problem, element, fvGeometry, fIdx, elemVolVars, onBoundary)
        {
            calculateDDt_(problem, element, elemVolVars);
            calculateNormalVelocity_(problem, element, elemVolVars);
        }

    /*!
     * \brief Default constructor
     * \note This can be removed when the deprecated constructor is removed.
     */
    DecoupledTwoPFluxVariables() = default;

    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const int fIdx,
                const ElementVolumeVariables &elemVolVars,
                const bool onBoundary = false)
    {
        ParentType::update(problem, element, fvGeometry, fIdx, elemVolVars, onBoundary);
//         calculateDDt_(problem, element, elemVolVars);
        calculateNormalVelocity_(problem, element, elemVolVars);

        onBoundary_ = onBoundary;
    }

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
        DimWorldMatrix Keff_i, Keff_j;

        Scalar factor_i, factor_j;

        if (GET_PROP_VALUE(TypeTag, ImplicitIsBox))
        {
//             Scalar exp_i, exp_j;
//             exp_i =
//                     22.2
//                             * (elemVolVars[this->face().i].effPorosity()
//                                 / elemVolVars[this->face().i].initialPorosity()
//                                     - 1);
//             exp_j =
//                     22.2
//                             * (elemVolVars[this->face().j].effPorosity()
//                                 / elemVolVars[this->face().j].initialPorosity()
//                                     - 1);

            factor_i = pow( (elemVolVars[this->face().i].effPorosity()
                                / elemVolVars[this->face().i].initialPorosity()), 15);

            factor_j = pow( (elemVolVars[this->face().j].effPorosity()
                                / elemVolVars[this->face().j].initialPorosity()), 15);

            Keff_i = spatialParams.intrinsicPermeability(element, this->fvGeometry_(), this->face().i);

            Keff_j = spatialParams.intrinsicPermeability(element, this->fvGeometry_(), this->face().j);

//             Keff_i *= exp(exp_i);
//             Keff_j *= exp(exp_j);

//             Keff_i *= factor_i;
//             Keff_j *= factor_j;

            spatialParams.meanK(Keff_,
                                Keff_i,
                                Keff_j);

//             spatialParams.meanK(K,
//                                 problem.getEffPermeability(element,
//                                                                     this->fvGeometry_(),
//                                                                     this->face().i),
//                                 problem.getEffPermeability(element,
//                                                                     this->fvGeometry_(),
//                                                                     this->face().j));
        }
        else
        {
            const Element& elementI = this->fvGeometry_().neighbors[this->face().i];
            FVElementGeometry fvGeometryI;
            fvGeometryI.subContVol[0].global = elementI.geometry().center();

            const Element& elementJ = this->fvGeometry_().neighbors[this->face().j];
            FVElementGeometry fvGeometryJ;
            fvGeometryJ.subContVol[0].global = elementJ.geometry().center();

            factor_i = pow( (elemVolVars[this->face().i].effPorosity()
                                / elemVolVars[this->face().i].initialPorosity()), 15);

            factor_j = pow( (elemVolVars[this->face().j].effPorosity()
                                / elemVolVars[this->face().j].initialPorosity()), 15);

            Keff_i = spatialParams.intrinsicPermeability(elementI, fvGeometryI, this->face().i);

            Keff_j = spatialParams.intrinsicPermeability(elementJ, fvGeometryJ, this->face().j);

//             Keff_i *= factor_i;
//             Keff_j *= factor_j;

            spatialParams.meanK(Keff_,
                                Keff_i,
                                Keff_j);

            int eIdx = problem.model().elementMapper().index(element);
            if (eIdx == 210)
            {
//                 std::cout << "phi_init_i = " << elemVolVars[this->face().i].initialPorosity() << std::endl;
//                 std::cout << "phi_init_j = " << elemVolVars[this->face().j].initialPorosity() << std::endl;

//                 std::cout << "elementI = " << problem.model().elementMapper().index(elementI) << std::endl;
//                 std::cout << "elementJ = " << problem.model().elementMapper().index(elementJ) << std::endl;
//
//                 std::cout << "this->face().i = " << this->face().i << std::endl;
//                 std::cout << "this->face().j = " << this->face().i << std::endl;
// //
//                 std::cout << "effPorosity[" << this->face().i <<"] = " << elemVolVars[this->face().i].effPorosity() << std::endl;
//                 std::cout << "effPorosity[" << this->face().j <<"]  = " << elemVolVars[this->face().j].effPorosity() << std::endl;
//
//                 std::cout << "Keff_i = " << Keff_i << std::endl;
//                 std::cout << "Keff_j = " << Keff_j << std::endl;
            }

//             spatialParams.meanK(K,
//                                 problem.getEffPermeability(elementI, fvGeometryI, 0),
//                                 problem.getEffPermeability(elementJ, fvGeometryJ, 0));
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


            Keff_.mv(this->potentialGrad_[phaseIdx], this->kGradP_[phaseIdx]);
            this->kGradPNormal_[phaseIdx] = this->kGradP_[phaseIdx]*this->face().normal;

//             if ((problem.coupled() == true) && (phaseIdx == 0))
//             {
//                 std::cout << "FluxVars at element " << eIdx << std::endl;
//                 std::cout << "potentialGrad_ = " << this->potentialGrad_[phaseIdx] << std::endl;
// //                 std::cout << "Keff = " << Keff << std::endl;
// //                 std::cout << "kGradP_[phaseIdx] = " << this->kGradP_[phaseIdx] << std::endl;
//                 std::cout << std::endl;
//             }

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

//             int eIdx = problem.model().elementMapper().index(element);
//             if ((problem.coupled() == true) && (eIdx == 2925))
//             {
// //                 std::cout << "upVolVars.mobility(phaseIdx) = " << this->upVolVars.mobility(phaseIdx) << std::endl;
// //                 std::cout << std::endl;
//                 std::cout << "upVolVars.kr(" << phaseIdx << ") at " << this->upstreamIdx_[phaseIdx] << " = " << upVolVars.mobility(phaseIdx) * upVolVars.fluidState().viscosity(phaseIdx) << std::endl;
//                 std::cout << std::endl;
//
//                 std::cout << "saturation(" << phaseIdx << ") at " << upVolVars.fluidState().saturation(phaseIdx) << std::endl;
//                 std::cout << std::endl;
//             }

            // set the volume flux
            this->volumeFlux_[phaseIdx] = this->velocity_[phaseIdx] * this->face().normal;
        }// loop all phases
    }

    Scalar timeDerivUNormal() const
    {
        return timeDerivUNormal_;
    }

    DimWorldMatrix Keff() const
    {
        return Keff_;
    }
    /*
     * \brief Return the gradient of the potential for each phase.
     */
    DimVector potentialGrad(int phaseIdx) const
    {
        return this->potentialGrad_[phaseIdx];
    }

    Scalar kGradPNormal(int phaseIdx) const
    {
        return this->kGradPNormal_[phaseIdx];
    }

    Scalar timeDerivUNormal_;
    //! change of solid displacement with time at integration point
    GlobalPosition dU_;
    DimWorldMatrix Keff_;

    /*!
     * \brief Indicates if a face is on a boundary. Used for in the
     *        face() method (e.g. for outflow boundary conditions).
     */
    bool onBoundary() const
    { return onBoundary_; }
private:
    const FVElementGeometry* fvGeometryPtr_; //!< Information about the geometry of discretization
    bool onBoundary_;


};

} // end namespace

#endif // DUMUX_MODELS_2P_FLUX_VARIABLES_HH
