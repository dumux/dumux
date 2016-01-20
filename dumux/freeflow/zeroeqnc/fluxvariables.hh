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
 *        the fluxes of the compositional ZeroEq model over a face of a finite volume.
 *
 * This means the methods to calculate the eddy diffusivity.
 */
#ifndef DUMUX_ZEROEQNC_FLUX_VARIABLES_HH
#define DUMUX_ZEROEQNC_FLUX_VARIABLES_HH

#include <dumux/common/math.hh>
#include <dumux/common/valgrind.hh>

#include <dumux/freeflow/zeroeq/fluxvariables.hh>
#include <dumux/freeflow/zeroeqnc/indices.hh>
#include <dumux/freeflow/zeroeqnc/properties.hh>

namespace Dumux
{

/*!
 * \ingroup BoxZeroEqncModel
 * \ingroup ImplicitFluxVariables
 * \brief This template class contains data which is required to
 *        calculate the component fluxes over a face of a finite
 *        volume for a compositional ZeroEq model.
 *
 * This means the methods to calculate the eddy diffusivity.
 */
template <class TypeTag>
class ZeroEqncFluxVariables : public ZeroEqFluxVariables<TypeTag>
{
    typedef ZeroEqFluxVariables<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { dim = GridView::dimension,
           phaseIdx = Indices::phaseIdx,
           transportCompIdx = Indices::transportCompIdx};

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dim> DimVector;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

public:
    ZeroEqncFluxVariables(const Problem &problem,
                          const Element &element,
                          const FVElementGeometry &fvGeometry,
                          const int fIdx,
                          const ElementVolumeVariables &elemVolVars,
                          const bool onBoundary = false)
        : ParentType(problem, element, fvGeometry, fIdx, elemVolVars, onBoundary)
        , flowNormal_(GET_PARAM_FROM_GROUP(TypeTag, int, ZeroEq, FlowNormal))
        , wallNormal_(GET_PARAM_FROM_GROUP(TypeTag, int, ZeroEq, WallNormal))
        , eddyDiffusivityModel_(GET_PARAM_FROM_GROUP(TypeTag, int, ZeroEq, EddyDiffusivityModel))
    {
        globalPos_ = this->face().ipGlobal;
        posIdx_ = problem.model().getPosIdx(globalPos_);
        wallIdx_ = problem.model().getWallIdx(globalPos_, posIdx_);
        velGrad_ = this->velocityGrad_[flowNormal_][wallNormal_];
        velGradWall_ = problem.model().wall[wallIdx_].wallVelGrad[posIdx_];
        yPlusReal_ = this->distanceToWallReal() * this->frictionVelocityWall()
                     / problem.model().wall[wallIdx_].wallKinematicViscosity[posIdx_];

        eddyDiffusivity_ = 0.0;
        mixingLengthDiffusivity_ = 0.0;

        DimVector tmp(0.0);
        densityGrad_ = 0.0;

        // calculate gradients and secondary variables at IPs
        for (int idx = 0; idx < this->fvGeometry_.numScv; idx++)
        {
            tmp = this->face().grad[idx];
            tmp *= elemVolVars[idx].density();
            densityGrad_ += tmp;
        }

        // Richardson number
        // Schlichting, Boundary Layer Theory, 1997, p472
        Scalar gravity = 0.0;
        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
            gravity = problem.gravity()[wallNormal_];
        else
            gravity = 9.81;
        richardsonNumber_ = -gravity / this->density()
                            * densityGrad()[wallNormal_] / (velGradWall_ * velGradWall_);

        // calculation of an eddy diffusivity only makes sense with Navier-Stokes equation
        if (GET_PROP_VALUE(TypeTag, EnableNavierStokes))
            calculateEddyDiffusivity_(problem, element, elemVolVars);
    };

protected:
    /*!
     * \brief This functions calculates the eddy diffusivity.
     *
     * The eddy diffusivity is added to the diffusivity in stokesnclocalresidual.hh
     * at each scv face.
     */
    void calculateEddyDiffusivity_(const Problem &problem,
                                   const Element &element,
                                   const ElementVolumeVariables &elemVolVars)
    {
        // no eddy diffusivity model
        if (eddyDiffusivityModel_ == EddyDiffusivityIndices::noEddyDiffusivityModel)
            return;

        // Reynolds analogy
        // Bird, Stewart and Lightfoot, Transport Phenomena, 2007, p 657
        else if (eddyDiffusivityModel_ == EddyDiffusivityIndices::reynoldsAnalogy)
        {
            eddyDiffusivity_ = this->kinematicEddyViscosity()
                               / GET_PARAM_FROM_GROUP(TypeTag, Scalar, ZeroEq, TurbulentSchmidtNumber);
        }

        // modified Van-Driest
        // e.g. Bird, Stewart, and Lightfoot, E. N. Transport phenomena, 2007
        else if (eddyDiffusivityModel_ == EddyDiffusivityIndices::modifiedVanDriest)
        {
            Scalar aPlus = 26.0;
            Scalar bPlus = 0.26;
            // eddy diffusivity can only be calculated correctly for non-zero distance to walls
            mixingLengthDiffusivity_ = 0.0;
            if (this->distanceToWallReal() > 0.0 && yPlusReal_ > 0.0)
                mixingLengthDiffusivity_ = this->karmanConstant() * this->distanceToWallReal()
                                           * (1.0 - std::exp(-yPlusReal_ / aPlus))
                                           / std::sqrt(1.0 - std::exp(-bPlus * yPlusReal_));

            eddyDiffusivity_ = mixingLengthDiffusivity_ * mixingLengthDiffusivity_ * std::abs(velGrad_);
        }

        // Deissler near wall law
        // Deissler, R. G. "Analysis of Turbulent Heat Transfer, Mass Transfer, and Friction in Smooth Tubes at High Prandtl and Schmidt Numbers"
        //   NACA Report, 1954, 1210, 69-82
        else if (eddyDiffusivityModel_ == EddyDiffusivityIndices::deissler)
        {
            const Scalar deisslerConstant = 0.124;
            const Scalar beta = this->density() * deisslerConstant * deisslerConstant
                                * std::abs(this->velocity()[flowNormal_])
                                * this->distanceToWallReal();
            eddyDiffusivity_ = beta * (1.0 - std::exp(-beta / this->dynamicViscosity()));
            eddyDiffusivity_ /= this->density();
        }

        // Meier and Rotta
        // Cebeci, Analysis of turbulent boundary layer, 1974, p 251f
        else if (eddyDiffusivityModel_ == EddyDiffusivityIndices::meier)
        {
            // Sc_t at flow = 0.86
            // Sc_t in wall = 1.34
            Scalar kappaMeier = this->karmanConstant() / std::sqrt(0.86);
            Scalar aPlusMeier = std::sqrt(1.34) / std::sqrt(0.86) * 26.0;
            mixingLengthDiffusivity_ = 0.0;
            if (this->distanceToWallReal() > 0.0 && yPlusReal_ > 0.0)
                mixingLengthDiffusivity_ = kappaMeier * this->distanceToWallReal()
                                           * (1.0 - std::exp(- yPlusReal_ / aPlusMeier));
            eddyDiffusivity_  = mixingLengthDiffusivity_ * mixingLengthDiffusivity_ * std::abs(velGrad_);
        }

        // exponential (after Mamayev)
        // Lehfeldt, Ein algebraisches Turbulenzmodell für Ästuare, 1991, p 65
        else if (eddyDiffusivityModel_ == EddyDiffusivityIndices::exponential)
        {
            Scalar m = 0.8; // 18.0=Perrels, 0.8=Mamayev
            if (velGradWall_ == 0) // means incredible high Ri
                eddyDiffusivity_ = 0.0;
            else
                eddyDiffusivity_ = std::exp(-m * richardsonNumber())
                                   * this->kinematicEddyViscosity();
        }

        else
        {
            DUNE_THROW(Dune::NotImplemented, "This eddy diffusivity model is not implemented.");
        }

        Valgrind::CheckDefined(eddyDiffusivity_);
    }


public:
    /*!
     * \brief Returns the eddy diffusivity \f$\mathrm{[m^2/s]}\f$.
     *
     * The eddy diffusivity does not depend on the component
     */
    const Scalar eddyDiffusivity() const
    { return eddyDiffusivity_; }

    /*!
     * \brief Returns the density gradient \f$\mathrm{[kg/m^4]}\f$.
     */
    const DimVector densityGrad() const
    { return densityGrad_; }

    /*!
     * \brief Returns the richardson number \f$\mathrm{[-]}\f$.
     */
    const Scalar richardsonNumber() const
    { return richardsonNumber_; }

    /*!
     * \brief Returns the Schmidt number (molecular) \f$\mathrm{[-]}\f$.
     */
    const Scalar schmidtNumber() const
    { return this->kinematicViscosity() / this->diffusionCoeff(transportCompIdx); }

    /*!
     * \brief Returns the turbulent Schmidt number \f$\mathrm{[-]}\f$.
     */
    const Scalar turbulentSchmidtNumber() const
    { return this->kinematicEddyViscosity() / eddyDiffusivity(); }

private:
        const int flowNormal_;
        const int wallNormal_;
        const int eddyDiffusivityModel_;
        int posIdx_;
        int wallIdx_;

        Scalar velGrad_;
        Scalar velGradWall_;
        Scalar yPlusReal_;
        DimVector globalPos_;

        Scalar eddyDiffusivity_;
        Scalar mixingLengthDiffusivity_;
        Scalar richardsonNumber_;
        DimVector densityGrad_;
};

} // end namespace

#endif // DUMUX_ZEROEQNC_FLUX_VARIABLES_HH
