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
 *        the fluxes of the non-isothermal compositional ZeroEq model over
 *        a face of a finite volume.
 *
 * This means the methods to calculate the eddy conductivity.
 */
#ifndef DUMUX_ZEROEQNCNI_FLUX_VARIABLES_HH
#define DUMUX_ZEROEQNCNI_FLUX_VARIABLES_HH

#include <dumux/common/math.hh>
#include <dumux/common/valgrind.hh>

#include <dumux/freeflow/zeroeqnc/fluxvariables.hh>
#include <dumux/freeflow/zeroeqncni/indices.hh>
#include <dumux/freeflow/zeroeqncni/properties.hh>

namespace Dumux
{

/*!
 * \ingroup BoxZeroEqncniModel
 * \ingroup ImplicitFluxVariables
 * \brief This template class contains data which is required to
 *        calculate the component fluxes over a face of a finite
 *        volume for a non-isothermal compositional ZeroEq model.
 *
 * This means the methods to calculate the eddy conductivity.
 */
template <class TypeTag>
class ZeroEqncniFluxVariables : public ZeroEqncFluxVariables<TypeTag>
{
    typedef ZeroEqncFluxVariables<TypeTag> ParentType;
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
    /*!
     * \brief Compute / update the flux variables
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param fIdx The local index of the SCV (sub-control-volume) face
     * \param elemVolVars The volume variables of the current element
     * \param onBoundary A boolean variable to specify whether the flux variables
     * are calculated for interior SCV faces or boundary faces, default=false
     */
    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const int fIdx,
                const ElementVolumeVariables &elemVolVars,
                const bool onBoundary = false)
    {
        ParentType::update(problem, element, fvGeometry, fIdx, elemVolVars, onBoundary);

        flowNormal_ = GET_PARAM_FROM_GROUP(TypeTag, int, ZeroEq, FlowNormal);
        wallNormal_ = GET_PARAM_FROM_GROUP(TypeTag, int, ZeroEq, WallNormal);
        eddyConductivityModel_ = GET_PARAM_FROM_GROUP(TypeTag, int, ZeroEq, EddyConductivityModel);

        globalPos_ = this->face().ipGlobal;
        posIdx_ = problem.model().getPosIdx(globalPos_);
        wallIdx_ = problem.model().getWallIdx(globalPos_, posIdx_);
        velGrad_ = this->velocityGrad_[flowNormal_][wallNormal_];
        velGradWall_ = problem.model().wall[wallIdx_].wallVelGrad[posIdx_];
        yPlusReal_ = this->distanceToWallReal() * this->frictionVelocityWall()
                     / problem.model().wall[wallIdx_].wallKinematicViscosity[posIdx_];

        temperatureEddyConductivity_ = 0.0;
        mixingLengthConductivity_ = 0.0;

        // calculation of an eddy conductivity only makes sense with Navier-Stokes equation
        if (GET_PROP_VALUE(TypeTag, EnableNavierStokes))
            calculateEddyConductivity_(problem, element, elemVolVars);
    }

protected:
    /*!
     * \brief This functions calculates the eddy conductivity.
     *
     * The eddy conductivity is added to the conductivity in stokesncnilocalresidual.hh
     * at each scv face.
     */
    void calculateEddyConductivity_(const Problem &problem,
                                    const Element &element,
                                    const ElementVolumeVariables &elemVolVars)
    {
        // IMPORTANT:
        // the temperatureEddyConductivity_ a_t [m^2/s] is converted to
        // thermalEddyConductivity \lambda_t [W/(m K)] by the convenience function

        // no eddy conductivity model
        if (eddyConductivityModel_ == EddyConductivityIndices::noEddyConductivityModel)
            return;

        // Reynolds analogy
        // Bird, Stewart and Lightfoot, Transport Phenomena, 2007, p 657
        else if (eddyConductivityModel_ == EddyConductivityIndices::reynoldsAnalogy)
        {
            temperatureEddyConductivity_ = this->kinematicEddyViscosity()
                                           / GET_PARAM_FROM_GROUP(TypeTag, Scalar, ZeroEq, TurbulentPrandtlNumber);
        }

        // modified Van-Driest
        // e.g. Bird, Stewart, and Lightfoot, E. N. Transport phenomena, 2007
        else if (eddyConductivityModel_ == EddyConductivityIndices::modifiedVanDriest)
        {
            Scalar aPlus = 26.0;
            Scalar bPlus = 0.26;
            // eddy conductivity can only be calculated correctly for non-zero distance to walls
            mixingLengthConductivity_ = 0.0;
            if (this->distanceToWallReal() > 0.0 && yPlusReal_ > 0.0)
                mixingLengthConductivity_ = this->karmanConstant() * this->distanceToWallReal()
                                            * (1.0 - std::exp(-yPlusReal_ / aPlus))
                                            / std::sqrt(1.0 - std::exp(-bPlus * yPlusReal_));

            temperatureEddyConductivity_ = mixingLengthConductivity_ * mixingLengthConductivity_
                                           * std::abs(velGrad_);
        }

        // Deissler near wall law
        // Deissler, R. G. "Analysis of Turbulent Heat Transfer, Mass Transfer, and Friction in Smooth Tubes at High Prandtl and Schmidt Numbers"
        //   NACA Report, 1954, 1210, 69-82
        else if (eddyConductivityModel_ == EddyConductivityIndices::deissler)
        {
            const Scalar deisslerConstant = 0.124;
            const Scalar beta = this->density() * deisslerConstant * deisslerConstant
                                * std::abs(this->velocity()[flowNormal_])
                                * this->distanceToWallReal();
            temperatureEddyConductivity_ = beta * (1.0 - std::exp(-beta / this->dynamicViscosity()));
            temperatureEddyConductivity_ /= this->density();
        }

        // Meier and Rotta
        // Cebeci, Analysis of turbulent boundary layer, 1974, p 251f
        else if (eddyConductivityModel_ == EddyConductivityIndices::meier)
        {
            // Pr_t at flow = 0.86
            // Pr_t in wall = 1.34
            Scalar kappaMeier = this->karmanConstant() / std::sqrt(0.86);
            Scalar aPlusMeier = std::sqrt(1.34) / std::sqrt(0.86) * 26.0;
            mixingLengthConductivity_ = 0.0;
            if (this->distanceToWallReal() > 0.0 && yPlusReal_ > 0.0)
                mixingLengthConductivity_ = kappaMeier * this->distanceToWallReal()
                                            * (1.0 - std::exp(- yPlusReal_ / aPlusMeier));
            temperatureEddyConductivity_  = mixingLengthConductivity_ * mixingLengthConductivity_ * std::abs(velGrad_);
        }

        else
        {
            DUNE_THROW(Dune::NotImplemented, "This eddy conductivity model is not implemented.");
        }

        Valgrind::CheckDefined(temperatureEddyConductivity_);
    }


public:
    /*!
     * \brief Returns the thermal eddy conductivity \f$\mathrm{[W/(m*K)]}\f$.
     */
    const Scalar thermalEddyConductivity() const
    { return temperatureEddyConductivity() * this->density() * this->heatCapacity(); }

    /*!
     * \brief Returns the temperature eddy conductivity \f$\mathrm{[m^2/s]}\f$.
     */
    const Scalar temperatureEddyConductivity() const
    { return temperatureEddyConductivity_; }

    /*!
     * \brief Returns Prandtl number (molecular) \f$\mathrm{[-]}\f$.
     */
    const Scalar prandtlNumber() const
    { return this->dynamicViscosity() * this->heatCapacity() / this->thermalConductivity(); }

    /*!
     * \brief Returns the turbulent Prandtl number \f$\mathrm{[-]}\f$.
     */
    const Scalar turbulentPrandtlNumber() const
    { return this->dynamicEddyViscosity() * this->heatCapacity() / thermalEddyConductivity(); }


private:
    int flowNormal_;
    int wallNormal_;
    int eddyConductivityModel_;
    int posIdx_;
    int wallIdx_;

    Scalar velGrad_;
    Scalar velGradWall_;
    DimVector globalPos_;
    Scalar yPlusReal_;

    Scalar temperatureEddyConductivity_;
    Scalar mixingLengthConductivity_;
};

} // end namespace

#endif // DUMUX_ZEROEQNCNI_FLUX_VARIABLES_HH
