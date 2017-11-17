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
 *
 *
 *
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the Stokes box model.
 */

#ifndef DUMUX_STOKES_LOCAL_RESIDUAL_BASE_HH
#define DUMUX_STOKES_LOCAL_RESIDUAL_BASE_HH

#include <dune/common/float_cmp.hh>
#include <dune/grid/common/grid.hh>

#include <dumux/implicit/model.hh>
#include "properties.hh"
#include "volumevariables.hh"
#include "fluxvariables.hh"

namespace Dumux
{
/*!
 * \ingroup FemModel, instead of BoxStokesModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the local Jacobian matrix for problems
 *        using the Stokes box model.
 *
 * This class is also used for the non-isothermal and the two-component Stokes
 * model (static polymorphism).
 */
template<class TypeTag>
class StokesLocalResidual : public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
protected:
    typedef typename GET_PROP_TYPE(TypeTag, BaseLocalResidual) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    //from elastic
    using IpData = typename GET_PROP_TYPE(TypeTag, FemIntegrationPointData);

    enum {
        dim = GridView::dimension,
        numEq = GET_PROP_VALUE(TypeTag, NumEq)
    };
    enum {
        massBalanceIdx = Indices::massBalanceIdx, //!< Index of the mass balance
        momentumXIdx = Indices::momentumXIdx, //!< Index of the x-component of the momentum balance
        lastMomentumIdx = Indices::lastMomentumIdx //!< Index of the last component of the momentum balance
    };
    enum { pressureIdx = Indices::pressureIdx }; //!< Index of the pressure in a solution vector

    typedef typename Dune::ReferenceElements<Scalar, dim> ReferenceElements;
    typedef typename Dune::ReferenceElement<Scalar, dim> ReferenceElement;

    typedef Dune::FieldVector<Scalar, dim> DimVector;

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    static const bool enableUnsymmetrizedVelocityGradient = GET_PROP_VALUE(TypeTag, EnableUnsymmetrizedVelocityGradient);
    static const bool calculateNavierStokes = GET_PROP_VALUE(TypeTag, EnableNavierStokes);
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
    static const bool enablePseudo3dWallFriction = GET_PROP_VALUE(TypeTag, EnablePseudoThreeDWallFriction);

 public:
    /*!
     * \brief Constructor. Sets the upwind weight and the stabilization parameters.
     */
    StokesLocalResidual()
    {
        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        massUpwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);
        stabilizationAlpha_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Stokes, StabilizationAlpha);
        stabilizationBeta_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Stokes, StabilizationBeta);
    }

    /*!
     * method signature taken from elastic (4 parameters)
     *
     * \brief Evaluates the amount of all conservation quantities
     *        (mass and momentum) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     *  \param storage The mass of the component within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    PrimaryVariables computeStorage(const Element& element, const IpData& ipData, const VolumeVariables& volVars,
            const ElementSolution& elemSol) const
    {
    PrimaryVariables storage(0.0);
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit Euler method.
//        const ElementVolumeVariables &elemVolVars = usePrevSol ? this->prevVolVars_()
//            : this->curVolVars_();
//        const VolumeVariables &volVars = elemVolVars[scvIdx];


        if(useMoles)
            // mass balance mole fraction based
            storage[massBalanceIdx] = volVars.molarDensity();  //
        else
            // mass balance mass fraction based
            storage[massBalanceIdx] = volVars.density();

        // momentum balance
        for (int momentumIdx = momentumXIdx; momentumIdx <= lastMomentumIdx; ++momentumIdx)
            storage[momentumIdx] = volVars.density()
                * volVars.velocity()[momentumIdx-momentumXIdx];

        return PrimaryVariables(0.0); //changed from void to PrimaryVariables as its required in geomechanics
    }

    /*!
     *  method signature taken from elastic (4 parameters)
     *
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a sub-control volume. The face may be within
     *        an element (SCV face) or on the boundary. The advective and
     *        the diffusive fluxes are computed.
     *
     * \param flux The flux over the SCV (sub-control-volume) face
     * \param fIdx The index of the SCV face (may also be a boundary face)
     * \param onBoundary Indicates, if the flux is evaluated on a boundary face. If it is true,
     *        the created fluxVars object contains boundary variables evaluated at the IP of the
     *        boundary face
     */
    PrimaryVariables computeFlux(const Element& element, const IpData& ipData, const VolumeVariables& volVars,
            const ElementSolution& elemSol) const
    {
        FluxVariables fluxVars;
        fluxVars.update(this->problem_(),
                        this->element_(),
                        this->fvGeometry_(),
                        fIdx,
                        this->curVolVars_(),
                        onBoundary);
        PrimaryVariables flux(0.0);    //Addition: declare 'flux' as PrimaryVariables

        asImp_()->computeAdvectiveFlux(flux, fluxVars);
        Valgrind::CheckDefined(flux);
        asImp_()->computeDiffusiveFlux(flux, fluxVars);
        Valgrind::CheckDefined(flux);

        return flux;
    }

    /*!
     * \brief Evaluates the advective fluxes over
     *        a face of a sub-control volume.
     *
     * \todo dilatation term has to be accounted for in outflow, coupling, neumann
     *
     * \param flux The advective flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current SCV/boundary face
     *
     * An additional wall friction term can be added to account for a dimensional reduction from 3d to 2d (Kunz et al., 2016) \cite Kunz2016 <BR>
     */
    void computeAdvectiveFlux(PrimaryVariables &flux,
                              const FluxVariables &fluxVars) const
    {
//        // if the momentum balance has a dirichlet b.c., the mass balance
//        // is replaced, thus we do not need to calculate outflow fluxes here
//        if (fluxVars.onBoundary() &&
//            momentumBalanceDirichlet_(this->bcTypes_(fluxVars.upstreamIdx())))
//        {
//            return;
//        }

//        // data attached to upstream and the downstream vertices
//        const VolumeVariables &up = this->curVolVars_(fluxVars.upstreamIdx());
//        const VolumeVariables &dn = this->curVolVars_(fluxVars.downstreamIdx());
//
          DimVector massBalanceResidual = fluxVars.velocity();
//
//        if(useMoles)
//        {
//            massBalanceResidual *= (massUpwindWeight_ * up.molarDensity()
//                                    + (1.-massUpwindWeight_) * dn.molarDensity());
//        }
//        else
//        {
//            massBalanceResidual *= (massUpwindWeight_ * up.density()
//                                    + (1.-massUpwindWeight_) * dn.density());
//        }

          // adapted
          if(useMoles)
                 {
                     massBalanceResidual *= volVars.molarDensity();
                 }
                 else
                 {
                     massBalanceResidual *= volVars.density();
                 }

//        if (!fluxVars.onBoundary())
//        {
//            // stabilization of the mass balance
//            // with 0.5*alpha*(V_i + V_j)*grad P
//            DimVector stabilizationTerm = fluxVars.pressureGrad();
//            stabilizationTerm *= stabilizationAlpha_*
//                fluxVars.averageSCVVolume();
//            massBalanceResidual += stabilizationTerm;
//        }

        flux[massBalanceIdx] +=
            massBalanceResidual*fluxVars.face().normal;

        // momentum balance - pressure is evaluated as volume term
        // at the center of the SCV in computeSource
        // dynamic viscosity is upwinded

        Dune::FieldMatrix<Scalar, dim, dim> velGrad = fluxVars.velocityGrad();
        if (enableUnsymmetrizedVelocityGradient)
        {
            // nothing has to be done in this case:
            // grad v
        }
        else
        {
            // compute symmetrized gradient for the momentum flux:
            // grad v + (grad v)^T
            for (int i=0; i<dim; ++i)
                for (int j=0; j<dim; ++j)
                    velGrad[i][j] += fluxVars.velocityGrad()[j][i];
        }

        DimVector velGradComp(0.0);
        for (int velIdx = 0; velIdx < dim; ++velIdx)
        {
            velGradComp = velGrad[velIdx];

            //            velGradComp[velIdx] += 2./3*fluxVars.velocityDiv;

            velGradComp *= fluxVars.dynamicViscosity() + fluxVars.dynamicEddyViscosity();

            flux[momentumXIdx + velIdx] -=
                velGradComp*fluxVars.face().normal;

            // gravity is accounted for in computeSource; alternatively:
            //            Scalar gravityTerm = fluxVars.density *
            //                    this->problem_().gravity()[dim-1] *
            //                    fluxVars.face().ipGlobal[dim-1]*
            //                    fluxVars.face().normal[velIdx];
            //            flux[momentumXIdx + velIdx] -=
            //                    gravityTerm;

        }

        // this term changes the Stokes equation to the Navier-Stokes equation
        // rho v (v*n)
        // rho and first v are upwinded, second v is evaluated at the face
        if (calculateNavierStokes)
        {
            for (int dimIndex = 0; dimIndex < dim; ++dimIndex)
                flux[momentumXIdx + dimIndex] +=
                        volVars.density() * volVars.velocity()[dimIndex] * fluxVars.normalVelocity();
        }
    }

    /*!
     * \brief Adds the diffusive flux to the flux vector over
     *        a SCV face or a boundary face.
     *
     * It doesn't do anything in the Stokes model but is used by the
     * transport and non-isothermal models to calculate diffusive and
     * conductive fluxes.
     *
     * \param flux The diffusive flux over the SCV face or boundary face for each component
     * \param fluxVars The flux variables at the current SCV/boundary face
     */
    void computeDiffusiveFlux(PrimaryVariables &flux,
                              const FluxVariables &fluxVars) const
    { }

    /*!
     *  method signature taken from elastic (4 parameters)
     *
     * \brief Calculate the source term of all equations.
     *        The pressure gradient at the center of a SCV is computed
     *        and the gravity term evaluated.
     *
     * \param source The source/sink in the sub control volume for each component
     * \param scvIdx The local index of the sub-control volume
     */
    PrimaryVariables computeSource(const Element& element,
            const IpData& ipData,
            const SecondaryVariables& secVars,
            const ElementSolution& elemSol)
    {   PrimaryVariables source(0.0);

        const ElementVolumeVariables &elemVolVars = this->curVolVars_();
        const VolumeVariables &volVars = elemVolVars[scvIdx];

        // retrieve the source term intrinsic to the problem
        this->problem_().solDependentSource(source,
                                     this->element_(),
                                     this->fvGeometry_(),
                                     scvIdx,
                                     this->curVolVars_());

//        // ATTENTION: The source term of the mass balance has to be chosen as
//        // div (q_momentum) in the problem file
//        const Scalar alphaH2 = stabilizationAlpha_*
//            this->fvGeometry_().subContVol[scvIdx].volume;
//        source[massBalanceIdx] *= alphaH2; // stabilization of the source term


        // pressure gradient at the center of the SCV,
        // the pressure is discretized as volume term,
        // while -mu grad v is calculated in computeFlux
        DimVector pressureGradAtSCVCenter(0.0);
        DimVector grad(0.0);

        for (int scvIdx2 = 0; scvIdx2 < this->fvGeometry_().numScv; scvIdx2++)
        {
            grad = this->fvGeometry_().subContVol[scvIdx].gradCenter[scvIdx2];
            Valgrind::CheckDefined(grad);
            grad *= elemVolVars[scvIdx2].pressure();

            pressureGradAtSCVCenter += grad;
        }

        // add the component of the pressure gradient to the respective part
        // of the momentum equation and take the gravity term into account
        // signs are inverted, since q is subtracted
        for (int dimIdx = 0; dimIdx < dim; ++dimIdx)
        {
            source[momentumXIdx + dimIdx] -= pressureGradAtSCVCenter[dimIdx];
            source[momentumXIdx + dimIdx] += volVars.density()*this->problem_().gravity()[dimIdx];

//            if(enablePseudo3dWallFriction)
//            {
//                // add a wall friction term to account for a dimensional reduction from 3d to 2d
//                const auto pos = this->element_().geometry().corner(scvIdx);
//                const Scalar height = this->problem_().extrusionFactorAtPos(pos);
//                const Scalar wallFriction = 12*volVars.velocity()[dimIdx]*volVars.dynamicViscosity()/(height*height);
//                source[momentumXIdx + dimIdx] -= wallFriction;
//            }
        }

        return source;
    }

    /*!
     * \brief The Stokes model needs a modified treatment of the boundary conditions as
     *        the common box models
     */
    void evalBoundary_(const Element& element,
            const ElementGeometry& geometry,
            const LocalView& localView,
            const LocalIndexSet& localIndexSet,
            const ElementSolutionVector& curElemSol,
            const ElementSolutionVector& prevElemSol)
    {   //TODO

    }

protected:
    //functions deleted from STOKES needed for evalBoundary_


    /*!
     * \brief Returns true, if all boundary conditions for the momentum balance
     *        at the considered vertex are Dirichlet.
     */
    bool momentumBalanceDirichlet_(const BoundaryTypes& bcTypes) const
    {
        for (int momentumIdx=momentumXIdx; momentumIdx<=lastMomentumIdx; ++momentumIdx)
            if (!bcTypes.isDirichlet(momentumIdx))
                return false;
        return true;
    }


    Implementation *asImp_()
    { return static_cast<Implementation *>(this); }
    const Implementation *asImp_() const
    { return static_cast<const Implementation *>(this); }
};

}

#endif
