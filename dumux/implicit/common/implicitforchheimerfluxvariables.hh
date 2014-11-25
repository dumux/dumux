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
 *        volume fluxes of fluid phases over a face of a finite volume by means
 *        of the Forchheimer approximation.
 */
#ifndef DUMUX_IMPLICIT_FORCHHEIMER_FLUX_VARIABLES_HH
#define DUMUX_IMPLICIT_FORCHHEIMER_FLUX_VARIABLES_HH

#include "implicitproperties.hh"

#include <dumux/common/parameters.hh>
#include <dumux/common/math.hh>
#include <dumux/implicit/common/implicitdarcyfluxvariables.hh>

namespace Dumux
{
    
namespace Properties
{
// forward declaration of properties 
NEW_PROP_TAG(MobilityUpwindWeight);
NEW_PROP_TAG(SpatialParams);
NEW_PROP_TAG(NumPhases);
NEW_PROP_TAG(ProblemEnableGravity);
}   

/*!
 * \ingroup ImplicitFluxVariables
 * \brief Evaluates the normal component of the Forchheimer velocity
 *        on a (sub)control volume face.
 *
 * The commonly used Darcy relation looses it's validity for \f$ Re > 1\f$.
 * If one encounters flow velocities in porous media above this Reynolds number,
 * the Forchheimer relation can be used. Like the Darcy relation, it relates
 * the gradient in potential to velocity.
 * However, this relation is not linear (as in the Darcy case) any more.
 *
 * Therefore, a Newton scheme is implemented in this class in order to calculate a
 * velocity from the current set of variables. This velocity can subsequently be used
 * e.g. by the localresidual.
 *
 * For Reynolds numbers above \f$ 500\f$ the (Standard) forchheimer relation also
 * looses it's validity.
 *
 * The Forchheimer equation looks as follows:
 * \f[ \nabla \left({p_\alpha + \rho_\alpha g z} \right)
 *     =  - \frac{\mu_\alpha}{k_{r \alpha} K} \underline{v_\alpha}
 *        - \frac{c_F}{\eta_{\alpha r} \sqrt{K}} \rho_\alpha
 *          |\underline{v_\alpha}| \underline{v_\alpha}
 * \f]
 *
 * For the formulation that is actually used in this class, see the documentation of the
 * function calculating the residual.
 *
 * - This algorithm does not find a solution if the fluid is incompressible and the
 *   initial pressure distribution is uniform.
 * - This algorithm assumes the relative passabilities to have the same functional
 *   form as the relative permeabilities.
 */
template <class TypeTag>
class ImplicitForchheimerFluxVariables
    : public ImplicitDarcyFluxVariables<TypeTag>
{
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
    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimWorldMatrix;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

public:
    /*!
     * \brief The constructor
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param fIdx The local index of the SCV (sub-control-volume) face
     * \param elemVolVars The volume variables of the current element
     * \param onBoundary A boolean variable to specify whether the flux variables
     * are calculated for interior SCV faces or boundary faces, default=false
     */
    ImplicitForchheimerFluxVariables(const Problem &problem,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const unsigned int fIdx,
                 const ElementVolumeVariables &elemVolVars,
                 const bool onBoundary = false)
        :   ImplicitDarcyFluxVariables<TypeTag>(problem, element, fvGeometry, 
                                                fIdx, elemVolVars, onBoundary)
    {
        calculateNormalVelocity_(problem, element, elemVolVars);
    }

protected:
    /*!
     * \brief Function for calculation of velocities.
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
                                spatialParams.intrinsicPermeability(element,
                                                                    this->fvGeometry_,
                                                                    this->face().i),
                                spatialParams.intrinsicPermeability(element,
                                                                    this->fvGeometry_,
                                                                    this->face().j));
        }
        else
        {
            const Element& elementI = *this->fvGeometry_.neighbors[this->face().i];
            FVElementGeometry fvGeometryI;
            fvGeometryI.subContVol[0].global = elementI.geometry().center();
            
            const Element& elementJ = *this->fvGeometry_.neighbors[this->face().j];
            FVElementGeometry fvGeometryJ;
            fvGeometryJ.subContVol[0].global = elementJ.geometry().center();
            
            spatialParams.meanK(K,
                                spatialParams.intrinsicPermeability(elementI, fvGeometryI, 0),
                                spatialParams.intrinsicPermeability(elementJ, fvGeometryJ, 0));
        }
        
        // obtain the Forchheimer coefficient from the spatial parameters
        const Scalar forchCoeff = spatialParams.forchCoeff(element,
                                                          this->fvGeometry_,
                                                          this->face().i);

        // Make sure the permeability matrix does not have off-diagonal entries
        assert( isDiagonal_(K) );

        DimWorldMatrix sqrtK(0.0);
        for (int i = 0; i < dim; ++i)
            sqrtK[i][i] = std::sqrt(K[i][i]);

        // loop over all phases
        for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
        {
            // initial guess of velocity: Darcy relation
            // first taken from base class, later overwritten in base class
            GlobalPosition velocity = this-> velocity(phaseIdx);

            GlobalPosition deltaV;           // the change in velocity between Newton iterations
            GlobalPosition residual(10e10);  // the residual (function value that is to be minimized)
            GlobalPosition tmp;              // temporary variable for numerical differentiation
            DimWorldMatrix    gradF;            // slope of equation that is to be solved

            // search by means of the Newton method for a root of Forchheimer equation
            for (int k = 0; residual.two_norm() > 1e-12 ; ++k) {

                if (k >= 30)
                    DUNE_THROW(NumericalProblem, "could not determine forchheimer velocity within "
                                                 << k << " iterations");

                // calculate current value of Forchheimer equation
                forchheimerResidual_(residual,
                                    forchCoeff,
                                    sqrtK,
                                    K,
                                    velocity,
                                    elemVolVars,
                                    this->potentialGrad_[phaseIdx],
                                    phaseIdx);

                // newton's method: slope (gradF) and current value (residual) of function is known,
                forchheimerDerivative_(gradF,
                                       forchCoeff,
                                       sqrtK,
                                       velocity,
                                       elemVolVars,
                                       phaseIdx);

                // solve for change in velocity ("x-Axis intercept")
                gradF.solve(deltaV, residual);
                velocity -= deltaV;
            }

            this->velocity_[phaseIdx]     =  velocity ; // store the converged velocity solution
            // velocity dot normal times cross sectional area is volume flux:
            this->volumeFlux_[phaseIdx]   =  velocity * this->face().normal;
        }// end loop over all phases
    }

    /*! \brief Function for calculation of Forchheimer relation.
     *
     * The relative passability \f$ \eta_r\f$ is the "Forchheimer-equivalent" to the relative
     * permeability \f$ k_r\f$.
     * We use the same function as for \f$ k_r \f$ (VG, BC, linear) other authors use a simple
     * power law e.g.: \f$\eta_{rw} = S_w^3\f$
     *
     * Some rearrangements have been made to the original Forchheimer relation:
     *
     * \f[ \underline{v_\alpha} + c_F \sqrt{K} \frac{\rho_\alpha}{\mu_\alpha }
     *     |\underline{v_\alpha}| \underline{v_\alpha}
     *     + \frac{k_{r \alpha}}{\mu_\alpha} K \nabla \left(p_\alpha
     *     + \rho_\alpha g z \right)=  0
     * \f]
     *
     * This already includes the assumption \f$ k_r(S_w) = \eta_r(S_w)\f$:
     * - \f$\eta_{rw} = S_w^x\f$ looks very similar to e.g. Van Genuchten relative permeabilities
     * - Fichot, et al. (2006), Nuclear Engineering and Design, state that several authors claim
     *   that \f$ k_r(S_w), \eta_r(S_w)\f$ can be chosen equal
     * - It leads to the equation not degenerating for the case of \f$S_w=1\f$, because I do not
     *   need to multiply with two different functions, and therefore there are terms not being
     *   zero.
     * - If this assumption is not to be made: Some regularization needs to be introduced ensuring
     *   that not all terms become zero for \f$S_w=1\f$.
     *
     *  As long as the correct velocity is not found yet, the above equation is not fulfilled, but
     *  the left hand side yields some residual. This function calculates the left hand side of the
     *  equation and returns the residual.
     *
     * \param residual The current function value for the given velocity
     * \param forchCoeff The Forchheimer coefficient
     * \param sqrtK A diagonal matrix whose entries are square roots of the permeabilities
     * \param K The permeability matrix
     * \param velocity The current velocity approximation.
     * \param elemVolVars The volume variables of the current element
     * \param potentialGrad The gradient in potential
     * \param phaseIdx The index of the currently considered phase
     */
     void forchheimerResidual_(GlobalPosition & residual,
    		 	 	 	 	 const Scalar forchCoeff,
    		 	 	 	 	 const DimWorldMatrix & sqrtK,
    		 	 	 	 	 const DimWorldMatrix & K,
    		 	 	 	 	 const GlobalPosition & velocity,
    		 	 	 	 	 const ElementVolumeVariables & elemVolVars,
    		 	 	 	 	 const GlobalPosition & potentialGrad,
    		 	 	 	 	 const unsigned int phaseIdx) const
     {
         const VolumeVariables upVolVars    = elemVolVars[this->upstreamIdx(phaseIdx)];
         const VolumeVariables downVolVars  = elemVolVars[this->downstreamIdx(phaseIdx)];

         // Obtaining the physical quantities
         Scalar alpha = this->mobilityUpwindWeight_;
         const Scalar mobility = alpha * upVolVars.mobility(phaseIdx)
                 + (1. - alpha)*  downVolVars.mobility(phaseIdx) ;

         const Scalar viscosity = alpha * upVolVars.fluidState().viscosity(phaseIdx)
                 + (1. - alpha)*  downVolVars.fluidState().viscosity(phaseIdx) ;

         const Scalar density = alpha * upVolVars.fluidState().density(phaseIdx)
                 + (1. - alpha) * downVolVars.fluidState().density(phaseIdx) ;

         // residual  = v
         residual = velocity;

         // residual += k_r/mu  K grad p
         K.usmv(mobility, potentialGrad, residual);

         // residual += c_F rho / mu abs(v) sqrt(K) v
         sqrtK.usmv(forchCoeff * density / viscosity * velocity.two_norm(), velocity, residual);
     }


     /*! \brief Function for calculation of the gradient of Forchheimer relation with respect
      * to velocity.
      *
      * This function already exploits that \f$ \sqrt{K}\f$ is a diagonal matrix. Therefore,
      * we only have to deal with the main entries.
      * The gradient of the Forchheimer relations looks as follows (mind that \f$ \sqrt{K}\f$ 
      * is a tensor):
      *
      * \f[  f\left(\underline{v_\alpha}\right) =
      * \left(
      * \begin{array}{ccc}
      * 1 & 0 &0 \\
      * 0 & 1 &0 \\
      * 0 & 0 &1 \\
      * \end{array}
      * \right)
      * +
      * c_F \frac{\rho_\alpha}{\mu_\alpha} |\underline{v}_\alpha| \sqrt{K}
      * +
      * c_F \frac{\rho_\alpha}{\mu_\alpha}\frac{1}{|\underline{v}_\alpha|} \sqrt{K}
      * \left(
      * \begin{array}{ccc}
      * v_x^2 & v_xv_y & v_xv_z \\
      * v_yv_x & v_{y}^2 & v_yv_z \\
      * v_zv_x & v_zv_y &v_{z}^2 \\
      * \end{array}
      * \right)
      *  \f]
      *
      * \param derivative The gradient of the forchheimer equation
      * \param forchCoeff Forchheimer coefficient
      * \param sqrtK A diagonal matrix whose entries are square roots of the permeabilities.
      * \param velocity The current velocity approximation.
      * \param elemVolVars The volume variables of the current element
      * \param phaseIdx The index of the currently considered phase
      */
     void forchheimerDerivative_(DimWorldMatrix & derivative,
    		 	 	 	 	 	 const Scalar forchCoeff,
    		 	 	 	 	 	 const DimWorldMatrix & sqrtK,
    		 	 	 	 	 	 const GlobalPosition & velocity,
    		 	 	 	 	 	 const ElementVolumeVariables & elemVolVars,
    		 	 	 	 	 	 const unsigned int phaseIdx) const
     {
         const VolumeVariables upVolVars    = elemVolVars[this->upstreamIdx(phaseIdx)];
         const VolumeVariables downVolVars  = elemVolVars[this->downstreamIdx(phaseIdx)];

         // Obtaining the physical quantities
         Scalar alpha = this->mobilityUpwindWeight_;
         const Scalar viscosity = alpha * upVolVars.fluidState().viscosity(phaseIdx)
                 + (1. - alpha)*  downVolVars.fluidState().viscosity(phaseIdx) ;

         const Scalar density = alpha * upVolVars.fluidState().density(phaseIdx)
                 + (1. - alpha) * downVolVars.fluidState().density(phaseIdx) ;

         // Initialize because for low velocities, we add and do not set the entries.
         derivative = 0.;

         const Scalar absV = velocity.two_norm() ;
         // This part of the derivative is only calculated if the velocity is sufficiently small:
         // do not divide by zero!
         // This should not be very bad because the derivative is only intended to give an
         // approximation of the gradient of the
         // function that goes into the Newton scheme.
         // This is important in the case of a one-phase region in two-phase flow. The non-existing
         // phase has a velocity of zero (kr=0).
         // f = sqrtK * vTimesV (this is  a matrix)
         // f *= forchCoeff density / |velocity| / viscosity (multiply all entries with scalars)
         if(absV > 1e-20)
             for (int i=0; i<dim; i++)
                 for (int k=0; k<dim; k++)
                     derivative[i][k]= sqrtK[i][i] * velocity[i] * velocity[k] * forchCoeff 
                                       * density /  absV  / viscosity;

         // add on the main diagonal:
         // 1 + sqrtK_i forchCoeff density |v| / viscosity
         for (int i=0; i<dim; i++)
             derivative[i][i] += (1.0 + ( sqrtK[i][i]*forchCoeff * density * absV / viscosity ) ) ;
     }

     /*!
      * \brief Check whether all off-diagonal entries of a tensor are zero.
      *
      * \param K the tensor that is to be checked.
      *
      * \return True if all off-diagonals are zero.
      *
     */
     const bool isDiagonal_(const DimWorldMatrix & K) const
     {
         for (int i = 0; i < dim; i++) {
             for (int k = 0; k < dim; k++) {
                 if ((i != k) && (std::abs(K[i][k]) >= 1e-25)) {
                   return false;
                 }
             }
         }
         return true;
     }
};

} // end namespace

#endif
