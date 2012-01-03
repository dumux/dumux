// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief This file contains the data which is required to calculate
 *        all fluxes of fluid phases over a face of a finite volume.
 *
 * This means pressure and temperature gradients, phase densities at
 * the integration point, etc.
 *
 */
#ifndef DUMUX_1P2C_FLUX_VARIABLES_HH
#define DUMUX_1P2C_FLUX_VARIABLES_HH

#include "1p2cproperties.hh"

#include <dumux/common/math.hh>

namespace Dumux
{

/*!
 * \ingroup OnePTwoCBoxModel
 * \ingroup BoxFluxVariables
 * \brief This template class contains the data which is required to
 *        calculate the fluxes of the fluid phases over a face of a
 *        finite volume for the one-phase, two-component model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the intergration point, etc.
 */
template <class TypeTag>
class OnePTwoCFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePTwoCIndices)) Indices;
    enum {
        phaseIdx = Indices::phaseIdx,
        comp0Idx = Indices::comp0Idx,
        comp1Idx = Indices::comp1Idx
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    enum { 
        dim = GridView::dimension,  
        dimWorld = GridView::dimensionworld 
    };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef Dune::FieldVector<Scalar, dimWorld> Vector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> Tensor;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

public:
    /*
     * \brief The constructor
     *
     * \param problem The problem
     * \param element The finite element
     * \param elemGeom The finite-volume geometry in the box scheme
     * \param scvfIdx The local index of the SCV (sub-control-volume) face
     * \param elemDat The volume variables of the current element
     */
    OnePTwoCFluxVariables(const Problem &problem,
                          const Element &element,
                          const FVElementGeometry &elemGeom,
                          int scvfIdx,
                          const ElementVolumeVariables &elemDat)
        : fvElemGeom_(elemGeom)
    {
        scvfIdx_ = scvfIdx;

        viscosityAtIP_ = Scalar(0);
        molarDensityAtIP_ = Scalar(0);
        densityAtIP_ = Scalar(0);
        potentialGrad_ = Scalar(0);
        concentrationGrad_ = Scalar(0);
        moleFracGrad_ = Scalar(0);
        massFracGrad_ = Scalar(0);

        calculateGradients_(problem, element, elemDat);
        calculateK_(problem, element, elemDat);
        calculateVelocities_(problem, element, elemDat);
        calculateDiffCoeffPM_(problem, element, elemDat);
        calculateDispersionTensor_(problem, element, elemDat);
    };

public:
    /*!
    * \brief Return the pressure potential multiplied with the
    *        intrinsic permeability which goes from vertex i to
    *        vertex j.
    *
    * Note that the length of the face's normal is the area of the
    * phase, so this is not the actual velocity by the integral of
    * the velocity over the face's area. Also note that the phase
    * mobility is not yet included here since this would require a
    * decision on the upwinding approach (which is done in the
    * actual model).
    */
   Scalar KmvpNormal() const
   { return KmvpNormal_; }

   /*!
    * \brief Return the pressure potential multiplied with the
    *        intrinsic permeability as vector (for velocity output).
    */
   Vector Kmvp() const
   { return Kmvp_; }

   /*!
    * \brief Return the subcontrol volume face.
    */
    const SCVFace &face() const
    { return fvElemGeom_.subContVolFace[scvfIdx_]; }

    /*!
     * \brief Return the intrinsic permeability tensor.
     */
    const Tensor &intrinsicPermeability() const
    { return K_; }

    /*!
     * \brief Return the dispersion tensor.
     */
    const Tensor &dispersionTensor() const
    { return dispersionTensor_; }

    /*!
     * \brief Return the pressure potential gradient.
     */
    const Vector &potentialGrad() const
    { return potentialGrad_; }

    /*!
     * \brief Return the concentration gradient.
     *
     * \param compIdx The index of the considered component
     */
    const Vector &concentrationGrad(int compIdx) const
    {
        if (compIdx != 1)
        { DUNE_THROW(Dune::InvalidStateException,
                     "The 1p2c model is supposed to need "
                     "only the concentration gradient of "
                     "the second component!"); }
        return concentrationGrad_;
    };

    /*!
        * \brief The molar concentration gradient of a component in a phase.
        *
        * \param compIdx The index of the considered component
        */
       const Vector &moleFracGrad(int compIdx) const
       {
           if (compIdx != 1)
           { DUNE_THROW(Dune::InvalidStateException,
                    "The 1p2c model is supposed to need "
                    "only the concentration gradient of "
                    "the second component!"); }
           return moleFracGrad_;
       };

       /*!
       * \brief The mass fraction gradient of a component in a phase.
       *
       * \param compIdx The index of the considered component
       */
      const Vector &massFracGrad(int compIdx) const
      {
          if (compIdx != 1)
          { DUNE_THROW(Dune::InvalidStateException,
                   "The 1p2c model is supposed to need "
                   "only the concentration gradient of "
                   "the second component!"); }
          return massFracGrad_;
      };

    /*!
    * \brief The binary diffusion coefficient for each fluid phase in the porous medium.
    */
    Scalar porousDiffCoeff() const
    {
        // TODO: tensorial diffusion coefficients
        return diffCoeffPM_;
    };

    /*!
    * \brief Return viscosity \f$\mathrm{[Pa s]}\f$ of a phase at the integration
    *        point.
    */
    Scalar viscosityAtIP() const
    { return viscosityAtIP_;}

    /*!
     * \brief Return molar density \f$\mathrm{[mol/m^3]}\f$ of a phase at the integration
     *        point.
     */
    Scalar molarDensityAtIP() const
    { return molarDensityAtIP_; }

    /*!
     * \brief Return density \f$\mathrm{[kg/m^3]}\f$ of a phase at the integration
     *        point.
     */
    Scalar densityAtIP() const
        { return densityAtIP_; }

    /*!
     * \brief Given the intrinisc permeability times the pressure
     *        potential gradient and SCV face normal for a phase,
     *        return the local index of the upstream control volume
     *        for a given phase.
     *
     *        \param normalFlux The flux over a face of the sub control volume
     */
    int upstreamIdx(Scalar normalFlux) const
    { return (normalFlux >= 0)?face().i:face().j; }

    /*!
     * \brief Given the intrinisc permeability times the pressure
     *        potential gradient and SCV face normal for a phase,
     *        return the local index of the downstream control volume
     *        for a given phase.
     *
     *        \param normalFlux The flux over a face of the sub control volume
     */
    int downstreamIdx(Scalar normalFlux) const
    { return (normalFlux > 0)?face().j:face().i; }

    /*!
    * \brief Return the local index of the upstream control volume
    *        for a given phase.
    */
   int upstreamIdx() const
   { return upstreamIdx_; }

   /*!
    * \brief Return the local index of the downstream control volume
    *        for a given phase.
    */
   int downstreamIdx() const
   { return downstreamIdx_; }

protected:

    /*!
     * \brief Calculation of the pressure and concentration gradients
     *
     *        \param problem The considered problem file
     *        \param element The considered element of the grid
     *        \param elemVolVars The parameters stored in the considered element
     */
    void calculateGradients_(const Problem &problem,
                             const Element &element,
                             const ElementVolumeVariables &elemVolVars)
    {
        const VolumeVariables &vVars_i = elemVolVars[face().i];
        const VolumeVariables &vVars_j = elemVolVars[face().j];

        Vector tmp;
        //The decision of the if-statement depends on the function useTwoPointGradient(const Element &elem,
        //int vertexI,int vertexJ) defined in test/tissue_tumor_spatialparameters.hh
        if (!problem.spatialParameters().useTwoPointGradient(element, face().i, face().j)) {
            // use finite-element gradients
            tmp = 0.0;
            for (int idx = 0;
                    idx < fvElemGeom_.numVertices;
                    idx++) // loop over adjacent vertices
            {
                // FE gradient at vertex idx
                const Vector &feGrad = face().grad[idx];

                // the pressure gradient
                tmp = feGrad;
                tmp *= elemVolVars[idx].pressure();
                potentialGrad_ += tmp;

                // the concentration gradient [mol/m^3/m]
                tmp = feGrad;
                tmp *= elemVolVars[idx].molarity(comp1Idx);
                concentrationGrad_ += tmp;

                tmp = feGrad;
                tmp *= elemVolVars[idx].moleFraction(comp1Idx);
                moleFracGrad_ += tmp;

                tmp = feGrad;
                tmp *= elemVolVars[idx].massFraction(comp1Idx);
                massFracGrad_ += tmp;
                // phase viscosity
                viscosityAtIP_ += elemVolVars[idx].viscosity()*face().shapeValue[idx];

                //phase moledensity
                molarDensityAtIP_ += elemVolVars[idx].molarDensity()*face().shapeValue[idx];

                //phase density
                densityAtIP_ += elemVolVars[idx].density()*face().shapeValue[idx];
            }
        }
        else {
            // use two-point gradients
            const GlobalPosition &posI = element.geometry().corner(face().i);
            const GlobalPosition &posJ = element.geometry().corner(face().j);
            for (int i = 0; i < Vector::size; ++ i)
                tmp[i] = posI[i] - posJ[i];
            Scalar dist = tmp.two_norm();

            tmp = face().normal;
            tmp /= face().normal.two_norm()*dist;

            potentialGrad_ = tmp;
            potentialGrad_ *= vVars_j.pressure() - vVars_i.pressure();
            concentrationGrad_ = tmp;
            concentrationGrad_ *= vVars_j.molarity(comp1Idx) - vVars_i.molarity(comp1Idx);
            moleFracGrad_ = tmp;
            moleFracGrad_ *= vVars_j.moleFraction(comp1Idx) - vVars_i.moleFraction(comp1Idx);
        }

        ///////////////
        // correct the pressure gradients by the gravitational acceleration
        ///////////////
        if (GET_PARAM(TypeTag, bool, EnableGravity)) {
            // calculate the phase density at the integration point. we
            // only do this if the wetting phase is present in both cells
            Scalar rhoI = elemVolVars[face().i].density();
            Scalar rhoJ = elemVolVars[face().j].density();
            Scalar density = (rhoI + rhoJ)/2;

            // estimate the gravitational acceleration at a given SCV face
            // using the arithmetic mean
            Vector f(problem.boxGravity(element, fvElemGeom_, face().i));
            f += problem.boxGravity(element, fvElemGeom_, face().j);
            f /= 2;

            // make it a force
            f *= density;

            // calculate the final potential gradient
            potentialGrad_ -= f;
        }
    }

    /*!
    * \brief Calculation of the harmonic mean of the intrinsic permeability
    *        uses the meanK function in the boxspatialparameters.hh file in the folder
    *        material/spatialparameters
    *
    *        \param problem The considered problem file
    *        \param element The considered element of the grid
    *        \param elemDat The parameters stored in the considered element
    */
    void calculateK_(const Problem &problem,
                     const Element &element,
                     const ElementVolumeVariables &elemDat)
    {
        const SpatialParameters &sp = problem.spatialParameters();
        sp.meanK(K_,
                 sp.intrinsicPermeability(element,
                                          fvElemGeom_,
                                          face().i),
                 sp.intrinsicPermeability(element,
                                          fvElemGeom_,
                                          face().j));

    }

    /*!
      * \brief Calculation of the velocity normal to face using Darcy's law.
      *     Tensorial permeability is multiplied with the potential Gradient and the face normal.
      *     Identify upstream node of face
      *
      *        \param problem The considered problem file
      *        \param element The considered element of the grid
      *        \param elemDat The parameters stored in the considered element
      */
    void calculateVelocities_(const Problem &problem,
                              const Element &element,
                              const ElementVolumeVariables &elemDat)
    {
        K_.mv(potentialGrad_, Kmvp_);
        KmvpNormal_ = 0;
        for (int i = 0; i < Vector::size; ++i)
            KmvpNormal_ += Kmvp_[i] * face().normal[i];
        KmvpNormal_ *= -1;

        // set the upstream and downstream vertices
        upstreamIdx_ = face().i;
        downstreamIdx_ = face().j;

        if (KmvpNormal_ < 0)
        {
            std::swap(upstreamIdx_,
                      downstreamIdx_);
        }
    }
    /*!
    * \brief Calculation of the effective diffusion coefficient
    *
    *        \param problem The considered problem file
    *        \param element The considered element of the grid
    *        \param elemDat The parameters stored in the considered element
    */
    void calculateDiffCoeffPM_(const Problem &problem,
                               const Element &element,
                               const ElementVolumeVariables &elemDat)
    {
        const VolumeVariables &vDat_i = elemDat[face().i];
        const VolumeVariables &vDat_j = elemDat[face().j];

        // Diffusion coefficient in the porous medium
        diffCoeffPM_
            = harmonicMean(vDat_i.porosity() * vDat_i.tortuosity() * vDat_i.diffCoeff(),
               vDat_j.porosity() * vDat_j.tortuosity() * vDat_j.diffCoeff());
//            = 1./2*(vDat_i.porosity() * vDat_i.tortuosity() * vDat_i.diffCoeff() +
//                    vDat_j.porosity() * vDat_j.tortuosity() * vDat_j.diffCoeff());
    }

    /*!
    * \brief Calculation of the dispersion
    *
    *        \param problem The considered problem file
    *        \param element The considered element of the grid
    *        \param elemDat The parameters stored in the considered element
    */
    void calculateDispersionTensor_(const Problem &problem,
                                    const Element &element,
                                    const ElementVolumeVariables &elemDat)
    {
        const VolumeVariables &vDat_i = elemDat[face().i];
        const VolumeVariables &vDat_j = elemDat[face().j];

        //calculate dispersivity at the interface: [0]: alphaL = longitudinal disp. [m], [1] alphaT = transverse disp. [m]
        Scalar dispersivity[2];
        dispersivity[0] = 0.5 * (vDat_i.dispersivity()[0] +  vDat_j.dispersivity()[0]);
        dispersivity[1] = 0.5 * (vDat_i.dispersivity()[1] +  vDat_j.dispersivity()[1]);

        //calculate velocity at interface: v = -1/mu * vDarcy = -1/mu * K * grad(p)
        Vector velocity;
        Valgrind::CheckDefined(potentialGrad());
        Valgrind::CheckDefined(K_);
        K_.mv(potentialGrad(), velocity);
        velocity /= - 0.5 * (vDat_i.viscosity() + vDat_j.viscosity());

        //matrix multiplication of the velocity at the interface: vv^T
        dispersionTensor_ = 0;
        for (int i=0; i<dim; i++)
            for (int j = 0; j<dim; j++)
                dispersionTensor_[i][j]=velocity[i]*velocity[j];

        //normalize velocity product --> vv^T/||v||, [m/s]
        Scalar vNorm = velocity.two_norm();

        dispersionTensor_ /= vNorm;
        if (vNorm < 1e-20)
            dispersionTensor_ = 0;

        //multiply with dispersivity difference: vv^T/||v||*(alphaL - alphaT), [m^2/s] --> alphaL = longitudinal disp., alphaT = transverse disp.
        dispersionTensor_ *= (dispersivity[0] - dispersivity[1]);

        //add ||v||*alphaT to the main diagonal:vv^T/||v||*(alphaL - alphaT) + ||v||*alphaT, [m^2/s]
        for (int i = 0; i<dim; i++)
            dispersionTensor_[i][i] += vNorm*dispersivity[1];
    }

    const FVElementGeometry &fvElemGeom_;
    int scvfIdx_;

    //! pressure potential gradient
    Vector potentialGrad_;
    //! concentratrion gradient
    Vector concentrationGrad_;
    //! molar concentratrion gradient
    Vector moleFracGrad_;
    Vector massFracGrad_;
    //! the effective diffusion coefficent in the porous medium
    Scalar diffCoeffPM_;

    //! the dispersion tensor in the porous medium
    Tensor dispersionTensor_;

    //! the intrinsic permeability tensor
    Tensor K_;
    // intrinsic permeability times pressure potential gradient
    Vector Kmvp_;
    // projected on the face normal
    Scalar KmvpNormal_;

    // local index of the upwind vertex for each phase
   int upstreamIdx_;
   // local index of the downwind vertex for each phase
   int downstreamIdx_;

    //! viscosity of the fluid at the integration point
    Scalar viscosityAtIP_;

    //! molar densities of the fluid at the integration point
    Scalar molarDensityAtIP_, densityAtIP_;
};

} // end namepace

#endif
