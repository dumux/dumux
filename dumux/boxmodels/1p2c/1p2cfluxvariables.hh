// $Id$
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
 * \ingroup OnePTwoCBoxModel
 */
#ifndef DUMUX_1P2C_FLUX_VARIABLES_HH
#define DUMUX_1P2C_FLUX_VARIABLES_HH

#include "1p2cproperties.hh"

#include <dumux/common/math.hh>

namespace Dumux
{

/*!
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
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    typedef Dune::FieldVector<Scalar, dimWorld> Vector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> Tensor;

    typedef typename FVElementGeometry::SubControlVolume SCV;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

    typedef typename GridView::template Codim<0>::Entity Element;

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

        calculateGradients_(problem, element, elemDat);
        calculateK_(problem, element, elemDat);
        calculateDiffCoeffPM_(problem, element, elemDat);
        calculateDispersionTensor_(problem, element, elemDat);
    };

public:
    const SCVFace &face() const
    { return fvElemGeom_.subContVolFace[scvfIdx_]; }

    /*!
     * \brief Return the intrinsic permeability.
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

    Scalar porousDiffCoeff() const
    {
        // TODO: tensorial diffusion coefficients
        return diffCoeffPM_;
    };

    Scalar viscosityAtIP() const
    { return viscosityAtIP_;}

    Scalar molarDensityAtIP() const
    { return molarDensityAtIP_; }


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

protected:

        /*!
         * \brief Calculation of the pressure and concentration gradient
         *
         *        \param problem The considered problem file
         *        \param element The considered element of the grid
         *        \param elemDat The parameters stored in the considered element
         */
    void calculateGradients_(const Problem &problem,
                             const Element &element,
                             const ElementVolumeVariables &elemDat)
    {
        const VolumeVariables &vVars_i = elemDat[face().i];
        const VolumeVariables &vVars_j = elemDat[face().j];

        potentialGrad_ = 0.0;
        concentrationGrad_ = 0.0;

        Vector tmp;
        //The decision of the if-statement depends on the function useTwoPointGradient(const Element &elem,
        //int vertexI,int vertexJ) defined in test/tissue_tumor_spatialparameters.hh
        if (!problem.spatialParameters().useTwoPointGradient(element, face().i, face().j)) {
            // use finite-element gradients
            tmp = 0.0;
            int n = element.template count<dim>();
            for (int idx = 0; idx < n; idx++) // loop over adjacent vertices
            {
                // FE gradient at vertex idx
                const Vector &feGrad = face().grad[idx];

                // the pressure gradient
                tmp = feGrad;
                tmp *= elemDat[idx].pressure();
                potentialGrad_ += tmp;

                // the concentration gradient
                tmp = feGrad;
                tmp *= elemDat[idx].concentration(1);
                concentrationGrad_ += tmp;

                // phase viscosity
                viscosityAtIP_ += elemDat[idx].viscosity()*face().shapeValue[idx];

                //phase moledensity
                molarDensityAtIP_ += elemDat[idx].molarDensity()*face().shapeValue[idx];
            }
        }
        else {
            // use two-point gradients
            tmp = element.geometry().corner(face().i);
            tmp -= element.geometry().corner(face().j);
            Scalar dist = tmp.two_norm();

            tmp = face().normal;
            tmp /= face().normal.two_norm()*dist;

            potentialGrad_ = tmp;
            potentialGrad_ *= vVars_j.pressure() - vVars_i.pressure();
            concentrationGrad_ = tmp;
            concentrationGrad_ *= vVars_j.moleFrac(1) - vVars_i.moleFrac(1);
        }

        // correct the pressure by the hydrostatic pressure due to
        // gravity
        if (GET_PROP_VALUE(TypeTag, PTAG(EnableGravity))) {
            tmp = problem.gravity();
            tmp *= 0.5*(vVars_i.density() + vVars_j.density());
            potentialGrad_ -= tmp;
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
            = 1./2*(vDat_i.porosity() * vDat_i.tortuosity() * vDat_i.diffCoeff() +
                    vDat_j.porosity() * vDat_j.tortuosity() * vDat_j.diffCoeff());
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

    //! the effective diffusion coefficent in the porous medium
    Scalar diffCoeffPM_;

    //! the dispersion tensor in the porous medium
    Tensor dispersionTensor_;

    //! the intrinsic permeability tensor
    Tensor K_;

    //! viscosity of the fluid at the integration point
    Scalar viscosityAtIP_;

    //! molar densities of the fluid at the integration point
    Scalar molarDensityAtIP_;
};

} // end namepace

#endif
