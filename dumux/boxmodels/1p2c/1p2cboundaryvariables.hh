// $Id$
/*****************************************************************************
 *   Copyright (C) 2010 by Katherina Baber, Klaus Mosthaf                    *
 *   Copyright (C) 2008-2009 by Bernd Flemisch, Andreas Lauser               *
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
 *        all fluxes of the fluid phase over the boundary of a finite volume.
 *
 * This means pressure and temperature gradients, phase densities at
 * the integration point of the boundary, etc.
 */
#ifndef DUMUX_1P2C_BOUNDARY_VARIABLES_HH
#define DUMUX_1P2C_BOUNDARY_VARIABLES_HH

#include <dumux/common/math.hh>

namespace Dumux
{

/*!
 * \ingroup OnePTwoCModel
 * \brief This template class contains the data which is required to
 *        calculate the fluxes of the fluid phases over the boundary of a
 *        finite volume for the 1p2c model.
 *
 * This means pressure and velocity gradients, phase density and viscosity at
 * the integration point of the boundary, etc.
 */
template <class TypeTag>
class OnePTwoCBoundaryVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePTwoCIndices)) Indices;

    enum {
        dim = GridView::dimension,
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dim> VelocityVector;
    typedef Dune::FieldVector<Scalar, dim> ScalarGradient;
    typedef Dune::FieldMatrix<Scalar, dim, dim> VectorGradient;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolume SCV;
    typedef typename FVElementGeometry::BoundaryFace BoundaryFace;

    enum {
        phaseIdx = Indices::phaseIdx,

        comp0Idx = Indices::comp0Idx,
        comp1Idx = Indices::comp1Idx,

        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases))
    };

public:
    OnePTwoCBoundaryVariables(const Problem &problem,
                     const Element &element,
                     const FVElementGeometry &elemGeom,
                     int boundaryFaceIdx,
                     const ElementVolumeVariables &elemDat,
                     int scvIdx)
        : fvElemGeom_(elemGeom), scvIdx_(scvIdx)
    {
        boundaryFace_ = &fvElemGeom_.boundaryFace[boundaryFaceIdx];

            densityAtIP_ = Scalar(0);
            viscosityAtIP_ = Scalar(0);
            pressureAtIP_ = Scalar(0);
            molarDensityAtIP_ = Scalar(0);
            potentialGrad_ = Scalar(0);
            concentrationGrad_ = Scalar(0);
            molarConcGrad_ = Scalar(0);
            K_= Scalar(0);

        calculateBoundaryValues_(problem, element, elemDat);
    };

    Scalar KmvpNormal() const
    { return KmvpNormal_; }

    /*!
     * \brief The binary diffusion coefficient for each fluid phase.
     */
    Scalar porousDiffCoeff() const
    { return porousDiffCoeff_; };

    /*!
     * \brief Return density \f$\mathrm{[kg/m^3]}\f$ of a phase at the integration
     *        point.
     */
    Scalar densityAtIP() const
    { return densityAtIP_; }

    /*!
    * \brief Return viscosity \f$\mathrm{[Pa s]}\f$ of a phase at the integration
    *        point.
    */
   Scalar viscosityAtIP() const
   { return viscosityAtIP_; }

    /*!
     * \brief Return molar density \f$\mathrm{[mol/m^3]}\f$ of a phase at the integration
     *        point.
     */
    Scalar molarDensityAtIP() const
    { return molarDensityAtIP_; }

    /*!
     * \brief The concentration gradient of a component in a phase.
     */
    const ScalarGradient &concentrationGrad() const
    { return concentrationGrad_; };

    /*!
     * \brief The molar concentration gradient of a component in a phase.
     */
    const ScalarGradient &molarConcGrad(int compIdx) const
    { return molarConcGrad_; };

    const FVElementGeometry &fvElemGeom() const
    { return fvElemGeom_; }

    const BoundaryFace& boundaryFace() const
    { return *boundaryFace_; }

protected:

    void calculateK_(Scalar k)
    {
        K_[0][0] = K_[1][1] = k;
    }

    void calculateK_(VectorGradient K)
    {
       K_ = K;
    }

    void calculateBoundaryValues_(const Problem &problem,
                             const Element &element,
                             const ElementVolumeVariables &elemDat)
    {
        ScalarGradient tmp(0.0);

        // calculate gradients and secondary variables at IPs of the boundary
        for (int idx = 0;
             idx < fvElemGeom_.numVertices;
             idx++) // loop over adjacent vertices
        {
            // FE gradient at vertex idx
            const ScalarGradient& feGrad = boundaryFace_->grad[idx];

            // compute sum of pressure gradients for each phase
                // the pressure gradient
                tmp = feGrad;
                tmp *= elemDat[idx].pressure();
                potentialGrad_ += tmp;

            // the concentration gradient of the non-wetting
            // component in the wetting phase
            tmp = feGrad;
            tmp *= elemDat[idx].fluidState().massFrac(phaseIdx, comp1Idx);
            concentrationGrad_ += tmp;

            tmp = feGrad;
            tmp *= elemDat[idx].fluidState().moleFrac(phaseIdx, comp1Idx);
            molarConcGrad_ += tmp;

                pressureAtIP_ += elemDat[idx].pressure()*boundaryFace_->shapeValue[idx];
                densityAtIP_ += elemDat[idx].density()*boundaryFace_->shapeValue[idx];
                viscosityAtIP_ += elemDat[idx].viscosity()*boundaryFace_->shapeValue[idx];
                molarDensityAtIP_ += elemDat[idx].molarDensity()*boundaryFace_->shapeValue[idx];
        }

            tmp = problem.gravity();
            tmp *= densityAtIP_;

            potentialGrad_ -= tmp;

            calculateK_(problem.spatialParameters().intrinsicPermeability(element, fvElemGeom_, scvIdx_));
            ScalarGradient Kmvp;
            K_.mv(potentialGrad_, Kmvp);
            KmvpNormal_ = - (Kmvp*boundaryFace_->normal);

            const VolumeVariables &vertDat = elemDat[scvIdx_];

            Scalar tau = problem.spatialParameters().tortuosity(element, fvElemGeom_, scvIdx_);

            porousDiffCoeff_ = vertDat.porosity()*tau*vertDat.diffCoeff();
    }

    const FVElementGeometry &fvElemGeom_;
    const BoundaryFace *boundaryFace_;

    // gradients
    ScalarGradient potentialGrad_;
    ScalarGradient concentrationGrad_;
    ScalarGradient molarConcGrad_;

    // density of each face at the integration point
    Scalar pressureAtIP_;
    Scalar densityAtIP_;
    Scalar viscosityAtIP_;
    Scalar molarDensityAtIP_;

    //intrinsic permeability tensor
    VectorGradient K_;
    // intrinsic permeability times pressure potential gradient
    // projected on the face normal
    Scalar KmvpNormal_;

    // the diffusion coefficient for the porous medium
    Scalar porousDiffCoeff_;

    int scvIdx_;
};

} // end namespace

#endif
