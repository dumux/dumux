// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
#include "1p2cproperties.hh"

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
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, OnePTwoCIndices) Indices;
    enum {
        phaseIdx = Indices::phaseIdx,
        comp0Idx = Indices::comp0Idx,
        comp1Idx = Indices::comp1Idx
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    enum { dim = GridView::dimension };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dim> Vector;
    typedef Dune::FieldMatrix<Scalar, dim, dim> Tensor;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::BoundaryFace BoundaryFace;

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

            pressureAtBIP_ = Scalar(0);
            massFractionAtBIP_ = Scalar(0);
            densityAtIP_ = Scalar(0);
            viscosityAtIP_ = Scalar(0);
            molarDensityAtIP_ = Scalar(0);
            potentialGrad_ = Scalar(0);
            concentrationGrad_ = Scalar(0);
            moleFracGrad_ = Scalar(0);
            massFracGrad_ = Scalar(0);
            K_= Scalar(0);

        calculateBoundaryValues_(problem, element, elemDat);
    };

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
     * \brief The binary diffusion coefficient for each fluid phase in the porous medium.
     */
    Scalar porousDiffCoeff() const
    { return porousDiffCoeff_; };

    /*!
     * \brief Return pressure \f$\mathrm{[Pa]}\f$ of a phase at the integration
     *        point.
     */
    Scalar pressureAtBIP() const
    { return pressureAtBIP_; }

    /*!
     * \brief Return massFraction \f$\mathrm{[-]}\f$ of component 1 at the integration
     *        point.
     */
    Scalar massFractionAtBIP() const
    { return massFractionAtBIP_; }

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

    const FVElementGeometry &fvElemGeom() const
    { return fvElemGeom_; }

    const BoundaryFace& boundaryFace() const
    { return *boundaryFace_; }

protected:

    /*!
    * \brief Transforms scalar permeability into tensor.
    *
    *        \param k scalar permeability
    *
    */
    void calculateK_(Scalar k)
    {
        K_[0][0] = K_[1][1] = k;
    }

    void calculateK_(Tensor K)
    {
       K_ = K;
    }
    /*!
   * \brief Calculation of the pressure and concentration gradients
   *
   *        \param problem The considered problem file
   *        \param element The considered element of the grid
   *        \param elemDat The parameters stored in the considered element
   */
    void calculateBoundaryValues_(const Problem &problem,
                             const Element &element,
                             const ElementVolumeVariables &elemDat)
    {
        Vector tmp(0.0);

        // calculate gradients and secondary variables at IPs of the boundary
        for (int idx = 0;
             idx < fvElemGeom_.numVertices;
             idx++) // loop over adjacent vertices
        {
            // FE gradient at vertex idx
            const Vector& feGrad = boundaryFace_->grad[idx];

            pressureAtBIP_ += elemDat[idx].pressure() *
                             boundaryFace_->shapeValue[idx];

            // compute sum of pressure gradients for each phase
                // the pressure gradient
                tmp = feGrad;
                tmp *= elemDat[idx].pressure();
                potentialGrad_ += tmp;

            massFractionAtBIP_ += 
                elemDat[idx].massFraction(comp1Idx)
                * boundaryFace_->shapeValue[idx];
            // the concentration gradient of the non-wetting
            // component in the wetting phase
            tmp = feGrad;
            tmp *= elemDat[idx].fluidState().molarity(phaseIdx, comp1Idx);
            concentrationGrad_ += tmp;

            tmp = feGrad;
            tmp *= elemDat[idx].fluidState().moleFraction(phaseIdx, comp1Idx);
            moleFracGrad_ += tmp;

            tmp = feGrad;
            tmp *= elemDat[idx].fluidState().massFraction(phaseIdx, comp1Idx);
            massFracGrad_ += tmp;


            densityAtIP_ += elemDat[idx].density()*boundaryFace_->shapeValue[idx];
            viscosityAtIP_ += elemDat[idx].viscosity()*boundaryFace_->shapeValue[idx];
            molarDensityAtIP_ += elemDat[idx].molarDensity()*boundaryFace_->shapeValue[idx];
        }
        
        tmp = problem.gravity();
        tmp *= densityAtIP_;

        potentialGrad_ -= tmp;
        
        calculateK_(problem.spatialParameters().intrinsicPermeability(element, fvElemGeom_, scvIdx_));
        Vector Kmvp;
        K_.mv(potentialGrad_, Kmvp);
        KmvpNormal_ = 0;
        for (int i = 0; i < dim; ++i)
            KmvpNormal_ += Kmvp[i]*boundaryFace_->normal[i];
        KmvpNormal_ *= -1;
        
        const VolumeVariables &vertDat = elemDat[scvIdx_];
        
        Scalar tau = problem.spatialParameters().tortuosity(element, fvElemGeom_, scvIdx_);
        
        porousDiffCoeff_ = vertDat.porosity()*tau*vertDat.diffCoeff();
    }
    
    const FVElementGeometry &fvElemGeom_;
    const BoundaryFace *boundaryFace_;
    
    // gradients
    Vector potentialGrad_;
    Vector concentrationGrad_;
    Vector moleFracGrad_;
    Vector massFracGrad_;

    // quanitities at the integration point
    Scalar pressureAtBIP_;
    Scalar massFractionAtBIP_;
    Scalar densityAtIP_;
    Scalar viscosityAtIP_;
    Scalar molarDensityAtIP_;

    //intrinsic permeability tensor
    Tensor K_;
    // intrinsic permeability times pressure potential gradient
    // projected on the face normal
    Scalar KmvpNormal_;

    // the diffusion coefficient for the porous medium
    Scalar porousDiffCoeff_;

    int scvIdx_;
};

} // end namespace

#endif
