// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
 */
#ifndef DUMUX_2P_DFM_FLUX_VARIABLES_HH
#define DUMUX_2P_DFM_FLUX_VARIABLES_HH

#include "2pdfmproperties.hh"
#include <dumux/boxmodels/2p/2pfluxvariables.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/math.hh>

namespace Dumux
{

/*!
 * \ingroup TwoPBoxModel
 * \ingroup BoxFluxVariables
 * \brief This template class contains the data which is required to
 *        calculate the fluxes of the fluid phases over a face of a
 *        finite volume for the two-phase model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the intergration point, etc.
 */
template <class TypeTag>
class TwoPDFMFluxVariables : public TwoPFluxVariables<TypeTag>
{
	typedef TwoPFluxVariables<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) SpatialParameters;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> Tensor;
    typedef Dune::FieldVector<Scalar, dimWorld> Vector;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;

public:
    /*
     * \brief The constructor
     *
     * \param problem The problem
     * \param element The finite element
     * \param elemGeom The finite-volume geometry in the box scheme
     * \param faceIdx The local index of the SCV (sub-control-volume) face
     * \param elemVolVars The volume variables of the current element
     * \param onBoundary A boolean variable to specify whether the flux variables
     * are calculated for interior SCV faces or boundary faces, default=false
     */
    TwoPDFMFluxVariables(const Problem &problem,
                 const Element &element,
                 const FVElementGeometry &elemGeom,
                 int faceIdx,
                 const ElementVolumeVariables &elemVolVars,
                 const bool onBoundary = false)
        : ParentType (problem, element, elemGeom, faceIdx, elemVolVars, onBoundary)
//        	fvElemGeom_(elemGeom), onBoundary_(onBoundary)
    {
        scvfIdx_ = faceIdx;
        faceSCV_ = &this->fvElemGeom_.subContVolFace[faceIdx];
        for (int phase = 0; phase < numPhases; ++phase) {
        	potentialGradFracture_[phase] = 0.0;
        }

        calculateGradientsInFractures_(problem, element, elemVolVars, faceIdx);
        calculateVelocitiesFracture_(problem, element, elemVolVars, faceIdx);
//        calculateK_(problem, element, elemVolVars);
    }

//public:
//    /*
//     * \brief Return the intrinsic permeability.
//     */
//    const Tensor &intrinsicPermeability() const
//    { return K_; }
//
//    /*!
//     * \brief Return the pressure potential gradient.
//     *
//     * \param phaseIdx The index of the fluid phase
//     */
//    const Vector &potentialGrad(int phaseIdx) const
//    { return potentialGrad_[phaseIdx]; }
//
//    /*!
//     * \brief Return the local index of the downstream control volume
//     *        for a given phase as a function of the normal flux.
//     *
//     * \param normalFlux The normal flux i.e. the given intrinsic permeability
//     *                   times the pressure potential gradient and SCV face normal.
//     */
//    int downstreamIdx(Scalar normalFlux) const
//    { return (normalFlux >= 0)?face().j:face().i; }
//
//    /*!
//     * \brief Return the local index of the upstream control volume
//     *        for a given phase as a function of the normal flux.
//     *
//     * \param normalFlux The normal flux i.e. the given intrinsic permeability
//     *                   times the pressure potential gradient and SCV face normal.
//     */
//    int upstreamIdx(Scalar normalFlux) const
//    { return (normalFlux > 0)?face().i:face().j; }
//
//    /*!
//     * \brief Return the SCV (sub-control-volume) face
//    */
//    const SCVFace &face() const
//    {
//        if (this->onBoundary_)
//            return this->fvElemGeom_.boundaryFace[scvfIdx_];
//        else
//            return this->fvElemGeom_.subContVolFace[scvfIdx_];
//    }
//

//
//private:
    void calculateGradientsInFractures_(const Problem &problem,
                             const Element &element,
                             const ElementVolumeVariables &elemVolVars,
                             int faceIdx)
    {
        // calculate gradients
        for (int idx = 0;
             idx < this->fvElemGeom_.numFAP;
             idx++) // loop over adjacent vertices
        {
            // FE gradient at vertex idx
            const Vector &feGrad = this->face().grad[idx];

            int i = faceSCV_->i;
            int j = faceSCV_->j;


			// index for the element volume variables
			int volVarsIdx = this->face().fapIndices[idx];

            // compute sum of pressure gradients for each phase
            for (int phase = 0; phase < numPhases; phase++)
            {
                // the pressure gradient
//                Vector tmp(feGrad);
//                tmp *= elemVolVars[volVarsIdx].pressure(phase);
//                potentialGrad_[phase] += tmp;

                isFracture_ = problem.spatialParameters().isEdgeFracture(element, faceIdx);
                if (isFracture_)
                {
                    const GlobalPosition localIdx_i = element.geometry().corner(i);
					const GlobalPosition localIdx_j = element.geometry().corner(j);
					GlobalPosition diff_ij = localIdx_j;
					diff_ij -=localIdx_i;
					//pressFracture = pressMatrix
					potentialGradFracture_[phase] =
					(elemVolVars[j].pressure(phase) - elemVolVars[i].pressure(phase))
					/
					diff_ij.two_norm();
//					std::cout<<"pressureGradFracture["<<phase<<"]"
//							<<potentialGradFracture_[phase]<<"\n";//TODO delete
                }
                else
                {
                	potentialGradFracture_[phase] = 0;
                }

            }
        }
      }
//
//    void calculateK_(const Problem &problem,
//                     const Element &element,
//                     const ElementVolumeVariables &elemVolVars)
//    {
//        const SpatialParameters &spatialParams = problem.spatialParameters();
//        // calculate the intrinsic permeability
//        spatialParams.meanK(K_,
//                            spatialParams.intrinsicPermeability(element,
//                                                                fvElemGeom_,
//                                                                face().i),
//                            spatialParams.intrinsicPermeability(element,
//                                                                fvElemGeom_,
//                                                                face().j));
//    }


    void calculateVelocitiesFracture_(const Problem &problem,
									  const Element &element,
									  const ElementVolumeVariables &elemVolVars,
									  int faceIdx)
    {
        isFracture_ = problem.spatialParameters().isEdgeFracture(element, faceIdx);
        fractureWidth_ = problem.spatialParameters().fractureWidth(element, faceIdx);

        Scalar KFracture, KFi, KFj;
        if (isFracture_)
        {
            KFi = problem.spatialParameters().
                    intrinsicPermeabilityFracture(element, this->fvElemGeom_, faceSCV_->i);
            KFj = problem.spatialParameters().
                    intrinsicPermeabilityFracture(element, this->fvElemGeom_, faceSCV_->j);
        }
        else
        {
            KFi = 0;
            KFj = 0;
        }

        KFracture = Dumux::harmonicMean(KFi, KFj);

        // temporary vector for the Darcy velocity
        Scalar vDarcyFracture;
        for (int phase=0; phase < numPhases; phase++)
        {
        	vDarcyFracture = - (KFracture*potentialGradFracture_[phase]);
        	Valgrind::CheckDefined(KFracture);

            Valgrind::CheckDefined(fractureWidth_);
            vDarcyFracture_[phase] = (vDarcyFracture * fractureWidth_);
            Valgrind::CheckDefined(vDarcyFracture_[phase]);
            if (vDarcyFracture_[phase]!=0){
//            std::cout<<"fractureWidth "<<fractureWidth_ <<std::endl;//TODO delete
//            std::cout<<"vDarcyFracture["<<phase<<"]="<<vDarcyFracture_[phase]<<std::endl;//TODO delete
            }

        }

        // set the upstream and downstream vertices
        for (int phase = 0; phase < numPhases; ++phase)
        {
            upstreamFractureIdx[phase] = faceSCV_->i;
            downstreamFractureIdx[phase] = faceSCV_->j;

            if (vDarcyFracture_[phase] < 0) {
            	std::swap(upstreamFractureIdx[phase],
                          downstreamFractureIdx[phase]);
            }
        }
    }


public:
//    const FVElementGeometry &fvElemGeom_;
    int scvfIdx_;

    // gradients
    Vector potentialGrad_[numPhases];
    Scalar potentialGradFracture_[numPhases];
    PhasesVector vDarcyFracture_;

    int upstreamFractureIdx[numPhases];
    int downstreamFractureIdx[numPhases];

    bool isFracture_;
    Scalar fractureWidth_;
    const SCVFace *faceSCV_;


};

} // end namepace

#endif
