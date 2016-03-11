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
 *
 * \brief Element-wise calculation of the residual for the finite volume in the
 *        two phase discrete fracture-matrix model.
 */
#ifndef DUMUX_MODELS_2PDFM_LOCAL_RESIDUAL_HH
#define DUMUX_MODELS_2PDFM_LOCAL_RESIDUAL_HH

#include <dumux/porousmediumflow/2p/implicit/localresidual.hh>
#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup TwoPDFMModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the two-phase discrete fracture fully implicit model.
 */
template<class TypeTag>
class TwoPDFMLocalResidual : public TwoPLocalResidual<TypeTag>
{
protected:
    typedef TwoPLocalResidual<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum
    {
        contiWEqIdx = Indices::contiWEqIdx,
        contiNEqIdx = Indices::contiNEqIdx,
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) GridType;
    typedef typename GridType::ctype DT;
    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename Dune::ReferenceElements<DT, dim> ReferenceElements;
    typedef typename Dune::ReferenceElement<DT, dim> ReferenceElement;

public:
    /*!
     * \brief Constructor. Sets the upwind weight.
     */
    TwoPDFMLocalResidual()
    {
        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        massUpwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);
    };

    /*!
     * \brief Evaluate the amount all conservation quantities
     *        (e.g. phase mass) within a finite sub-control volume.
     *
     *  \param storage The phase mass within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    void computeStorage(PrimaryVariables &storage, int scvIdx, bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit Euler method.
        ParentType::computeStorage(storage,scvIdx,usePrevSol);
        computeStorageFracture(storage,scvIdx,usePrevSol);
    }

    /*!
     * \brief Evaluate the storage term of the current solution in a
     *        lower-dimensional fracture.
     *
     *  \param storage The phase mass within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    void computeStorageFracture(PrimaryVariables &storage, int scvIdx, bool usePrevSol) const
    {
        const ElementVolumeVariables &elemVolVars = usePrevSol ? this->prevVolVars_()
                        : this->curVolVars_();
        const VolumeVariables &volVars = elemVolVars[scvIdx];
        const Element &elem = this->element_();
        bool isFracture = this->problem_().spatialParams().isVertexFracture(elem, scvIdx);
        /*
         * sn = wf * SnF + wm * SnM
         * First simple case before determining the real wf is to assume that it is 0
         * and wm = 1
         *
         */
        ///////////////////////////////////////////////////////////////////////
        Scalar wf, wm; //volumetric fractions of fracture and matrix;
        Scalar fractureVolume = 0.0;
        wf = 0.0;
        /*
         * Calculate the fracture volume fraction wf = 0.5 * Fwidth * 0.5 * Length
         */
        const auto geometry = elem.geometry();
        Dune::GeometryType geomType = geometry.type();
        const ReferenceElement &refElement = ReferenceElements::general(geomType);

        Scalar vol; //subcontrol volume
        FVElementGeometry fvelem = this->fvGeometry_();
        vol = fvelem.subContVol[scvIdx].volume;
        for (int fIdx=0; fIdx<refElement.size(1); fIdx++)
        {
            SCVFace face = fvelem.subContVolFace[fIdx];
            int i=face.i;
            int j=face.j;

            if (this->problem_().spatialParams().isEdgeFracture(elem, fIdx)
                    && (i == scvIdx || j == scvIdx))
            {
                Scalar fracture_width = this->problem_().spatialParams().fractureWidth(elem, fIdx);

                const GlobalPosition global_i = geometry.corner(i);
                const GlobalPosition global_j = geometry.corner(j);
                GlobalPosition diff_ij = global_j;
                diff_ij -=global_i;
                Scalar fracture_length = 0.5*diff_ij.two_norm();
                //half is taken from the other element
                fractureVolume += 0.5 * fracture_length * fracture_width;
            }
        }
        wf = fractureVolume/vol;
        wm = 1-wf;
        ///////////////////////////////////////////////////////////////////////
        Scalar storageFracture[numPhases];
        Scalar storageMatrix [numPhases];
        storageFracture[wPhaseIdx]    = 0.0;
        storageFracture[nPhaseIdx]    = 0.0;
        storageMatrix[wPhaseIdx]    = 0.0;
        storageMatrix[nPhaseIdx]    = 0.0;
        //        const GlobalPosition &globalPos = geometry.corner(scvIdx);

        if (isFracture)
        {
            for (int phaseIdx = 0; phaseIdx<2; phaseIdx++)
            {
                storageFracture[phaseIdx] = volVars.density(phaseIdx)
                                        * volVars.porosityFracture()
                                        * wf
                                        * volVars.saturationFracture(phaseIdx);
                storageMatrix[phaseIdx] = volVars.density(phaseIdx)
                                        * volVars.porosity()
                                        * wm
                                        * volVars.saturationMatrix(phaseIdx);
            }
        }
        else
        {
            for (int phaseIdx = 0; phaseIdx < 2; phaseIdx++)
            {
                storageFracture[phaseIdx] = 0.0;
                storageMatrix[phaseIdx] = volVars.density(phaseIdx)
                                        * volVars.porosity()
                                        * volVars.saturation(phaseIdx);
            }

        }
        // wetting phase mass
        storage[contiWEqIdx] =  storageFracture[wPhaseIdx]
                            + storageMatrix[wPhaseIdx];

        // non-wetting phase mass
        storage[contiNEqIdx] =  storageFracture[nPhaseIdx]
                            + storageMatrix[nPhaseIdx];
    }

    /*!
     * \brief Evaluates the mass flux over a face of a sub-control
     *        volume.
     *
     * \param flux The flux over the SCV (sub-control-volume) face for each phase
     * \param fIdx The index of the SCV face
     * \param onBoundary A boolean variable to specify whether the flux variables
     *        are calculated for interior SCV faces or boundary faces, default=false
     */
    void computeFlux(PrimaryVariables &flux, int fIdx, const bool onBoundary=false) const
    {
        FluxVariables fluxVars;
        fluxVars.update(this->problem_(),
                        this->element_(),
                        this->fvGeometry_(),
                        fIdx,
                        this->curVolVars_(),
                        onBoundary);
        flux = 0;
        asImp_()->computeAdvectiveFlux(flux, fluxVars);
        asImp_()->computeAdvectiveFluxFracture(flux, fluxVars);
        asImp_()->computeDiffusiveFlux(flux, fluxVars);
        asImp_()->computeDiffusiveFluxFracture(flux, fluxVars);
    }

    /*!
     * \brief Adds the flux vector in the lower dimensional fracture to the
     *        flux vector over the face of a sub-control volume.
     *
     * This method is called by compute flux and is mainly there for
     * derived models to ease adding equations selectively.
     *
     * \param flux The advective flux over the sub-control-volume face for each phase
     * \param fluxVars The flux variables at the current SCV
     */
    void computeAdvectiveFluxFracture(PrimaryVariables &flux, const FluxVariables &fluxVars) const
    {
        ////////
        // advective fluxes of all components in all phases
        ////////
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // calculate the flux in direction of the
            // current fracture
            //
            // v = - (K grad p) * n
            //
            // (the minus comes from the Darcy law which states that
            // the flux is from high to low pressure potentials.)

            // data attached to upstream and the downstream vertices
            // of the current phase
            const VolumeVariables &upFracture =
                    this->curVolVars_(fluxVars.upstreamFractureIdx[phaseIdx]);
            const VolumeVariables &dnFracture =
                    this->curVolVars_(fluxVars.downstreamFractureIdx[phaseIdx]);

            Scalar fractureFlux =
                          0.5 * fluxVars.vDarcyFracture_[phaseIdx]
                              * ( massUpwindWeight_ // upstream vertex
                                 *  (upFracture.density(phaseIdx)
                                    * upFracture.mobilityFracture(phaseIdx))
                                 + (1 - massUpwindWeight_) // downstream vertex
                                    * (dnFracture.density(phaseIdx)
                                      * dnFracture.mobilityFracture(phaseIdx)));

            flux[phaseIdx] += fractureFlux;
        }
    }

    /*!
     * \brief Adds the diffusive flux of the fracture to the flux vector over
     *        the face of a sub-control volume.
     *
     * \param flux The diffusive flux over the sub-control-volume face for each phase
     * \param fluxVars The flux variables at the current SCV
     *
     * It doesn't do anything in two-phase model but is used by the
     * non-isothermal two-phase models to calculate diffusive heat
     * fluxes
     */
    void computeDiffusiveFluxFracture(PrimaryVariables &flux, const FluxVariables &fluxVars) const
    {
        // diffusive fluxes
        flux += 0.0;
    }

protected:
    Implementation *asImp_()
    {
        return static_cast<Implementation *> (this);
    }
    const Implementation *asImp_() const
    {
        return static_cast<const Implementation *> (this);
    }

private:
    Scalar massUpwindWeight_;
};

} // end namespace

#endif // DUMUX_MODELS_2PDFM_LOCAL_RESIDUAL_HH
