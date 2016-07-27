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
 *        all fluxes of fluid phases over a face of a finite volume.
 *
 * This means pressure and mole-fraction gradients, phase densities at
 * the integration point, etc.
 *
 */
#ifndef DUMUX_1P2C_FLUX_VARIABLES_HH
#define DUMUX_1P2C_FLUX_VARIABLES_HH

#include <dune/common/exceptions.hh>
#include <dune/common/deprecated.hh>

#include "properties.hh"

#include <dumux/common/math.hh>
#include <dumux/common/valgrind.hh>

namespace Dumux
{

/*!
 * \ingroup OnePTwoCModel
 * \ingroup ImplicitFluxVariables
 * \brief This template class contains the data which is required to
 *        calculate the fluxes of the fluid phases over a face of a
 *        finite volume for the one-phase, two-component model.
 *
 * This means pressure and mole-fraction gradients, phase densities at
 * the intergration point, etc.
 */
template <class TypeTag>
class OnePTwoCFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename GET_PROP_TYPE(TypeTag, EffectiveDiffusivityModel) EffectiveDiffusivityModel;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        transportCompIdx = Indices::transportCompIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dim> DimVector;
    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimWorldMatrix;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

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
     * \todo The fvGeometry should be better initialized, passed and stored as an std::shared_ptr
     */
    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const int fIdx,
                const ElementVolumeVariables &elemVolVars,
                const bool onBoundary = false)
    {
        fvGeometryPtr_ = &fvGeometry;
        onBoundary_ = onBoundary;
        faceIdx_ = fIdx;

        mobilityUpwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MobilityUpwindWeight);

        viscosity_ = 0.0;
        molarDensity_ = 0.0;
        density_ = 0.0;
        potentialGrad_ = 0.0;
        moleFractionGrad_ = 0.0;

        asImp_().calculateGradients_(problem, element, elemVolVars);
        asImp_().calculateK_(problem, element, elemVolVars);
        asImp_().calculateVelocities_(problem, element, elemVolVars);
        asImp_().calculatePorousDiffCoeff_(problem, element, elemVolVars);
        asImp_().calculateDispersionTensor_(problem, element, elemVolVars);
    }

    /*!
    * \brief Return the pressure potential multiplied with the
    *        intrinsic permeability  and the face normal which
    *        goes from vertex i to vertex j.
    *
    * Note that the length of the face's normal is the area of the
    * phase, so this is not the actual velocity but the integral of
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
   GlobalPosition Kmvp() const
   { return Kmvp_; }

   /*!
    * \brief The face of the current sub-control volume. This may be either
    *        an inner sub-control-volume face or a SCV face on the boundary.
    */
   const SCVFace &face() const
   {
       if (onBoundary_)
           return fvGeometry_().boundaryFace[faceIdx_];
       else
           return fvGeometry_().subContVolFace[faceIdx_];
   }

    /*!
     * \brief Return the intrinsic permeability tensor \f$\mathrm{[m^2]}\f$.
     */
    const DimWorldMatrix &intrinsicPermeability() const
    { return K_; }

    /*!
     * \brief Return the dispersion tensor \f$\mathrm{[m^2/s]}\f$.
     */
    const DimWorldMatrix &dispersionTensor() const
    { return dispersionTensor_; }

    /*!
     * \brief Return the pressure potential gradient \f$\mathrm{[Pa/m]}\f$.
     */
    const GlobalPosition &potentialGrad() const
    { return potentialGrad_; }


    /*!
     * \brief Return the mole-fraction gradient of a component in a phase \f$\mathrm{[mol/mol/m)]}\f$.
     *
     * \param compIdx The index of the considered component
     */
    const GlobalPosition &moleFractionGrad(int compIdx) const
    {
       if (compIdx != 1)
       { DUNE_THROW(Dune::InvalidStateException,
                "The 1p2c model is supposed to need "
                "only the concentration gradient of "
                "the second component!"); }
       return moleFractionGrad_;
    };

    /*!
    * \brief The binary diffusion coefficient for each fluid phase in the porous medium \f$\mathrm{[m^2/s]}\f$.
    */
    Scalar porousDiffCoeff() const
    {
        // TODO: tensorial porousDiffCoeff_usion coefficients
        return porousDiffCoeff_;
    };

    /*!
    * \brief Return viscosity \f$\mathrm{[Pa s]}\f$ of a phase at the integration
    *        point.
    */
    Scalar viscosity() const
    { return viscosity_;}

    /*!
     * \brief Return molar density \f$\mathrm{[mol/m^3]}\f$ of a phase at the integration
     *        point.
     */
    Scalar molarDensity() const
    { return molarDensity_; }

    /*!
     * \brief Return density \f$\mathrm{[kg/m^3]}\f$ of a phase at the integration
     *        point.
     */
    Scalar density() const
    { return density_; }

    /*!
     * \brief Given the intrinsic permeability times the pressure
     *        potential gradient and SCV face normal for a phase,
     *        return the local index of the upstream control volume
     *        for a given phase.
     *
     *        \param normalFlux The flux over a face of the sub-control volume
     */
    int upstreamIdx(Scalar normalFlux) const
    { return (normalFlux >= 0)?face().i:face().j; }

    /*!
     * \brief Given the intrinsic permeability times the pressure
     *        potential gradient and SCV face normal for a phase,
     *        return the local index of the downstream control volume
     *        for a given phase.
     *
     *        \param normalFlux The flux over a face of the sub-control volume
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

    /*!
    * \brief Return the local index of the upstream control volume
    *        for a given phase.
    */
    int upstreamIdx(int phaseIdx) const
    { return upstreamIdx_; }

    /*!
     * \brief Return the local index of the downstream control volume
     *        for a given phase.
     */
    int downstreamIdx(int phaseIdx) const
    { return downstreamIdx_; }

    /*!
     * \brief Return the volumetric flux over a face of a given phase.
     *
     *        This is the calculated velocity multiplied by the unit normal
     *        and the area of the face.
     *        face().normal
     *        has already the magnitude of the area.
     *
     * \param phaseIdx index of the phase
     */
    Scalar volumeFlux(const unsigned int phaseIdx) const
    {
        assert (phaseIdx == Indices::phaseIdx);
        return volumeFlux_;
    }

protected:
    //! Returns the implementation of the flux variables (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    /*!
     * \brief Calculation of the pressure and mole-/mass-fraction gradients.
     *
     *        \param problem The considered problem file
     *        \param element The considered element of the grid
     *        \param elemVolVars The parameters stored in the considered element
     */
    void calculateGradients_(const Problem &problem,
                             const Element &element,
                             const ElementVolumeVariables &elemVolVars)
    {
        // loop over flux approximation points
        for (unsigned int idx = 0; idx < face().numFap; idx++)
        {
            // FE gradient at vertex idx
            const GlobalPosition &feGrad = face().grad[idx];

            // index for the element volume variables
            int volVarsIdx = face().fapIndices[idx];

            // the pressure gradient
            GlobalPosition tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].pressure();
            potentialGrad_ += tmp;

            // the mole-fraction gradient
            tmp = feGrad;
            tmp *= elemVolVars[volVarsIdx].moleFraction(transportCompIdx);
            moleFractionGrad_ += tmp;

            // phase viscosity
            viscosity_ += elemVolVars[volVarsIdx].viscosity()*face().shapeValue[idx];

            //phase molar density
            molarDensity_ += elemVolVars[volVarsIdx].molarDensity()*face().shapeValue[idx];

            //phase density
            density_ += elemVolVars[volVarsIdx].density()*face().shapeValue[idx];
        }


        ///////////////
        // correct the pressure gradients by the gravitational acceleration
        ///////////////
        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity)) {
            // calculate the phase density at the integration point. we
            // only do this if the wetting phase is present in both cells
            Scalar rhoI = elemVolVars[face().i].density();
            Scalar rhoJ = elemVolVars[face().j].density();
            Scalar density = (rhoI + rhoJ)/2;

            // ask for the gravitational acceleration at the given SCV face
            GlobalPosition g(problem.gravityAtPos(face().ipGlobal));

            // make it a force
            g *= density;

            // calculate the final potential gradient
            potentialGrad_ -= g;
        }
    }

    /*!
    * \brief Calculation of the harmonic mean of the intrinsic permeability
    *        uses the meanK function in the boxspatialparameters.hh file in the folder
    *        material/spatialparameters
    *
    *        \param problem The considered problem file
    *        \param element The considered element of the grid
    *        \param elemVolVars The parameters stored in the considered element
    */
    void calculateK_(const Problem &problem,
                     const Element &element,
                     const ElementVolumeVariables &elemVolVars)
    {
        const SpatialParams &sp = problem.spatialParams();
        if (GET_PROP_VALUE(TypeTag, ImplicitIsBox))
        {
            sp.meanK(K_,
                     sp.intrinsicPermeability(element,
                                              fvGeometry_(),
                                              face().i),
                     sp.intrinsicPermeability(element,
                                              fvGeometry_(),
                                              face().j));
        }
        else
        {
            const Element& elementI = fvGeometry_().neighbors[face().i];
            FVElementGeometry fvGeometryI;
            fvGeometryI.subContVol[0].global = elementI.geometry().center();

            const Element& elementJ = fvGeometry_().neighbors[face().j];
            FVElementGeometry fvGeometryJ;
            fvGeometryJ.subContVol[0].global = elementJ.geometry().center();

            sp.meanK(K_,
                     sp.intrinsicPermeability(elementI, fvGeometryI, 0),
                     sp.intrinsicPermeability(elementJ, fvGeometryJ, 0));
        }
    }

    /*!
      * \brief Calculation of the velocity normal to face using Darcy's law.
      *     Tensorial permeability is multiplied with the potential gradient and the face normal.
      *     Identify upstream node of face.
      *
      *        \param problem The considered problem file
      *        \param element The considered element of the grid
      *        \param elemVolVars The parameters stored in the considered element
      */
    void calculateVelocities_(const Problem &problem,
                              const Element &element,
                              const ElementVolumeVariables &elemVolVars)
    {
        K_.mv(potentialGrad_, Kmvp_);
        KmvpNormal_ = -(Kmvp_*face().normal);

        // set the upstream and downstream vertices
        upstreamIdx_ = face().i;
        downstreamIdx_ = face().j;

        if (KmvpNormal_ < 0)
        {
            std::swap(upstreamIdx_,
                      downstreamIdx_);
        }

        volumeFlux_ = KmvpNormal_;
        volumeFlux_ *= mobilityUpwindWeight_/elemVolVars[upstreamIdx_].viscosity()
                    + (1.0 - mobilityUpwindWeight_)/elemVolVars[downstreamIdx_].viscosity();
    }
    /*!
    * \brief Calculation of the effective diffusion coefficient.
    *
    *        \param problem The considered problem file
    *        \param element The considered element of the grid
    *        \param elemVolVars The parameters stored in the considered element
    */
    void calculatePorousDiffCoeff_(const Problem &problem,
                                   const Element &element,
                                   const ElementVolumeVariables &elemVolVars)
    {
        const VolumeVariables &volVarsI = elemVolVars[face().i];
        const VolumeVariables &volVarsJ = elemVolVars[face().j];

        const Scalar diffCoeffI = EffectiveDiffusivityModel::effectiveDiffusivity(volVarsI.porosity(),
                                                                     /*sat=*/1.0,
                                                                     volVarsI.diffCoeff());

        const Scalar diffCoeffJ = EffectiveDiffusivityModel::effectiveDiffusivity(volVarsJ.porosity(),
                                                                     /*sat=*/1.0,
                                                                     volVarsJ.diffCoeff());

        // -> harmonic mean
        porousDiffCoeff_ = harmonicMean(diffCoeffI, diffCoeffJ);
    }

    /*!
    * \brief Calculation of the dispersion.
    *
    *        \param problem The considered problem file
    *        \param element The considered element of the grid
    *        \param elemVolVars The parameters stored in the considered element
    */
    void calculateDispersionTensor_(const Problem &problem,
                                    const Element &element,
                                    const ElementVolumeVariables &elemVolVars)
    {
        const VolumeVariables &volVarsI = elemVolVars[face().i];
        const VolumeVariables &volVarsJ = elemVolVars[face().j];

        //calculate dispersivity at the interface: [0]: alphaL = longitudinal disp. [m], [1] alphaT = transverse disp. [m]
        Scalar dispersivity[2];
        dispersivity[0] = 0.5 * (volVarsI.dispersivity()[0] +  volVarsJ.dispersivity()[0]);
        dispersivity[1] = 0.5 * (volVarsI.dispersivity()[1] +  volVarsJ.dispersivity()[1]);

        //calculate velocity at interface: v = -1/mu * vDarcy = -1/mu * K * grad(p)
        GlobalPosition velocity;
        Valgrind::CheckDefined(potentialGrad());
        Valgrind::CheckDefined(K_);
        K_.mv(potentialGrad(), velocity);
        velocity /= - 0.5 * (volVarsI.viscosity() + volVarsJ.viscosity());

        //matrix multiplication of the velocity at the interface: vv^T
        dispersionTensor_ = 0;
        for (int i=0; i<dim; i++)
            for (int j = 0; j<dim; j++)
                dispersionTensor_[i][j] = velocity[i]*velocity[j];

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

    // return const reference to fvGeometry
    const FVElementGeometry& fvGeometry_() const
    { return *fvGeometryPtr_; }

    int faceIdx_;
    bool onBoundary_;

    //! pressure potential gradient
    GlobalPosition potentialGrad_;
    //! mole-fraction gradient
    GlobalPosition moleFractionGrad_;
    //! the effective diffusion coefficent in the porous medium
    Scalar porousDiffCoeff_;

    //! the dispersion tensor in the porous medium
    DimWorldMatrix dispersionTensor_;

    //! the intrinsic permeability tensor
    DimWorldMatrix K_;
    // intrinsic permeability times pressure potential gradient
    GlobalPosition Kmvp_;
    // projected on the face normal
    Scalar KmvpNormal_;

    // local index of the upwind vertex for each phase
    int upstreamIdx_;
    // local index of the downwind vertex for each phase
    int downstreamIdx_;

    //! viscosity of the fluid at the integration point
    Scalar viscosity_;

    //! molar densities of the fluid at the integration point
    Scalar molarDensity_, density_;

    Scalar volumeFlux_; //!< Velocity multiplied with normal (magnitude=area)
    Scalar mobilityUpwindWeight_; //!< Upwind weight for mobility. Set to one for full upstream weighting

private:
    const FVElementGeometry* fvGeometryPtr_; //!< Information about the geometry of discretization
};

} // end namespace

#endif
