// $Id: 2plocalresidual.hh 3794 2010-06-25 16:04:52Z melanie $
/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
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
 * \brief Element-wise calculation of the residual for the two-phase box model.
 */
#ifndef DUMUX_TWOP_LOCAL_RESIDUAL_BASE_HH
#define DUMUX_TWOP_LOCAL_RESIDUAL_BASE_HH

#include <dumux/boxmodels/common/boxmodel.hh>

#include "2pproperties.hh"

namespace Dumux
{
/*!
 * \ingroup TwoPBoxModel
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the two-phase box model.
 *
 * This class is also used for the non-isothermal model, which means
 * that it uses static polymorphism.
 */
template<class TypeTag>
class TwoPLocalResidual: public BoxLocalResidual<TypeTag>
{
protected:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalResidual)) Implementation;
    typedef TwoPLocalResidual<TypeTag> ThisType;
    typedef BoxLocalResidual<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    enum
    {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),

        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,

        contiWEqIdx = Indices::contiWEqIdx,
        contiNEqIdx = Indices::contiNEqIdx,

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVariables)) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;

    typedef Dune::FieldVector<Scalar, dimWorld> Vector;
    typedef Dune::FieldMatrix<Scalar, dim, dim> Tensor;

    static const Scalar mobilityUpwindAlpha =
            GET_PROP_VALUE(TypeTag, PTAG(MobilityUpwindAlpha));

public:
    /*!
     * \brief Evaluate the amount all conservation quantities
     *        (e.g. phase mass) within a finite sub-control volume.
     *
     *  \param result The phase mass within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    void computeStorage(PrimaryVariables &result, int scvIdx, bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const ElementVolumeVariables &elemDat = usePrevSol ? this->prevVolVars_() : this->curVolVars_();
        const VolumeVariables &vertDat = elemDat[scvIdx];

        // wetting phase mass
        result[contiWEqIdx] = vertDat.density(wPhaseIdx) * vertDat.porosity()
                * vertDat.saturation(wPhaseIdx);

        // non-wetting phase mass
        result[contiNEqIdx] = vertDat.density(nPhaseIdx) * vertDat.porosity()
                * vertDat.saturation(nPhaseIdx);
    }

    /*!
     * \brief Evaluates the mass flux over a face of a sub-control
     *        volume.
     *
     * \param flux The flux over the SCV (sub-control-volume) face for each phase
     * \param faceIdx The index of the SCV face
     */
    void computeFlux(PrimaryVariables &flux, int faceIdx) const
    {
        FluxVariables vars(this->problem_(),
                           this->elem_(),
                           this->fvElemGeom_(),
                           faceIdx,
                           this->curVolVars_());
        flux = 0;
        asImp_()->computeAdvectiveFlux(flux, vars);
        asImp_()->computeDiffusiveFlux(flux, vars);
    }

    /*!
     * \brief Evaluates the advective mass flux of all components over
     *        a face of a sub-control volume.
     *
     * \param flux The advective flux over the sub-control-volume face for each phase
     * \param fluxVars The flux variables at the current SCV
     *
     * This method is called by compute flux and is mainly there for
     * derived models to ease adding equations selectively.
     */
    void computeAdvectiveFlux(PrimaryVariables &flux, const FluxVariables &fluxVars) const
    {
        ////////
        // advective fluxes of all components in all phases
        ////////
        Vector tmpVec;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // calculate the flux in the normal direction of the
            // current sub control volume face:
            //
            // v = - (K grad p) * n
            //
            // (the minus comes from the Darcy law which states that
            // the flux is from high to low pressure potentials.)
            fluxVars.intrinsicPermeability().mv(fluxVars.potentialGrad(phaseIdx), tmpVec);
            Scalar normalFlux = - (tmpVec*fluxVars.face().normal);

            // data attached to upstream and the downstream vertices
            // of the current phase
            const VolumeVariables &up = this->curVolVars_(fluxVars.upstreamIdx(normalFlux));
            const VolumeVariables &dn = this->curVolVars_(fluxVars.downstreamIdx(normalFlux));

            // add advective flux of current component in current
            // phase
            int eqIdx = (phaseIdx == wPhaseIdx) ? contiWEqIdx : contiNEqIdx;
            flux[eqIdx] +=
                normalFlux
                *
                ((    mobilityUpwindAlpha)*up.density(phaseIdx)*up.mobility(phaseIdx)
                 +
                 (1 - mobilityUpwindAlpha)*dn.density(phaseIdx)*dn.mobility(phaseIdx));
        }
    }

    /*!
     * \brief Adds the diffusive flux to the flux vector over
     *        the face of a sub-control volume.
     *
     * \param flux The diffusive flux over the sub-control-volume face for each phase
     * \param fluxData The flux variables at the current SCV
     *
     * It doesn't do anything in two-phase model but is used by the
     * non-isothermal two-phase models to calculate diffusive heat
     * fluxes
     */
    void computeDiffusiveFlux(PrimaryVariables &flux, const FluxVariables &fluxData) const
    {
        // diffusive fluxes
        flux += 0.0;
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param q The source/sink in the SCV for each phase
     * \param localVertexIdx The index of the SCV
     *
     */
    void computeSource(PrimaryVariables &q, int localVertexIdx)
    {
        // retrieve the source term intrinsic to the problem
        this->problem_().source(q,
                                this->elem_(),
                                this->fvElemGeom_(),
                                localVertexIdx);
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
};

}

#endif
