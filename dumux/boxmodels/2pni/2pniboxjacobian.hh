/*****************************************************************************
 *   Copyright (C) 2008-2009 by Melanie Darcis                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUMUX_NEW_2PNI_BOX_JACOBIAN_HH
#define DUMUX_NEW_2PNI_BOX_JACOBIAN_HH

#include <dumux/boxmodels/2p/2pboxjacobian.hh>

#include <dumux/boxmodels/2pni/2pnielementdata.hh>
#include <dumux/boxmodels/2pni/2pnivertexdata.hh>
#include <dumux/boxmodels/2pni/2pnifluxdata.hh>

#include <dumux/boxmodels/2pni/2pniproperties.hh>

namespace Dune
{
/** 
 * \brief The local jacobian operator for the non-isothermal
 *        two-phase model.
 */
template<class TypeTag>
class TwoPNIBoxJacobian 
    : public TwoPBoxJacobianBase<TypeTag, TwoPNIBoxJacobian<TypeTag> >
{
    typedef TwoPNIBoxJacobian<TypeTag>               ThisType;
    typedef TwoPBoxJacobianBase<TypeTag, ThisType>   ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))   Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView))  GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))    Scalar;

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVarVector        PrimaryVarVector;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    enum {
        dim              = GridView::dimension,
        dimWorld         = GridView::dimensionworld,

        numPhases        = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),
        temperatureIdx   = Indices::temperatureIdx,

        wPhase   = Indices::wPhase,
        nPhase   = Indices::nPhase
    };


    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexData))   VertexData;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxData))     FluxData;
    typedef std::vector<VertexData> VertexDataArray;

    typedef Dune::FieldVector<Scalar, dim>       LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld>  GlobalPosition;

    static const Scalar mobilityUpwindAlpha = GET_PROP_VALUE(TypeTag, PTAG(MobilityUpwindAlpha));

public:
    TwoPNIBoxJacobian(Problem &problem)
        : ParentType(problem)
    {
    };

    /*!
     * \brief Evaluate the amount all conservation quantites
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     */
    void computeStorage(PrimaryVarVector &result, int scvIdx, bool usePrevSol) const
    {
        // compute the storage term for phase mass
        ParentType::computeStorage(result, scvIdx, usePrevSol);

        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const VertexDataArray &elemDat = usePrevSol ? this->prevElemDat_  : this->curElemDat_;
        const VertexData  &vertDat = elemDat[scvIdx];

        result[temperatureIdx] =
            vertDat.porosity*(vertDat.density[wPhase] *
                              vertDat.intEnergy[wPhase] *
                              vertDat.satW
                              +
                              vertDat.density[nPhase] *
                              vertDat.intEnergy[nPhase] *
                              vertDat.satN)
            +
            vertDat.temperature *
            this->problem_.soil().heatCap(this->curElementGeom_.elementGlobal,
                                          this->curElement_(),
                                          this->curElementGeom_.elementLocal);
    }

    /*!
     * \brief Sets the temperature term of the flux vector to the
     *        heat flux due to advection of the fluids.
     */
    void computeAdvectiveFlux(PrimaryVarVector &flux,
                              const FluxData &fluxData) const
    {
        // advective mass flux
        ParentType::computeAdvectiveFlux(flux, fluxData);

        // advective heat flux in all phases
        flux[temperatureIdx] = 0;
        for (int phase = 0; phase < numPhases; ++phase) {
            // vertex data of the upstream and the downstream vertices
            const VertexData &up = this->curElemDat_[fluxData.upstreamIdx[phase]];
            const VertexData &dn = this->curElemDat_[fluxData.downstreamIdx[phase]];

            flux[temperatureIdx] +=
                fluxData.vDarcyNormal[phase] * (
                    mobilityUpwindAlpha * // upstream vertex
                    (  up.density[phase] *
                       up.mobility[phase] *
                       up.enthalpy[phase])
                    +
                    (1 - mobilityUpwindAlpha) * // downstream vertex
                    (  dn.density[phase] *
                       dn.mobility[phase] *
                       dn.enthalpy[phase]) );
        }
    }

    /*!
     * \brief Adds the diffusive heat flux to the flux vector over
     *        the face of a sub-control volume.
     */
    void computeDiffusiveFlux(PrimaryVarVector &flux,
                              const FluxData &fluxData) const
    {
        // diffusive mass flux
        ParentType::computeDiffusiveFlux(flux, fluxData);

        // diffusive heat flux
        flux[temperatureIdx] += (fluxData.temperatureGrad*fluxData.face->normal)*fluxData.heatCondAtIp;
    }

    // internal method!
    template <class PrimaryVarVector>
    Scalar temperature(const PrimaryVarVector &sol)
    { return sol[temperatureIdx]; }
};

}

#endif
