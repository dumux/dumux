/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2009 by Onur Dogan                                        *
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
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the one-phase box model.
 */
#ifndef DUMUX_1P_LOCAL_RESIDUAL_HH
#define DUMUX_1P_LOCAL_RESIDUAL_HH

#include <dumux/boxmodels/common/boxlocalresidual.hh>

#include "1pvolumevariables.hh"
#include "1pfluxvariables.hh"
#include "1pproperties.hh"

namespace Dumux
{
/*!
 * \ingroup OnePBoxModel
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the one-phase box model.
 */
template<class TypeTag>
class OnePLocalResidual : public BoxLocalResidual<TypeTag>
{
    typedef OnePLocalResidual<TypeTag> ThisType;
    typedef BoxLocalResidual<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePIndices)) Indices;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),

        pressureIdx = Indices::pressureIdx,
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVariables)) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;

    typedef Dune::FieldVector<Scalar, dimWorld> Vector;

public:

    /*!
     * \brief Constructor. Sets the upwind weight.
     */
    OnePLocalResidual()
    {
        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        upwindWeight_ = GET_PARAM(TypeTag, Scalar, UpwindWeight);
    };

    /*!
     * \brief Evaluate the rate of change of all conservation
     *        quantites (e.g. phase mass) within a sub control
     *        volume of a finite volume element for the OneP
     *        model.
     *
     * This function should not include the source and sink terms.
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
        const ElementVolumeVariables &elemVars = usePrevSol ? this->prevVolVars_() : this->curVolVars_();
        const VolumeVariables &volVars = elemVars[scvIdx];

        // partial time derivative of the wetting phase mass
        result[pressureIdx] =  volVars.density() * volVars.porosity();
    }


    /*!
     * \brief Evaluates the mass flux over a face of a subcontrol
     *        volume.
     *
     * \param flux The flux over the SCV (sub-control-volume) face
     * \param faceIdx The index of the SCV face
     */
    void computeFlux(PrimaryVariables &flux, int faceIdx) const
    {
        FluxVariables fluxVars(this->problem_(),
                               this->elem_(),
                               this->fvElemGeom_(),
                               faceIdx,
                               this->curVolVars_());
        Vector tmpVec;
        fluxVars.intrinsicPermeability().mv(fluxVars.potentialGrad(),
                                            tmpVec);

        Scalar normalFlux = 0;
        for (int i = 0; i < Vector::size; ++i)
            normalFlux += tmpVec[i]*fluxVars.face().normal[i];
        normalFlux *= -1;

        const VolumeVariables &up = this->curVolVars_(fluxVars.upstreamIdx(normalFlux));
        const VolumeVariables &dn = this->curVolVars_(fluxVars.downstreamIdx(normalFlux));
        flux[pressureIdx] =
            ((    upwindWeight_)*(up.density()/up.viscosity())
             +
             (1 - upwindWeight_)*(dn.density()/dn.viscosity()))
            *
            normalFlux;
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param q The source/sink in the SCV
     * \param localVertexIdx The index of the SCV
     *
     */
    void computeSource(PrimaryVariables &q, int localVertexIdx)
    {
        this->problem_().boxSDSource(q,
                                     this->elem_(),
                                     this->fvElemGeom_(),
                                     localVertexIdx,
                                     this->curVolVars_());
    }

    /*!
     * \brief Return the temperature given the solution vector of a
     *        finite volume.
     */
    template <class PrimaryVariables>
    Scalar temperature(const PrimaryVariables &sol)
    { return this->problem_.temperature(); /* constant temperature */ }

private:
    ThisType &asImp_()
    { return *static_cast<ThisType *>(this); }

    const ThisType &asImp_() const
    { return *static_cast<const ThisType *>(this); }

    Scalar upwindWeight_;
};

};

#endif
