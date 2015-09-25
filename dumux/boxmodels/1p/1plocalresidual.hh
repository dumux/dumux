// $Id: 1plocalresidual.hh 3738 2010-06-15 14:01:09Z lauser $
/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
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

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;


    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementSolutionVector)) ElementSolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePIndices)) Indices;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        pressureIdx = Indices::pressureIdx,
    };

    static const Scalar upwindWeight = GET_PROP_VALUE(TypeTag, PTAG(UpwindWeight));

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVariables)) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;

    typedef Dune::FieldVector<Scalar, dimWorld> Vector;

public:


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
        Scalar normalFlux = - (tmpVec*fluxVars.face().normal);

        const VolumeVariables &up = this->curVolVars_(fluxVars.upstreamIdx(normalFlux));
        const VolumeVariables &dn = this->curVolVars_(fluxVars.downstreamIdx(normalFlux));
        flux[pressureIdx] =
            ((    upwindWeight)*(up.density()/up.viscosity())
             +
             (1 - upwindWeight)*(dn.density()/dn.viscosity()))
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
        this->problem_().source(q,
                                this->elem_(),
                                this->fvElemGeom_(),
                                localVertexIdx);
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
};

};

#endif
