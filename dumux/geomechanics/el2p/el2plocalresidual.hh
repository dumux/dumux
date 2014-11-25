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
 *
 * \brief Element-wise calculation of the residual for the linear elastic,
 * two-phase model in the fully implicit scheme.
 */
#ifndef DUMUX_ELASTIC2P_LOCAL_RESIDUAL_HH
#define DUMUX_ELASTIC2P_LOCAL_RESIDUAL_HH

#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dumux/implicit/box/boxlocalresidual.hh>
#include "el2pproperties.hh"

namespace Dumux {
/*!
 * \ingroup ElTwoPModel
 * \ingroup ImplicitLocalResidual
 *
 * \brief Element-wise calculation of the Jacobian matrix for problems
 * using the two-phase linear-elasticity fully implicit model.
 */
template<class TypeTag>
class ElTwoPLocalResidual: public BoxLocalResidual<TypeTag> {
protected:
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    enum {
        dim = GridView::dimension
    };

    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;
    typedef Dune::FieldVector<Scalar, dim> DimVector;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum {
        numFluidPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases))
    };
    enum {
        contiWEqIdx = Indices::contiWEqIdx,
        contiNEqIdx = Indices::contiNEqIdx,
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };
    //TODO: delete this if not required
    // only effective porosity update in element variables doesn't work
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridFunctionSpace)) GridFunctionSpace;
    typedef Dune::PDELab::LocalFunctionSpace<GridFunctionSpace> LocalFunctionSpace;

public:
    /*!
     * \brief Constructor. Sets the upwind weight.
     */
    ElTwoPLocalResidual() {
        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        massUpwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit,
                MassUpwindWeight);
    }
    ;

    /*!
     * \brief Evaluate the amount all conservation quantities
     *        (e.g. phase mass) within a finite sub-control volume.
     *
     *  \param storage The phase mass within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    void computeStorage(PrimaryVariables &storage, int scvIdx,
            bool usePrevSol) const {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit Euler method.
        const ElementVolumeVariables &elemVolVars =
                usePrevSol ? this->prevVolVars_() : this->curVolVars_();
        const VolumeVariables &volVars = elemVolVars[scvIdx];

        storage = Scalar(0);

        // wetting phase mass
        storage[contiWEqIdx] = volVars.density(wPhaseIdx)
                * volVars.saturation(wPhaseIdx) * volVars.effPorosity;
        // non-wetting phase mass
        storage[contiNEqIdx] = volVars.density(nPhaseIdx)
                * volVars.saturation(nPhaseIdx) * volVars.effPorosity;
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
    void computeFlux(PrimaryVariables &flux, int fIdx,
            const bool onBoundary = false) const {
        //TODO: delete this if not required
        //      adapts the effective porosity node-wise for evaluation of derivatives.
        //      At the moment computeFlux is called before computeStorage so effPorosity in
        //      computeStorage should also be correct.

//        if (fIdx == 0)
//        {
//           int numScv = this->element_().template count<dim> ();
//
//            LocalFunctionSpace localFunctionSpace(this->problem_().model().jacobianAssembler().gridFunctionSpace());
//            localFunctionSpace.bind(this->element_());
//            std::vector<Scalar> values;
//            localFunctionSpace.vread(this->problem_().model().curSol(), values);
//
//
//            typedef typename LocalFunctionSpace::template Child<1>::Type DisplacementLFS;
//            const DisplacementLFS& displacementLFS = localFunctionSpace.template child<1>();
//            const unsigned int dispSize = displacementLFS.child(0).size();
//
//                    typedef typename DisplacementLFS::template Child<0>::Type ScalarDispLFS;
//            typedef typename ScalarDispLFS::Traits::FiniteElementType::
//                            Traits::LocalBasisType::Traits::JacobianType JacobianType_V;
//            typedef typename ScalarDispLFS::Traits::FiniteElementType::
//                            Traits::LocalBasisType::Traits::RangeFieldType RF;
//            typedef typename ScalarDispLFS::Traits::FiniteElementType::
//                            Traits::LocalBasisType::Traits::DomainFieldType DF;
//
//           // calculate the divergence of the displacement
//            for (int scvIdx = 0; scvIdx < numScv; scvIdx++)
//            {
//                    const DimVector& scvCenter = this->fvGeometry_().subContVol[scvIdx].localCenter;
//
//                    // evaluate gradient of displacement shape functions
//                    std::vector<JacobianType_V> vRefShapeGradient(dispSize);
//                    displacementLFS.child(0).finiteElement().localBasis().evaluateJacobian(scvCenter, vRefShapeGradient);
//
//                    // transform gradient to physical element
//                    const Dune::FieldMatrix<DF,dim,dim> jacInvT = this->element_().geometry().jacobianInverseTransposed(scvCenter);
//                    std::vector<Dune::FieldVector<RF,dim> > vShapeGradient(dispSize);
//                    for (size_t i = 0; i < dispSize; i++)
//                    {
//                            vShapeGradient[i] = 0.0;
//                            jacInvT.umv(vRefShapeGradient[i][0],vShapeGradient[i]);
//                    }
//
//                    // calculate gradient of current displacement
//                    typedef Dune::FieldMatrix<RF, dim, dim> Tensor;
//                    Tensor uGradient(0.0);
//                    for(int comp = 0; comp < dim; ++comp){
//                            const ScalarDispLFS & scalarDispLFS = displacementLFS.child(comp);
//
//                            for (size_t i = 0; i < 8; i++)
//                                    uGradient[comp].axpy(this->curVolVars_()[i].displacement(comp), vShapeGradient[i]);
//
//                            for (size_t i = 8; i < scalarDispLFS.size(); i++)
//                                    uGradient[comp].axpy(values[scalarDispLFS.localIndex(i)], vShapeGradient[i]);
//                    }
//
//                    this->curVolVars_()[scvIdx].divU = 0.0;
//                    for (int comp = 0; comp < dim; comp++)
//                            this->curVolVars_()[scvIdx].divU += uGradient[comp][comp];
//                    if(this->problem_().coupled() == true){
//                        if (this->curVolVars_()[scvIdx].divU < - this->curVolVars_()[scvIdx].porosity()){
//                            this->curVolVars_()[scvIdx].effPorosity = this->curVolVars_()[scvIdx].porosity();
//                            std::cout<<"volume change too large"<<std::endl;
//                            }
//                        else{
//                            this->curVolVars_()[scvIdx].effPorosity = (this->curVolVars_()[scvIdx].porosity()
//                                            + this->curVolVars_()[scvIdx].divU)/(1.0 + this->curVolVars_()[scvIdx].divU);}
//                    }
//                    else
//                        this->curVolVars_()[scvIdx].effPorosity = this->curVolVars_()[scvIdx].porosity();
//            }
//        }

        FluxVariables fluxVars(this->problem_(), this->element_(),
                this->fvGeometry_(), fIdx, this->curVolVars_());

        flux = 0;
        this->computeAdvectiveFlux(flux, fluxVars);
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
    void computeAdvectiveFlux(PrimaryVariables &flux,
            const FluxVariables &fluxVars) const {
        // calculate effective permeability based on effective porosity
        // according to the relation given in Rutqvist and Tsang (2002)
        // this evaluation should be moved to another location
        DimVector tmpVec;
        DimMatrix Keff, Keff_i, Keff_j;

        Scalar exp_i, exp_j;
        // evaluate effective permeabilities for nodes i and j based on the
        // effective porosities, the initial porosities and the initial permeabilities
        exp_i =
                22.2
                        * (this->curVolVars_()[fluxVars.face().i].effPorosity
                                / this->curVolVars_()[fluxVars.face().i].porosity()
                                - 1);
        exp_j =
                22.2
                        * (this->curVolVars_()[fluxVars.face().j].effPorosity
                                / this->curVolVars_()[fluxVars.face().j].porosity()
                                - 1);
        Keff_i = fluxVars.intrinsicPermeability();
        Keff_i *= exp(exp_i);
        Keff_j = fluxVars.intrinsicPermeability();
        Keff_j *= exp(exp_j);

        // calculate the mean effective permeability at integration point
        this->problem_().spatialParams().meanK(Keff, Keff_i, Keff_j);

        // loop over all phases
        for (int phaseIdx = 0; phaseIdx < numFluidPhases; ++phaseIdx) {
            // data attached to upstream and the downstream vertices
            // of the current phase
            // calculate the flux in the normal direction of the
            // current sub control volume face

            // if geomechanical feedback on flow is taken into account the effective permeability is
            // applied for the flux calculations
            if (this->problem_().coupled() == true) {
                Keff.mv(fluxVars.potentialGrad(phaseIdx), tmpVec);
            } else {
                fluxVars.intrinsicPermeability().mv(
                        fluxVars.potentialGrad(phaseIdx), tmpVec);
            }
            Scalar normalFlux = -(tmpVec * fluxVars.face().normal);

            // data attached to upstream and the downstream vertices
            // of the current phase
            const VolumeVariables &up = this->curVolVars_(
                    fluxVars.upstreamIdx(phaseIdx));
            const VolumeVariables &dn = this->curVolVars_(
                    fluxVars.downstreamIdx(phaseIdx));

            // add advective flux of current phase
            int eqIdx = (phaseIdx == wPhaseIdx) ? contiWEqIdx : contiNEqIdx;
            flux[eqIdx] += normalFlux
                    * ((massUpwindWeight_) * up.density(phaseIdx)
                            * up.mobility(phaseIdx)
                            + (massUpwindWeight_) * dn.density(phaseIdx)
                                    * dn.mobility(phaseIdx));

            // if geomechanical feedback on flow is taken into account add the flux contribution
            // of the displacement velocity

            if (this->problem_().coupled() == true) {
                // use upwind displacement velocity to calculate phase transport (?)
                flux[eqIdx] += up.effPorosity * up.saturation(phaseIdx)
                        * up.density(phaseIdx) * fluxVars.timeDerivUNormal();
            }

        }
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param q The source/sink in the SCV for each phase
     * \param scvIdx The index of the SCV
     *
     */
    void computeSource(PrimaryVariables &q, int scvIdx) {
        // retrieve the source term intrinsic to the problem
        this->problem_().source(q, this->element_(), this->fvGeometry_(),
                scvIdx);
    }

protected:
    Implementation *asImp_() {
        return static_cast<Implementation *>(this);
    }
    const Implementation *asImp_() const {
        return static_cast<const Implementation *>(this);
    }

private:
    Scalar massUpwindWeight_;
};
}
#endif
