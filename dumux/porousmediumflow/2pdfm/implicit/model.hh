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
* \brief Adaption of the fully implicit scheme to the two-phase flow in discrete fracture-matrix
*        model.
*/

#ifndef DUMUX_MODELS_2PDFM_MODEL_HH
#define DUMUX_MODELS_2PDFM_MODEL_HH

#include <dumux/porousmediumflow/2p/implicit/model.hh>
#include "properties.hh"

namespace Dumux
{

/*!
 * \ingroup TwoPDFMModel
 * \brief A two-phase, isothermal flow model using the fully implicit scheme.
 *
 * This model implements two-phase flow of two immiscible fluids
 * \f$\alpha \in \{ w, n \}\f$ using a standard multiphase Darcy
 * approach as the equation for the conservation of momentum, i.e.
 \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \textbf{K}
 \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} {\textbf g} \right)
 \f]
 *
 * By inserting this into the equation for the conservation of the
 * phase mass, one gets
 \f[
 \phi \frac{\partial \varrho_\alpha S_\alpha}{\partial t}
 -
 \text{div} \left\{
 \varrho_\alpha \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K} \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 \right\} - q_\alpha = 0 \;,
 \f]
 *
 * These equations are discretized by a fully-coupled vertex centered finite volume
 * (box) scheme as spatial and the implicit Euler method as time
 * discretization.
 *
 * By using constitutive relations for the capillary pressure \f$p_c =
 * p_n - p_w\f$ and relative permeability \f$k_{r\alpha}\f$ and taking
 * advantage of the fact that \f$S_w + S_n = 1\f$, the number of
 * unknowns can be reduced to two. Currently the model supports
 * choosing either \f$p_w\f$ and \f$S_n\f$ or \f$p_n\f$ and \f$S_w\f$
 * as primary variables. The formulation which ought to be used can be
 * specified by setting the <tt>Formulation</tt> property to either
 * <tt>TwoPCommonIndices::pWsN</tt> or <tt>TwoPCommonIndices::pNsW</tt>. By
 * default, the model uses \f$p_w\f$ and \f$S_n\f$.
 */
template<class TypeTag >
class TwoPDFMModel : public TwoPModel<TypeTag>
{
typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        nPhaseIdx = Indices::nPhaseIdx,
        wPhaseIdx = Indices::wPhaseIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:
    /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     *
     *        \param sol The global solution vector
     *        \param writer The writer for multi-file VTK datasets
     */
    template<class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;
        typedef Dune::BlockVector<Dune::FieldVector<double, dimWorld> > VectorField;

        // get the number of degrees of freedom
        unsigned numDofs = this->numDofs();
        unsigned numElements = this->gridView_().size(0);

        // create the required scalar fields
        ScalarField *pw = writer.allocateManagedBuffer(numDofs);
        ScalarField *pn = writer.allocateManagedBuffer(numDofs);
        ScalarField *pc = writer.allocateManagedBuffer(numDofs);
        ScalarField *swFracture = writer.allocateManagedBuffer(numDofs);
        ScalarField *snFracture = writer.allocateManagedBuffer(numDofs);
        ScalarField *swMatrix = writer.allocateManagedBuffer(numElements);
        ScalarField *snMatrix = writer.allocateManagedBuffer(numElements);
        ScalarField *rhoW = writer.allocateManagedBuffer(numDofs);
        ScalarField *rhoN = writer.allocateManagedBuffer(numDofs);
        ScalarField *mobW = writer.allocateManagedBuffer(numDofs);
        ScalarField *mobN = writer.allocateManagedBuffer(numDofs);
        ScalarField *poro = writer.allocateManagedBuffer(numDofs);
        ScalarField *Te = writer.allocateManagedBuffer(numDofs);
        VectorField *velocityN = writer.template allocateManagedBuffer<double, dimWorld>(numDofs);
        VectorField *velocityW = writer.template allocateManagedBuffer<double, dimWorld>(numDofs);
        ImplicitVelocityOutput<TypeTag> velocityOutput(this->problem_());

        if (velocityOutput.enableOutput()) // check if velocity output is demanded
        {
            // initialize velocity fields
            for (unsigned int i = 0; i < numDofs; ++i)
            {
                (*velocityN)[i] = Scalar(0);
                (*velocityW)[i] = Scalar(0);
            }
        }

        ScalarField *rank = writer.allocateManagedBuffer(numElements);

        for (const auto& element : elements(this->gridView_()))
        {
            if(element.partitionType() == Dune::InteriorEntity)
            {
                int eIdx = this->elementMapper().index(element);

                (*rank)[eIdx] = this->gridView_().comm().rank();

                FVElementGeometry fvGeometry;
                fvGeometry.update(this->gridView_(), element);

                ElementVolumeVariables elemVolVars;
                elemVolVars.update(this->problem_(),
                                   element,
                                   fvGeometry,
                                   false /* oldSol? */);
                (*swMatrix)[eIdx] = 0;
                (*snMatrix)[eIdx] = 0;

                for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
                {
                    int dofIdxGlobal = this->dofMapper().subIndex(element, scvIdx, dofCodim);

                    (*pw)[dofIdxGlobal] = elemVolVars[scvIdx].pressure(wPhaseIdx);
                    (*pn)[dofIdxGlobal] = elemVolVars[scvIdx].pressure(nPhaseIdx);
                    (*pc)[dofIdxGlobal] = elemVolVars[scvIdx].capillaryPressure();
                    (*swMatrix)[eIdx] += elemVolVars[scvIdx].saturationMatrix(wPhaseIdx)*fvGeometry.subContVol[scvIdx].volume;
                    (*snMatrix)[eIdx] += elemVolVars[scvIdx].saturationMatrix(nPhaseIdx)*fvGeometry.subContVol[scvIdx].volume;
                    (*swFracture)[dofIdxGlobal] = elemVolVars[scvIdx].saturationFracture(wPhaseIdx);
                    (*snFracture)[dofIdxGlobal] = elemVolVars[scvIdx].saturationFracture(nPhaseIdx);
                    (*rhoW)[dofIdxGlobal] = elemVolVars[scvIdx].density(wPhaseIdx);
                    (*rhoN)[dofIdxGlobal] = elemVolVars[scvIdx].density(nPhaseIdx);
                    (*mobW)[dofIdxGlobal] = elemVolVars[scvIdx].mobility(wPhaseIdx);
                    (*mobN)[dofIdxGlobal] = elemVolVars[scvIdx].mobility(nPhaseIdx);
                    (*poro)[dofIdxGlobal] = elemVolVars[scvIdx].porosity();
                    (*Te)[dofIdxGlobal] = elemVolVars[scvIdx].temperature();
                }
                (*swMatrix)[eIdx] /= fvGeometry.elementVolume;
                (*snMatrix)[eIdx] /= fvGeometry.elementVolume;
                // velocity output
                velocityOutput.calculateVelocity(*velocityW, elemVolVars, fvGeometry, element, wPhaseIdx);
                velocityOutput.calculateVelocity(*velocityN, elemVolVars, fvGeometry, element, nPhaseIdx);
            }
        }

        writer.attachDofData(*snFracture, "SnFracture", isBox);
        writer.attachDofData(*swFracture, "SwFracture", isBox);
        writer.attachCellData(*snMatrix, "SnMatrix", isBox);
        writer.attachCellData(*swMatrix, "SwMatrix", isBox);
        writer.attachDofData(*pn, "pn", isBox);
        writer.attachDofData(*pw, "pw", isBox);
        writer.attachDofData(*pc, "pc", isBox);
        writer.attachDofData(*rhoW, "rhoW", isBox);
        writer.attachDofData(*rhoN, "rhoN", isBox);
        writer.attachDofData(*mobW, "mobW", isBox);
        writer.attachDofData(*mobN, "mobN", isBox);
        writer.attachDofData(*poro, "porosity", isBox);
        writer.attachDofData(*Te, "temperature", isBox);

        if (velocityOutput.enableOutput()) // check if velocity output is demanded
        {
            writer.attachDofData(*velocityW,  "velocityW", isBox, dim);
            writer.attachDofData(*velocityN,  "velocityN", isBox, dim);
        }

        writer.attachCellData(*rank, "process rank");
    }

};
} // end namespace

#include "propertydefaults.hh"

#endif // DUMUX_MODELS_2PDFM_MODEL_HH
