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
* \brief Adaption of the fully implicit box scheme to the two-phase n-component flow model.
*/

#ifndef DUMUX_2PNCMIN_MODEL_HH
#define DUMUX_2PNCMIN_MODEL_HH

#include "properties.hh"
#include "indices.hh"

#include <dumux/material/constants.hh>
#include <dumux/porousmediumflow/2pnc/implicit/model.hh>
#include "localresidual.hh"
#include <dumux/porousmediumflow/implicit/velocityoutput.hh>

namespace Dumux
{
/*!
 * \ingroup TwoPNCMinModel
 * \brief Adaption of the fully implicit scheme to the
 *        two-phase n-component fully implicit model with additional solid/mineral phases.
 *
 * This model implements two-phase n-component flow of two compressible and
 * partially miscible fluids \f$\alpha \in \{ w, n \}\f$ composed of the n components
 * \f$\kappa \in \{ w, n,\cdots \}\f$ in combination with mineral precipitation and dissolution.
 * The solid phases. The standard multiphase Darcy
 * approach is used as the equation for the conservation of momentum:
 * \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K}
 \left(\text{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 * \f]
 *
 * By inserting this into the equations for the conservation of the
 * components, one gets one transport equation for each component
 * \f{eqnarray}
 && \frac{\partial (\sum_\alpha \varrho_\alpha X_\alpha^\kappa \phi S_\alpha )}
 {\partial t}
 - \sum_\alpha  \text{div} \left\{ \varrho_\alpha X_\alpha^\kappa
 \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K}
 (\text{grad}\, p_\alpha - \varrho_{\alpha}  \mbox{\bf g}) \right\}
 \nonumber \\ \nonumber \\
    &-& \sum_\alpha \text{div} \left\{{\bf D_{\alpha, pm}^\kappa} \varrho_{\alpha} \text{grad}\, X^\kappa_{\alpha} \right\}
 - \sum_\alpha q_\alpha^\kappa = 0 \qquad \kappa \in \{w, a,\cdots \} \, ,
 \alpha \in \{w, g\}
 \f}
 *
 * The solid or mineral phases are assumed to consist of a single component.
 * Their mass balance consist only of a storage and a source term:
 *  \f$\frac{\partial \varrho_\lambda \phi_\lambda )} {\partial t}
 *  = q_\lambda\f$
 *
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme (this is not done for 2pnc approach yet, however possible) as
 * spatial and the implicit Euler method as time discretization.
 *
 * By using constitutive relations for the capillary pressure \f$p_c =
 * p_n - p_w\f$ and relative permeability \f$k_{r\alpha}\f$ and taking
 * advantage of the fact that \f$S_w + S_n = 1\f$ and \f$X^\kappa_w + X^\kappa_n = 1\f$, the number of
 * unknowns can be reduced to number of components.
 *
 * The used primary variables are, like in the two-phase model, either \f$p_w\f$ and \f$S_n\f$
 * or \f$p_n\f$ and \f$S_w\f$. The formulation which ought to be used can be
 * specified by setting the <tt>Formulation</tt> property to either
 * TwoPTwoCIndices::pWsN or TwoPTwoCIndices::pNsW. By
 * default, the model uses \f$p_w\f$ and \f$S_n\f$.
 *
 * Moreover, the second primary variable depends on the phase state, since a
 * primary variable switch is included. The phase state is stored for all nodes
 * of the system. The model is uses mole fractions.
 *Following cases can be distinguished:
 * <ul>
 *  <li> Both phases are present: The saturation is used (either \f$S_n\f$ or \f$S_w\f$, dependent on the chosen <tt>Formulation</tt>),
 *      as long as \f$ 0 < S_\alpha < 1\f$</li>.
 *  <li> Only wetting phase is present: The mole fraction of, e.g., air in the wetting phase \f$x^a_w\f$ is used,
 *      as long as the maximum mole fraction is not exceeded (\f$x^a_w<x^a_{w,max}\f$)</li>
 *  <li> Only non-wetting phase is present: The mole fraction of, e.g., water in the non-wetting phase, \f$x^w_n\f$, is used,
 *      as long as the maximum mole fraction is not exceeded (\f$x^w_n<x^w_{n,max}\f$)</li>
 * </ul>
 *
 * For the other components, the mole fraction \f$x^\kappa_w\f$ is the primary variable.
 * The primary variable of the solid phases is the volume fraction \f$\phi_\lambda = \frac{V_\lambda}{V_{total}}\f$.
 *
 * The source an sink terms link the mass balances of the n-transported component to the solid phases.
 * The porosity \f$\phi\f$ is updated according to the reduction of the initial (or solid-phase-free porous medium) porosity \f$\phi_0\f$
 * by the accumulated volume fractions of the solid phases:
 * \f$ \phi = \phi_0 - \sum (\phi_\lambda)\f$
 * Additionally, the permeability is updated depending on the current porosity.
 */

template<class TypeTag>
class TwoPNCMinModel: public TwoPNCModel<TypeTag>
{
    typedef TwoPNCMinModel<TypeTag> ThisType;
    typedef TwoPNCModel<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef Constants<Scalar> Constant;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numSPhases = GET_PROP_VALUE(TypeTag, NumSPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),
        numSecComponents = GET_PROP_VALUE(TypeTag, NumSecComponents),
        numMajorComponents = GET_PROP_VALUE(TypeTag, NumMajorComponents),

        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,

        wCompIdx = FluidSystem::wCompIdx,
        nCompIdx = FluidSystem::nCompIdx,

        wPhaseOnly = Indices::wPhaseOnly,
        nPhaseOnly = Indices::nPhaseOnly,
        bothPhases = Indices::bothPhases,

        pwsn = TwoPNCFormulation::pwsn,
        pnsw = TwoPNCFormulation::pnsw,
        formulation = GET_PROP_VALUE(TypeTag, Formulation)
    };

    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<0>::Entity Element;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldMatrix<CoordScalar, dimWorld, dimWorld> Tensor;
    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:

    /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     *
     * \param sol The solution vector
     * \param writer The writer for multi-file VTK datasets
     */
    template<class MultiWriter>
    //additional output of the permeability and the precipitate volume fractions
    void addOutputVtkFields(const SolutionVector &sol, MultiWriter &writer)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;
        typedef Dune::BlockVector<Dune::FieldVector<double, dim> > VectorField;

        // get the number of degrees of freedom
        auto numDofs = this->numDofs();

        // create the required scalar fields
        ScalarField *sn           = writer.allocateManagedBuffer(numDofs);
        ScalarField *sw           = writer.allocateManagedBuffer(numDofs);
        ScalarField *pn           = writer.allocateManagedBuffer (numDofs);
        ScalarField *pw           = writer.allocateManagedBuffer (numDofs);
        ScalarField *pc           = writer.allocateManagedBuffer (numDofs);
        ScalarField *rhoW         = writer.allocateManagedBuffer (numDofs);
        ScalarField *rhoN         = writer.allocateManagedBuffer (numDofs);
        ScalarField *mobW         = writer.allocateManagedBuffer (numDofs);
        ScalarField *mobN         = writer.allocateManagedBuffer (numDofs);
        ScalarField *phasePresence = writer.allocateManagedBuffer (numDofs);
        ScalarField *temperature  = writer.allocateManagedBuffer (numDofs);
        ScalarField *poro         = writer.allocateManagedBuffer (numDofs);
        ScalarField *permeabilityFactor = writer.allocateManagedBuffer (numDofs);
        ScalarField *precipitateVolumeFraction[numSPhases];

        for (int i = 0; i < numSPhases; ++i)
            precipitateVolumeFraction[i] = writer.allocateManagedBuffer(numDofs);

        ScalarField *massFraction[numPhases][numComponents];
        for (int i = 0; i < numPhases; ++i)
            for (int j = 0; j < numComponents; ++j)
                massFraction[i][j] = writer.allocateManagedBuffer(numDofs);

        ScalarField *moleFraction[numPhases][numComponents];
        for (int i = 0; i < numPhases; ++i)
            for (int j = 0; j < numComponents; ++j)
                moleFraction[i][j] = writer.allocateManagedBuffer(numDofs);

        ScalarField *molarity[numComponents];
        for (int j = 0; j < numComponents ; ++j)
            molarity[j] = writer.allocateManagedBuffer(numDofs);

        ScalarField *Perm[dim];
        for (int j = 0; j < dim; ++j) //Permeability only in main directions xx and yy
            Perm[j] = writer.allocateManagedBuffer(numDofs);

        VectorField *velocityN = writer.template allocateManagedBuffer<double, dim>(numDofs);
        VectorField *velocityW = writer.template allocateManagedBuffer<double, dim>(numDofs);
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

        auto numElements = this->gridView_().size(0);
        ScalarField *rank = writer.allocateManagedBuffer(numElements);

        for (const auto& element : elements(this->gridView_()))
        {
            auto eIdxGlobal = this->problem_().elementMapper().index(element);
            (*rank)[eIdxGlobal] = this->gridView_().comm().rank();
            FVElementGeometry fvGeometry;
            fvGeometry.update(this->gridView_(), element);

            ElementVolumeVariables elemVolVars;
            elemVolVars.update(this->problem_(),
                               element,
                               fvGeometry,
                               false /* oldSol? */);

            for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
            {
                auto dofIdxGlobal = this->dofMapper().subIndex(element, scvIdx, dofCodim);

                (*sn)[dofIdxGlobal] = elemVolVars[scvIdx].saturation(nPhaseIdx);
                (*sw)[dofIdxGlobal] = elemVolVars[scvIdx].saturation(wPhaseIdx);
                (*pn)[dofIdxGlobal] = elemVolVars[scvIdx].pressure(nPhaseIdx);
                (*pw)[dofIdxGlobal] = elemVolVars[scvIdx].pressure(wPhaseIdx);
                (*pc)[dofIdxGlobal] = elemVolVars[scvIdx].capillaryPressure();
                (*rhoW)[dofIdxGlobal] = elemVolVars[scvIdx].density(wPhaseIdx);
                (*rhoN)[dofIdxGlobal] = elemVolVars[scvIdx].density(nPhaseIdx);
                (*mobW)[dofIdxGlobal] = elemVolVars[scvIdx].mobility(wPhaseIdx);
                (*mobN)[dofIdxGlobal] = elemVolVars[scvIdx].mobility(nPhaseIdx);
                (*poro)[dofIdxGlobal] = elemVolVars[scvIdx].porosity();

                for (int sPhaseIdx = 0; sPhaseIdx < numSPhases; ++sPhaseIdx)
                    (*precipitateVolumeFraction[sPhaseIdx])[dofIdxGlobal] = elemVolVars[scvIdx].precipitateVolumeFraction(sPhaseIdx + numPhases);

                (*temperature)[dofIdxGlobal] = elemVolVars[scvIdx].temperature();
                (*permeabilityFactor)[dofIdxGlobal] = elemVolVars[scvIdx].permeabilityFactor();
                (*phasePresence)[dofIdxGlobal] = this->staticDat_[dofIdxGlobal].phasePresence;

                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                    for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                        (*massFraction[phaseIdx][compIdx])[dofIdxGlobal]= elemVolVars[scvIdx].massFraction(phaseIdx,compIdx);

                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                    for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                        (*moleFraction[phaseIdx][compIdx])[dofIdxGlobal]= elemVolVars[scvIdx].moleFraction(phaseIdx,compIdx);

                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                    (*molarity[compIdx])[dofIdxGlobal] = (elemVolVars[scvIdx].molarity(wPhaseIdx, compIdx));

                Tensor K = this->perm_(this->problem_().spatialParams().intrinsicPermeability(element, fvGeometry, scvIdx));

                for (int j = 0; j<dim; ++j)
                    (*Perm[j])[dofIdxGlobal] = K[j][j] * elemVolVars[scvIdx].permeabilityFactor();
            };

            // velocity output
            if(velocityOutput.enableOutput()){
                velocityOutput.calculateVelocity(*velocityW, elemVolVars, fvGeometry, element, wPhaseIdx);
                velocityOutput.calculateVelocity(*velocityN, elemVolVars, fvGeometry, element, nPhaseIdx);
            }
        }  // loop over element

        writer.attachDofData(*sn, "Sn", isBox);
        writer.attachDofData(*sw, "Sw", isBox);
        writer.attachDofData(*pn, "pn", isBox);
        writer.attachDofData(*pw, "pw", isBox);
        writer.attachDofData(*pc, "pc", isBox);
        writer.attachDofData(*rhoW, "rhoW", isBox);
        writer.attachDofData(*rhoN, "rhoN", isBox);
        writer.attachDofData(*mobW, "mobW", isBox);
        writer.attachDofData(*mobN, "mobN", isBox);
        writer.attachDofData(*poro, "porosity", isBox);
        writer.attachDofData(*permeabilityFactor, "permeabilityFactor", isBox);
        writer.attachDofData(*temperature, "temperature", isBox);
        writer.attachDofData(*phasePresence, "phase presence", isBox);

        for (int i = 0; i < numSPhases; ++i)
        {
            std::ostringstream oss;
            oss << "precipitateVolumeFraction_" << FluidSystem::phaseName(numPhases + i);
            writer.attachDofData(*precipitateVolumeFraction[i], oss.str(), isBox);
        }

        writer.attachDofData(*Perm[0], "Kxx", isBox);
        if (dim >= 2)
            writer.attachDofData(*Perm[1], "Kyy", isBox);
        if (dim == 3)
            writer.attachDofData(*Perm[2], "Kzz", isBox);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                std::ostringstream oss;
                oss << "X_" << FluidSystem::phaseName(phaseIdx) << "^" << FluidSystem::componentName(compIdx);
                writer.attachDofData(*massFraction[phaseIdx][compIdx], oss.str(), isBox);
            }
        }

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                std::ostringstream oss;
                oss << "x_" << FluidSystem::phaseName(phaseIdx) << "^" << FluidSystem::componentName(compIdx);
                writer.attachDofData(*moleFraction[phaseIdx][compIdx], oss.str(), isBox);
            }
        }

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            std::ostringstream oss;
            oss << "m_" << FluidSystem::phaseName(wPhaseIdx) << "^" << FluidSystem::componentName(compIdx);
            writer.attachDofData(*molarity[compIdx], oss.str(), isBox);
        }

        if (velocityOutput.enableOutput()) // check if velocity output is demanded
        {
            writer.attachDofData(*velocityW,  "velocityW", isBox, dim);
            writer.attachDofData(*velocityN,  "velocityN", isBox, dim);
        }

        writer.attachCellData(*rank, "process rank");
    }

    /*!
     * \brief Update the static data of all vertices in the grid.
     *
     * \param curGlobalSol The current global solution
     * \param oldGlobalSol The previous global solution
     */
    void updateStaticData(SolutionVector &curGlobalSol,
                          const SolutionVector &oldGlobalSol)
    {
        bool wasSwitched = false;

        for (unsigned i = 0; i < this->staticDat_.size(); ++i)
            this->staticDat_[i].visited = false;

        for (const auto& element : elements(this->gridView_()))
        {
            FVElementGeometry fvGeometry;
            fvGeometry.update(this->gridView_(), element);
            for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
            {
                auto dofIdxGlobal = this->dofMapper().subIndex(element, scvIdx, dofCodim);

                if (this->staticDat_[dofIdxGlobal].visited)
                    continue;

                this->staticDat_[dofIdxGlobal].visited = true;
                VolumeVariables volVars;
                volVars.update(curGlobalSol[dofIdxGlobal],
                               this->problem_(),
                               element,
                               fvGeometry,
                               scvIdx,
                               false);
                auto global = element.geometry().corner(scvIdx);
                if (primaryVarSwitch_(curGlobalSol, volVars, dofIdxGlobal, global))
                    wasSwitched = true;
            }
        }

        // make sure that if there was a variable switch in an
        // other partition we will also set the switch flag
        // for our partition.
        if (this->gridView_().comm().size() > 1)
            wasSwitched = this->gridView_().comm().max(wasSwitched);

        this->setSwitched_(wasSwitched);
    }
protected:

    /*!
     * \brief Set whether there was a primary variable switch after in
     *        the last timestep.
     */
    bool primaryVarSwitch_(SolutionVector &globalSol,
                           const VolumeVariables &volVars,
                           int dofIdxGlobal,
                           const GlobalPosition &globalPos)
    {
        // evaluate primary variable switch
        bool wouldSwitch = false;
        int phasePresence = this->staticDat_[dofIdxGlobal].phasePresence;
        int newPhasePresence = phasePresence;

        //check if a primary variable switch is necessary
        if (phasePresence == bothPhases)
        {
            Scalar Smin = 0.0; //saturation threshold
            if (this->staticDat_[dofIdxGlobal].wasSwitched)
                Smin = -0.01;

            //if saturation of wetting phase is smaller 0 switch
            if (volVars.saturation(wPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                //wetting phase has to disappear
                std::cout << "Wetting Phase disappears at vertex " << dofIdxGlobal
                            << ", coordinated: " << globalPos << ", Sw: "
                            << volVars.saturation(wPhaseIdx) << std::endl;
                newPhasePresence = nPhaseOnly;

                //switch saturation to xnH20 (not depending on formulation)
                globalSol[dofIdxGlobal][switchIdx]
                        = volVars.moleFraction(nPhaseIdx, wCompIdx /*H2O*/);
                //Here unlike 2pnc model we do not switch all components to to mole fraction in nonwetting phase
            }
            //if saturation of nonwetting phase is smaller than 0 switch
            else if (volVars.saturation(nPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                //nonwetting phase has to disappear
                std::cout << "Nonwetting Phase disappears at vertex " << dofIdxGlobal
                            << ", coordinated: " << globalPos << ", Sn: "
                            << volVars.saturation(nPhaseIdx) << std::endl;
                newPhasePresence = wPhaseOnly;

                //switch saturation to xwN2 (not depending on formulation)
                globalSol[dofIdxGlobal][switchIdx] = volVars.moleFraction(wPhaseIdx, nCompIdx /*N2*/);
            }
        }
        else if (phasePresence == nPhaseOnly)
        {
            Scalar sumxw = 0;
            //Calculate sum of mole fractions in the hypothetical wetting phase
            for (int compIdx = 0; compIdx < numComponents; compIdx++)
            {
                sumxw += volVars.moleFraction(wPhaseIdx, compIdx);
            }
            Scalar xwmax = 1.0;
            if (sumxw > xwmax)
                wouldSwitch = true;
            if (this->staticDat_[dofIdxGlobal].wasSwitched)
                xwmax *=1.02;

            //if the sum of the mole fractions would be larger than
            //1, wetting phase appears
            if (sumxw/*sum of mole fractions*/ > xwmax/*1*/)
            {
                // wetting phase appears
                std::cout << "Wetting Phase appears at vertex " << dofIdxGlobal
                          << ", coordinated: " << globalPos << ", sumxw: "
                          << sumxw << std::endl;
                newPhasePresence = bothPhases;
                if (formulation == pnsw)
                    globalSol[dofIdxGlobal][switchIdx] = 0.0;
                else if (formulation == pwsn)
                    globalSol[dofIdxGlobal][switchIdx] = 1.0;
                //Here unlike 2pnc model we do not switch all components to to mole fraction in nonwetting phase
            }
        }
        else if (phasePresence == wPhaseOnly)
        {
            Scalar xnmax = 1;
            Scalar sumxn = 0;
            //Calculate sum of mole fractions in the hypothetical nonwetting phase
            for (int compIdx = 0; compIdx < numComponents; compIdx++)
            {
                sumxn += volVars.moleFraction(nPhaseIdx, compIdx);
            }
            if (sumxn > xnmax)
                wouldSwitch = true;
            if (this->staticDat_[dofIdxGlobal].wasSwitched)
                xnmax *=1.02;
            //nonwetting phase appears if sum is larger than one
            if (sumxn > xnmax)
            {
                std::cout << "Nonwetting Phase appears at vertex " << dofIdxGlobal
                          << ", coordinated: " << globalPos << ", sumxn: "
                          << sumxn << std::endl;
                newPhasePresence = bothPhases;
                //saturation of the wetting phase set to 0.9999 (if formulation pnsw and vice versa)
                if (formulation == pnsw)
                    globalSol[dofIdxGlobal][switchIdx] = 0.999;
                else if (formulation == pwsn)
                    globalSol[dofIdxGlobal][switchIdx] = 0.001;

            }
        }
        this->staticDat_[dofIdxGlobal].phasePresence = newPhasePresence;
        this->staticDat_[dofIdxGlobal].wasSwitched = wouldSwitch;
        return phasePresence != newPhasePresence;
    }
};

}

#include "propertydefaults.hh"

#endif
