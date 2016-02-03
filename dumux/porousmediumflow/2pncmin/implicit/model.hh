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
 *        two-phase n-component fully implicit model.
 *
 * This model implements two-phase n-component flow of two compressible and
 * partially miscible fluids \f$\alpha \in \{ w, n \}\f$ composed of the n components
 * \f$\kappa \in \{ w, a,\cdots \}\f$. The standard multiphase Darcy
 * approach is used as the equation for the conservation of momentum:
 * \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K}
 \left(\text{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 * \f]
 *
 * By inserting this into the equations for the conservation of the
 * components, one gets one transport equation for each component
 * \f{eqnarray}
 && \phi \frac{\partial (\sum_\alpha \varrho_\alpha X_\alpha^\kappa S_\alpha )}
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
 *  <li> Only wetting phase is present: The mass fraction of, e.g., air in the wetting phase \f$X^a_w\f$ is used,
 *      as long as the maximum mass fraction is not exceeded (\f$X^a_w<X^a_{w,max}\f$)</li>
 *  <li> Only non-wetting phase is present: The mass fraction of, e.g., water in the non-wetting phase, \f$X^w_n\f$, is used,
 *      as long as the maximum mass fraction is not exceeded (\f$X^w_n<X^w_{n,max}\f$)</li>
 * </ul>
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
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, ElementMapper) ElementMapper;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef Dumux::Constants<Scalar> Constant;

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

        plSg = TwoPNCFormulation::plSg,
        pgSl = TwoPNCFormulation::pgSl,
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
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;
        typedef Dune::BlockVector<Dune::FieldVector<double, dim> > VectorField;

        // get the number of degrees of freedom
        unsigned numDofs = this->numDofs();

        // create the required scalar fields
        ScalarField *Sg           = writer.allocateManagedBuffer (numDofs);
        ScalarField *Sl        = writer.allocateManagedBuffer (numDofs);
        ScalarField *pg           = writer.allocateManagedBuffer (numDofs);
        ScalarField *pl           = writer.allocateManagedBuffer (numDofs);
        ScalarField *pc        = writer.allocateManagedBuffer (numDofs);
        ScalarField *rhoL       = writer.allocateManagedBuffer (numDofs);
        ScalarField *rhoG       = writer.allocateManagedBuffer (numDofs);
        ScalarField *mobL       = writer.allocateManagedBuffer (numDofs);
        ScalarField *mobG        = writer.allocateManagedBuffer (numDofs);
        ScalarField *phasePresence = writer.allocateManagedBuffer (numDofs);
        ScalarField *temperature   = writer.allocateManagedBuffer (numDofs);
        ScalarField *poro          = writer.allocateManagedBuffer (numDofs);
        ScalarField *boxVolume     = writer.allocateManagedBuffer (numDofs);
        ScalarField *cellNum        = writer.allocateManagedBuffer (numDofs);
        ScalarField *permeabilityFactor    = writer.allocateManagedBuffer (numDofs);
        ScalarField *precipitateVolumeFraction[numSPhases] ;

        for (int i = 0; i < numSPhases; ++i)
        {
            precipitateVolumeFraction[i]= writer.allocateManagedBuffer (numDofs);
        }

        ScalarField *massFraction[numPhases][numComponents];
        for (int i = 0; i < numPhases; ++i)
            for (int j = 0; j < numComponents; ++j)
                massFraction[i][j] = writer.allocateManagedBuffer(numDofs);

        ScalarField *molarity[numComponents];
        for (int j = 0; j < numComponents ; ++j)
            molarity[j] = writer.allocateManagedBuffer(numDofs);

        ScalarField *Perm[dim];
        for (int j = 0; j < dim; ++j) //Permeability only in main directions xx and yy
            Perm[j] = writer.allocateManagedBuffer(numDofs);

        *boxVolume = 0;

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
                (*cellNum)[i] = Scalar(0.0);
            }
        }

        unsigned numElements = this->gridView_().size(0);
        ScalarField *rank =
                writer.allocateManagedBuffer (numElements);

        FVElementGeometry fvGeometry;
        VolumeVariables volVars;
        ElementVolumeVariables elemVolVars;

        for (const auto& element : elements(this->gridView_()))
        {
            int idx = this->problem_().elementMapper().index(element);
            (*rank)[idx] = this->gridView_().comm().rank();
            fvGeometry.update(this->gridView_(), element);

            elemVolVars.update(this->problem_(),
                    element,
                    fvGeometry,
                    false /* oldSol? */);

            int numVerts = element.subEntities(dim);

            for (int i = 0; i < numVerts; ++i)
            {
                int globalIdx = this->vertexMapper().subIndex(element, i, dim);
                volVars.update(sol[globalIdx],
                               this->problem_(),
                               element,
                               fvGeometry,
                               i,
                               false);

                (*Sg)[globalIdx]              = volVars.saturation(nPhaseIdx);
                (*Sl)[globalIdx]              = volVars.saturation(wPhaseIdx);
                (*pg)[globalIdx]                = volVars.pressure(nPhaseIdx);
                (*pl)[globalIdx]                 = volVars.pressure(wPhaseIdx);
                (*pc)[globalIdx]              = volVars.capillaryPressure();
                (*rhoL)[globalIdx]               = volVars.density(wPhaseIdx);
                (*rhoG)[globalIdx]               = volVars.density(nPhaseIdx);
                (*mobL)[globalIdx]               = volVars.mobility(wPhaseIdx);
                (*mobG)[globalIdx]               = volVars.mobility(nPhaseIdx);
                (*boxVolume)[globalIdx]        += fvGeometry.subContVol[i].volume;
                (*poro)[globalIdx]              = volVars.porosity();

                for (int sPhaseIdx = 0; sPhaseIdx < numSPhases; ++sPhaseIdx)
                {
                    (*precipitateVolumeFraction[sPhaseIdx])[globalIdx]   = volVars.precipitateVolumeFraction(sPhaseIdx + numPhases);
                }
                (*temperature)[globalIdx]      = volVars.temperature();
                (*permeabilityFactor)[globalIdx]      = volVars.permeabilityFactor();
                (*phasePresence)[globalIdx]  = this->staticDat_[globalIdx].phasePresence;

                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                    for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                    {
                        (*massFraction[phaseIdx][compIdx])[globalIdx]= volVars.massFraction(phaseIdx,compIdx);

                        Valgrind::CheckDefined((*massFraction[phaseIdx][compIdx])[globalIdx]);

                    }
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                    (*molarity[compIdx])[globalIdx] = (volVars.molarity(wPhaseIdx, compIdx));

                Tensor K = this->perm_(this->problem_().spatialParams().intrinsicPermeability(element, fvGeometry, i));

                for (int j = 0; j<dim; ++j)
                    (*Perm[j])[globalIdx] = K[j][j] * volVars.permeabilityFactor();
            };

            // velocity output
            if(velocityOutput.enableOutput()){
                velocityOutput.calculateVelocity(*velocityW, elemVolVars, fvGeometry, element, wPhaseIdx);
                velocityOutput.calculateVelocity(*velocityN, elemVolVars, fvGeometry, element, nPhaseIdx);
            }
        }  // loop over element

        writer.attachVertexData(*Sg, "Sg");
        writer.attachVertexData(*Sl, "Sl");
        writer.attachVertexData(*pg, "pg");
        writer.attachVertexData(*pl, "pl");
        writer.attachVertexData(*pc, "pc");
        writer.attachVertexData(*rhoL, "rhoL");
        writer.attachVertexData(*rhoG, "rhoG");
        writer.attachVertexData(*mobL, "mobL");
        writer.attachVertexData(*mobG, "mobG");
        writer.attachVertexData(*poro, "porosity");
        writer.attachVertexData(*permeabilityFactor, "permeabilityFactor");
        writer.attachVertexData(*temperature, "temperature");
        writer.attachVertexData(*phasePresence, "phase presence");
        writer.attachVertexData(*boxVolume, "boxVolume");


        for (int i = 0; i < numSPhases; ++i)
        {
            std::ostringstream oss;
            oss << "precipitateVolumeFraction_"
                << FluidSystem::phaseName(numPhases + i);
            writer.attachDofData(*precipitateVolumeFraction[i], oss.str().c_str(), isBox);
        }

        writer.attachVertexData(*Perm[0], "Kxx");
        if (dim >= 2)
            writer.attachVertexData(*Perm[1], "Kyy");
        if (dim == 3)
            writer.attachVertexData(*Perm[2], "Kzz");

        for (int i = 0; i < numPhases; ++i)
        {
            for (int j = 0; j < numComponents; ++j)
            {
                std::ostringstream oss;
                oss << "X^"
                << FluidSystem::phaseName(i)
                << "_"
                << FluidSystem::componentName(j);
                writer.attachVertexData(*massFraction[i][j], oss.str().c_str());
            }
        }

        for (int j = 0; j < numComponents; ++j)
        {
            std::ostringstream oss;
            oss << "m^w_"
                << FluidSystem::componentName(j);
            writer.attachVertexData(*molarity[j], oss.str().c_str());
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

        FVElementGeometry fvGeometry;
        static VolumeVariables volVars;
        for (const auto& element : elements(this->gridView_()))
        {
            fvGeometry.update(this->gridView_(), element);
            for (int i = 0; i < fvGeometry.numScv; ++i)
            {
                int globalIdx = this->vertexMapper().subIndex(element, i, dim);

                if (this->staticDat_[globalIdx].visited)
                    continue;

                this->staticDat_[globalIdx].visited = true;
                volVars.update(curGlobalSol[globalIdx],
                               this->problem_(),
                               element,
                               fvGeometry,
                               i,
                               false);
                const GlobalPosition &global = element.geometry().corner(i);
                if (primaryVarSwitch_(curGlobalSol,
                                      volVars,
                                      globalIdx,
                                      global))
                { wasSwitched = true;
                }
            }
        }

        // make sure that if there was a variable switch in an
        // other partition we will also set the switch flag
        // for our partition.
        if (this->gridView_().comm().size() > 1)
        wasSwitched = this->gridView_().comm().max(wasSwitched);

        setSwitched_(wasSwitched);
    }
protected:

    /*!
     * \brief Set whether there was a primary variable switch after in
     *        the last timestep.
     */
    void setSwitched_(bool yesno)
    {
        switchFlag_ = yesno;
    }

    /*!
     * \brief Set whether there was a primary variable switch after in
     *        the last timestep.
     */
    bool primaryVarSwitch_(SolutionVector &globalSol,
                           const VolumeVariables &volVars, int globalIdx,
                           const GlobalPosition &globalPos)
    {
            // evaluate primary variable switch
            bool wouldSwitch = false;
            int phasePresence = this->staticDat_[globalIdx].phasePresence;
            int newPhasePresence = phasePresence;

            //check if a primary variable switch is necessary
            if (phasePresence == bothPhases)
            {
                Scalar Smin = 0.0; //saturation threshold
                if (this->staticDat_[globalIdx].wasSwitched)
                    Smin = -0.01;

                //if saturation of liquid phase is smaller 0 switch
                if (volVars.saturation(wPhaseIdx) <= Smin)
                {
                    wouldSwitch = true;
                    //liquid phase has to disappear
                    std::cout << "Liquid Phase disappears at vertex " << globalIdx
                                << ", coordinated: " << globalPos << ", Sl: "
                                << volVars.saturation(wPhaseIdx) << std::endl;
                    newPhasePresence = nPhaseOnly;

                    //switch not depending on formulation
                    //switch "Sl" to "xgH20"
                    globalSol[globalIdx][switchIdx]
                            = volVars.moleFraction(nPhaseIdx, wCompIdx /*H2O*/);
                    //Here unlike 2pnc model we do not switch all components to to mole fraction in gas phase
                }
                //if saturation of gas phase is smaller than 0 switch
                else if (volVars.saturation(nPhaseIdx) <= Smin)
                {
                    wouldSwitch = true;
                    //gas phase has to disappear
                    std::cout << "Gas Phase disappears at vertex " << globalIdx
                                << ", coordinated: " << globalPos << ", Sg: "
                                << volVars.saturation(nPhaseIdx) << std::endl;
                    newPhasePresence = wPhaseOnly;

                    //switch "Sl" to "xlN2"
                    globalSol[globalIdx][switchIdx]
                            = volVars.moleFraction(wPhaseIdx, nCompIdx /*N2*/);
                }
            }
            else if (phasePresence == nPhaseOnly)
            {
            Scalar sumxl = 0;
            //Calculate sum of mole fractions (water and air) in the hypothetical liquid phase
            for (int compIdx = 0; compIdx < numComponents; compIdx++)
                {
                    sumxl += volVars.moleFraction(wPhaseIdx, compIdx);
                }
                    Scalar xlmax = 1.0;
                    if (sumxl > xlmax)
                wouldSwitch = true;
                    if (this->staticDat_[globalIdx].wasSwitched)
                xlmax *=1.02;

            //if the sum of the mole fractions would be larger than
            //1, wetting phase appears
                    if (sumxl/*sum of mole fractions*/ > xlmax/*1*/)
                    {
                        // liquid phase appears
                        std::cout << "Liquid Phase appears at vertex " << globalIdx
                                << ", coordinated: " << globalPos << ", sumxl: "
                                << sumxl << std::endl;
                        newPhasePresence = bothPhases;
                        if (formulation == pgSl)
                            globalSol[globalIdx][switchIdx] = 0.0;
                        else if (formulation == plSg)
                            globalSol[globalIdx][switchIdx] = 1.0;
                    //Here unlike 2pnc model we do not switch all components to to mole fraction in gas phase
                    }
            }
            else if (phasePresence == wPhaseOnly)
            {
                Scalar xgmax = 1;
                Scalar sumxg = 0;
                //Calculate sum of mole fractions in the hypothetical gas phase
                for (int compIdx = 0; compIdx < numComponents; compIdx++)
                {
                    sumxg += volVars.moleFraction(nPhaseIdx, compIdx);
                }
                if (sumxg > xgmax)
                    wouldSwitch = true;
                if (this->staticDat_[globalIdx].wasSwitched)
                    xgmax *=1.02;
                //liquid phase appears if sum is larger than one
                if (sumxg > xgmax)
                {
                    std::cout << "Gas Phase appears at vertex " << globalIdx
                            << ", coordinated: " << globalPos << ", sumxg: "
                            << sumxg << std::endl;
                    newPhasePresence = bothPhases;
                    //saturation of the liquid phase set to 0.9999 (if formulation pgSl and vice versa)
                    if (formulation == pgSl)
                        globalSol[globalIdx][switchIdx] = 0.999;
                    else if (formulation == plSg)
                        globalSol[globalIdx][switchIdx] = 0.001;

                }
            }
            this->staticDat_[globalIdx].phasePresence = newPhasePresence;
            this->staticDat_[globalIdx].wasSwitched = wouldSwitch;
            return phasePresence != newPhasePresence;
        }
        // parameters given in constructor
        bool switchFlag_;
};

}

#include "propertydefaults.hh"

#endif
