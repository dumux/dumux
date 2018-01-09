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
 * \brief Base class for all models which use the single-phase,
 *        two-component fully implicit model.
 *        Adaption of the fully implicit scheme to the one-phase two-component flow model.
 */

#ifndef DUMUX_ONEP_TWOC_ADSORPTION_MODEL_HH
#define DUMUX_ONEP_TWOC_ADSORPTION_MODEL_HH

#include <dumux/porousmediumflow/implicit/velocityoutput.hh>
#include "properties.hh"

namespace Dumux
{

/*!
 * \ingroup OnePTwoCAdsorptionModel
 * \brief Adaption of the fully implicit scheme to the one-phase two-component flow model.
 *
 * This model implements a one-phase flow of a compressible fluid, that consists of two components,
 * using a standard Darcy
 * approach as the equation for the conservation of momentum:
 \f[
 v = - \frac{\textbf K}{\mu}
 \left(\textbf{grad}\, p - \varrho {\textbf g} \right)
 \f]
 *
 * Gravity can be enabled or disabled via the property system.
 * By inserting this into the continuity equation, one gets
 \f[
 \phi\frac{\partial \varrho}{\partial t}
 - \text{div} \left\{
   \varrho \frac{\textbf K}{\mu}  \left(\textbf{grad}\, p - \varrho {\textbf g} \right)
   + \sum_\kappa \varrho D^\kappa_\text{pm} \frac{M^\kappa}{M_\alpha} \textbf{grad} x^\kappa
 \right\} = q \;,
 \f]
 *
 * The transport of the components \f$\kappa \in \{ w, a \}\f$ is described by the following equation:
 \f[
 \phi \frac{ \partial \varrho X^\kappa}{\partial t}
 - \text{div} \left\lbrace \varrho X^\kappa \frac{{\textbf K}}{\mu} \left( \textbf{grad}\, p -
 \varrho {\textbf g} \right)
 + \varrho D^\kappa_\text{pm} \frac{M^\kappa}{M_\alpha} \textbf{grad} x^\kappa \right\rbrace = q.
 \f]
 *
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial
 * and the implicit Euler method as time discretization.
 * The model is able to use either mole or mass fractions. The property useMoles can be set to either true or false in the
 * problem file. Make sure that the according units are used in the problem setup. useMoles is set to true by default.
 *
 * The primary variables are the pressure \f$p\f$ and the mole or mass fraction of dissolved component \f$x\f$.
 */

template<class TypeTag >
class OnePTwoCAdsorptionModel : public GET_PROP_TYPE(TypeTag, BaseModel)
{
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum { phaseIdx = Indices::phaseIdx };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:
    /*!
     * \brief \copybrief ImplicitModel::addOutputVtkFields
     *
     * Specialization for the OnePTwoCModel, adding pressure,
     * mass and mole fractions, and the process rank to the VTK writer.
     */
    template<class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;
        typedef Dune::BlockVector<Dune::FieldVector<double, dimWorld> > VectorField;

        // create the required scalar fields
        unsigned numDofs = this->numDofs();
        ScalarField &pressure = *writer.allocateManagedBuffer(numDofs);
        ScalarField &delp     = *writer.allocateManagedBuffer(numDofs);
        ScalarField &moleFraction0 = *writer.allocateManagedBuffer(numDofs);
        ScalarField &moleFraction1 = *writer.allocateManagedBuffer(numDofs);
        ScalarField &massFraction0 = *writer.allocateManagedBuffer(numDofs);
        ScalarField &massFraction1 = *writer.allocateManagedBuffer(numDofs);
        ScalarField &rho      = *writer.allocateManagedBuffer(numDofs);
        ScalarField &rhoMol   = *writer.allocateManagedBuffer(numDofs);
        ScalarField &mu       = *writer.allocateManagedBuffer(numDofs);
        VectorField *velocity = writer.template allocateManagedBuffer<double, dimWorld>(numDofs);
        ImplicitVelocityOutput<TypeTag> velocityOutput(this->problem_());
        ScalarField *S_n_CH4    = writer.allocateManagedBuffer (numDofs);
        ScalarField *S_Ads_CH4  = writer.allocateManagedBuffer (numDofs);
        ScalarField *S_nAds_CH4 = writer.allocateManagedBuffer (numDofs);
        ScalarField *S_n_CO2    = writer.allocateManagedBuffer (numDofs);
        ScalarField *S_Ads_CO2  = writer.allocateManagedBuffer (numDofs);
        ScalarField *S_nAds_CO2 = writer.allocateManagedBuffer (numDofs);
        ScalarField *boxVolume        = writer.allocateManagedBuffer (numDofs);
        ScalarField *porosity        = writer.allocateManagedBuffer (numDofs);

        *boxVolume = 0;

        if (velocityOutput.enableOutput())
        {
            // initialize velocity field
            for (unsigned int i = 0; i < numDofs; ++i)
            {
                (*velocity)[i] = Scalar(0);
            }
        }

        unsigned numElements = this->gridView_().size(0);
        ScalarField &rank = *writer.allocateManagedBuffer(numElements);

        for (const auto& element : elements(this->gridView_(), Dune::Partitions::interior))
        {
           int eIdx = this->problem_().model().elementMapper().index(element);

            rank[eIdx] = this->gridView_().comm().rank();

            FVElementGeometry fvGeometry;
            fvGeometry.update(this->gridView_(), element);

            ElementVolumeVariables elemVolVars;
            elemVolVars.update(this->problem_(),
                               element,
                               fvGeometry,
                               false /* oldSol? */);

            for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
            {
                int dofIdxGlobal = this->dofMapper().subIndex(element, scvIdx, dofCodim);

                pressure[dofIdxGlobal] = elemVolVars[scvIdx].pressure();
                delp[dofIdxGlobal]     = elemVolVars[scvIdx].pressure() - 1e5;
                moleFraction0[dofIdxGlobal] = elemVolVars[scvIdx].moleFraction(0);
                moleFraction1[dofIdxGlobal] = elemVolVars[scvIdx].moleFraction(1);
                massFraction0[dofIdxGlobal] = elemVolVars[scvIdx].massFraction(0);
                massFraction1[dofIdxGlobal] = elemVolVars[scvIdx].massFraction(1);
                rho[dofIdxGlobal]    = elemVolVars[scvIdx].density();
                rhoMol[dofIdxGlobal] = elemVolVars[scvIdx].molarDensity();
                mu[dofIdxGlobal]     = elemVolVars[scvIdx].viscosity();
                (*boxVolume)[dofIdxGlobal]  = fvGeometry.subContVol[scvIdx].volume;
                (*porosity)[dofIdxGlobal]   = elemVolVars[scvIdx].porosity();
                (*S_n_CH4)[dofIdxGlobal]    = elemVolVars[scvIdx].molarDensity()
                                              *elemVolVars[scvIdx].moleFraction(0)
                                              *elemVolVars[scvIdx].porosity()
                                              *fvGeometry.subContVol[scvIdx].volume;
                (*S_Ads_CH4)[dofIdxGlobal]  = elemVolVars[scvIdx].adsorption(0)
                                               *fvGeometry.subContVol[scvIdx].volume;
                (*S_nAds_CH4)[dofIdxGlobal] = (elemVolVars[scvIdx].molarDensity()
                                               *elemVolVars[scvIdx].moleFraction(0)
                                               *elemVolVars[scvIdx].porosity()
                                               +elemVolVars[scvIdx].adsorption(0))
                                               *fvGeometry.subContVol[scvIdx].volume;
                (*S_n_CO2)[dofIdxGlobal]    = elemVolVars[scvIdx].molarDensity()
                                               *elemVolVars[scvIdx].moleFraction(1)
                                               *elemVolVars[scvIdx].porosity()
                                               *fvGeometry.subContVol[scvIdx].volume;
                (*S_Ads_CO2)[dofIdxGlobal]  = elemVolVars[scvIdx].adsorption(1)
                                               *fvGeometry.subContVol[scvIdx].volume;
                (*S_nAds_CO2)[dofIdxGlobal] = (elemVolVars[scvIdx].molarDensity()
                                               *elemVolVars[scvIdx].moleFraction(1)
                                               *elemVolVars[scvIdx].porosity()
                                               +elemVolVars[scvIdx].adsorption(1))
                                               *fvGeometry.subContVol[scvIdx].volume;
                        }

            // velocity output
            velocityOutput.calculateVelocity(*velocity, elemVolVars, fvGeometry, element, phaseIdx);
        }

        writer.attachDofData(pressure, "P", isBox);
        writer.attachDofData(delp, "delp", isBox);
        if (velocityOutput.enableOutput())
        {
            writer.attachDofData(*velocity,  "velocity", isBox, dim);
        }

        writer.attachDofData(moleFraction0, "x_" + FluidSystem::componentName(0), isBox);
        writer.attachDofData(moleFraction1, "x_" + FluidSystem::componentName(1), isBox);
        writer.attachDofData(massFraction0, "X_" + FluidSystem::componentName(0), isBox);
        writer.attachDofData(massFraction1, "X_" + FluidSystem::componentName(1), isBox);

        writer.attachDofData(rho, "rho", isBox);
        writer.attachDofData(rhoMol, "rhoMol", isBox);
        writer.attachDofData(mu, "mu", isBox);
        writer.attachVertexData(*boxVolume, "boxVolume", isBox);
        writer.attachVertexData(*porosity, "porosity", isBox);
        writer.attachDofData(*S_n_CH4, "S_n_CH4", isBox);
        writer.attachDofData(*S_Ads_CH4, "S_Ads_CH4", isBox);
        writer.attachDofData(*S_nAds_CH4, "S_nAds_CH4", isBox);
        writer.attachDofData(*S_n_CO2, "S_n_CO2", isBox);
        writer.attachDofData(*S_Ads_CO2, "S_Ads_CO2", isBox);
        writer.attachDofData(*S_nAds_CO2, "S_nAds_CO2", isBox);
        writer.attachCellData(rank, "process rank");
    }
};
}

#include "propertydefaults.hh"

#endif
