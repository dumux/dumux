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

#ifndef DUMUX_1PNCMIN_MODEL_HH
#define DUMUX_1PNCMIN_MODEL_HH

#include "properties.hh"
#include "indices.hh"
#include "localresidual.hh"

// #include <dumux/material/constants.hh>
#include <dumux/porousmediumflow/1pnc/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/velocityoutput.hh>

namespace Dumux
{
/*!
 * \ingroup OnePNCMinModel
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
class OnePNCMinModel: public OnePNCModel<TypeTag>
{
    typedef Dumux::OnePNCMinModel<TypeTag> ThisType;
    typedef Dumux::OnePNCModel<TypeTag> ParentType;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    //old
//     typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
//     typedef Dumux::Constants<Scalar> Constant;
//     typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

//         numEq = GET_PROP_VALUE(TypeTag, NumEq),
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numSPhases = GET_PROP_VALUE(TypeTag, NumSPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),
        numSComponents = FluidSystem::numSComponents,

        pressureIdx = Indices::pressureIdx,
        firstMoleFracIdx = Indices::firstMoleFracIdx,

        phaseIdx = Indices::phaseIdx,
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using CoordScalar = typename GridView::ctype;
    using Tensor = Dune::FieldMatrix<CoordScalar, dimWorld, dimWorld>;

    //old
//     typedef typename GridView::template Codim<dim>::Entity Vertex;
//     typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:

    /*!
     * \brief Apply the initial conditions to the model.
     *
     * \param problem The object representing the problem which needs to
     *             be simulated.
     */
    void init(Problem& problem)
    {
        ParentType::init(problem);

       // register standardized vtk output fields
        auto& vtkOutputModule = problem.vtkOutputModule();
        vtkOutputModule.addSecondaryVariable("p", [](const VolumeVariables& v){ return v.pressure(phaseIdx); });
        vtkOutputModule.addSecondaryVariable("pw", [](const VolumeVariables& v){ return v.pressure(phaseIdx); });
        vtkOutputModule.addSecondaryVariable("rho", [](const VolumeVariables& v){ return v.density(phaseIdx); });
        vtkOutputModule.addSecondaryVariable("porosity", [](const VolumeVariables& v){ return v.porosity(); });
        vtkOutputModule.addSecondaryVariable("permeabilityFactor", [](const VolumeVariables& v){ return v.permeabilityFactor(); });

        vtkOutputModule.addSecondaryVariable("temperature", [](const VolumeVariables& v){ return v.temperature(); });

        for (int sPhaseIdx = 0; sPhaseIdx < numSPhases; ++sPhaseIdx)
            vtkOutputModule.addSecondaryVariable("precipitateVolumeFraction_" + FluidSystem::phaseName(numPhases + sPhaseIdx),
                                                 [sPhaseIdx](const VolumeVariables& v)
                                                 { return v.precipitateVolumeFraction(numPhases + sPhaseIdx); });

        vtkOutputModule.addSecondaryVariable("Kxx",
                                             [this](const VolumeVariables& v){ return this->perm_(v.permeability())[0][0]; });
        if (dim >= 2)
            vtkOutputModule.addSecondaryVariable("Kyy",
                                                 [this](const VolumeVariables& v){ return this->perm_(v.permeability())[1][1]; });
        if (dim >= 3)
            vtkOutputModule.addSecondaryVariable("Kzz",
                                                 [this](const VolumeVariables& v){ return this->perm_(v.permeability())[2][2]; });

        for (int i = 0; i < numComponents; ++i)
                vtkOutputModule.addSecondaryVariable("x_" + FluidSystem::componentName(i),[i](const VolumeVariables& v){ return     v.moleFraction(phaseIdx,i); });

//         for (int i = 0; i < numComponents; ++i)
//            vtkOutputModule.addSecondaryVariable("m^w_" + FluidSystem::componentName(i),
//                                                  [i](const VolumeVariables& v){ return v.molarity(phaseIdx,i); });
    }

    /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     *
     * \param sol The solution vector
     * \param writer The writer for multi-file VTK datasets
     */
//     template<class MultiWriter>
//     //additional output of the permeability and the precipitate volume fractions
//     void addOutputVtkFields(const SolutionVector &sol, MultiWriter &writer)
//     {
//         typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;
//         typedef Dune::BlockVector<Dune::FieldVector<double, dim> > VectorField;
//
//         // get the number of degrees of freedom
//         auto numDofs = this->numDofs();
//         ScalarField *p           = writer.allocateManagedBuffer (numDofs);
//         ScalarField *rho         = writer.allocateManagedBuffer (numDofs);
//         ScalarField *temperature  = writer.allocateManagedBuffer (numDofs);
//         ScalarField *poro         = writer.allocateManagedBuffer (numDofs);
//         ScalarField *permeabilityFactor = writer.allocateManagedBuffer (numDofs);
//         ScalarField *precipitateVolumeFraction[numSPhases];
//
//         for (int i = 0; i < numSPhases; ++i)
//             precipitateVolumeFraction[i] = writer.allocateManagedBuffer(numDofs);
//
// //         ScalarField *massFraction[numPhases][numComponents];
// //         for (int i = 0; i < numPhases; ++i)
// //             for (int j = 0; j < numComponents; ++j)
// //                 massFraction[i][j] = writer.allocateManagedBuffer(numDofs);
//
//         ScalarField *molarity[numComponents];
//         for (int j = 0; j < numComponents ; ++j)
//             molarity[j] = writer.allocateManagedBuffer(numDofs);
//
//         ScalarField *moleFraction[numComponents];
//             for (int j = 0; j < numComponents; ++j)
//                 moleFraction[j] = writer.allocateManagedBuffer(numDofs);
//
//        ScalarField *moleFractionSolid[numSComponents];
//             for (int j = numComponents; j < numSComponents; ++j)
//                 moleFractionSolid[j] = writer.allocateManagedBuffer(numDofs);
//
//         ScalarField *Perm[dim];
//         for (int j = 0; j < dim; ++j) //Permeability only in main directions xx and yy
//             Perm[j] = writer.allocateManagedBuffer(numDofs);
//
//         VectorField *velocity = writer.template allocateManagedBuffer<double, dim>(numDofs);
//
//         ImplicitVelocityOutput<TypeTag> velocityOutput(this->problem_());
//
//         if (velocityOutput.enableOutput()) // check if velocity output is demanded
//         {
//             // initialize velocity fields
//             for (unsigned int i = 0; i < numDofs; ++i)
//             {
//                 (*velocity)[i] = Scalar(0);
//             }
//         }
//
//         auto numElements = this->gridView_().size(0);
//         ScalarField *rank = writer.allocateManagedBuffer(numElements);
//
//         for (const auto& element : elements(this->gridView_()))
//         {
//             auto eIdxGlobal = this->problem_().elementMapper().index(element);
//             (*rank)[eIdxGlobal] = this->gridView_().comm().rank();
//             FVElementGeometry fvGeometry;
//             fvGeometry.update(this->gridView_(), element);
//
//             ElementVolumeVariables elemVolVars;
//             elemVolVars.update(this->problem_(),
//                                element,
//                                fvGeometry,
//                                false /* oldSol? */);
//
//             for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
//             {
//                 auto dofIdxGlobal = this->dofMapper().subIndex(element, scvIdx, dofCodim);
//
//                 (*p)[dofIdxGlobal] = elemVolVars[scvIdx].pressure(phaseIdx);
//                 (*rho)[dofIdxGlobal] = elemVolVars[scvIdx].density(phaseIdx);
//                 (*poro)[dofIdxGlobal] = elemVolVars[scvIdx].porosity();
//
//                 for (int sPhaseIdx = 0; sPhaseIdx < numSPhases; ++sPhaseIdx){
//                     (*precipitateVolumeFraction[sPhaseIdx])[dofIdxGlobal] = elemVolVars[scvIdx].precipitateVolumeFraction(sPhaseIdx + numPhases);
//                  }
//
//                 (*temperature)[dofIdxGlobal] = elemVolVars[scvIdx].temperature();
//                 (*permeabilityFactor)[dofIdxGlobal] = elemVolVars[scvIdx].permeabilityFactor();
//
// //                     for (int compIdx = 0; compIdx < numComponents; ++compIdx)
// //                         (*massFraction[phaseIdx][compIdx])[dofIdxGlobal]= elemVolVars[scvIdx].massFraction(phaseIdx,compIdx);
//
//                 for (int compIdx = 0; compIdx < numComponents; ++compIdx)
//                     (*molarity[compIdx])[dofIdxGlobal] = (elemVolVars[scvIdx].molarity(compIdx));
//
//                 for (int compIdx = 0; compIdx < numComponents; ++compIdx)
//                     (*moleFraction[compIdx])[dofIdxGlobal]= elemVolVars[scvIdx].moleFraction(compIdx);
//
//                  for (int compIdx = numComponents; compIdx < numSComponents; ++compIdx)
//                     (*moleFractionSolid[compIdx])[dofIdxGlobal]= elemVolVars[scvIdx].moleFraction(compIdx);
//
//                 Tensor K = this->perm_(this->problem_().spatialParams().intrinsicPermeability(element, fvGeometry, scvIdx));
//
//                 for (int j = 0; j<dim; ++j)
//                     (*Perm[j])[dofIdxGlobal] = K[j][j] * elemVolVars[scvIdx].permeabilityFactor();
//             };
//
//             // velocity output
//             if(velocityOutput.enableOutput()){
//                 velocityOutput.calculateVelocity(*velocity, elemVolVars, fvGeometry, element, phaseIdx);
//             }
//         }  // loop over element
//
//         writer.attachDofData(*p, "pressure", isBox);
//         writer.attachDofData(*rho, "rho", isBox);
//         writer.attachDofData(*poro, "porosity", isBox);
//         writer.attachDofData(*permeabilityFactor, "permeabilityFactor", isBox);
//         writer.attachDofData(*temperature, "temperature", isBox);
//
//         for (int i = 0; i < numSPhases; ++i)
//         {
//             std::cout << "test1" <<  "\n";
//             std::ostringstream oss;
//             oss << "precipitateVolumeFraction_" << FluidSystem::componentName(numComponents + i);
//             writer.attachDofData(*precipitateVolumeFraction[i], oss.str(), isBox);
//         }
//
//         writer.attachDofData(*Perm[0], "Kxx", isBox);
//         if (dim >= 2)
//             writer.attachDofData(*Perm[1], "Kyy", isBox);
//         if (dim == 3)
//             writer.attachDofData(*Perm[2], "Kzz", isBox);
//
//
// //             for (int j = 0; j < numComponents; ++j)
// //             {
// //                 std::ostringstream oss;
// //                 oss << "X^" << FluidSystem::phaseName(i) << "_" << FluidSystem::componentName(j);
// //                 writer.attachDofData(*massFraction[i][j], oss.str(), isBox);
// //             }
//
//         for (int j = 0; j < numComponents; ++j)
//         {
//             std::ostringstream oss;
//             oss << "m^w_" << FluidSystem::componentName(j);
//             writer.attachDofData(*molarity[j], oss.str(), isBox);
//         }
//
//         for (int j = 0; j < numComponents; ++j)
//             {
//                 std::ostringstream oss;
//                 oss << "x"
//                     << FluidSystem::componentName(j);
//                 writer.attachDofData(*moleFraction[j], oss.str(), isBox);
//             }
//
//         for (int j = numComponents; j < numSComponents; ++j)
//         {
//             std::ostringstream oss;
//             oss << "m^w_" << FluidSystem::componentName(j);
//             writer.attachDofData(*molarity[j], oss.str(), isBox);
//         }
//
//         for (int j = numComponents; j < numSComponents; ++j)
//             {
//                 std::ostringstream oss;
//                 oss << "x"
//                     << FluidSystem::componentName(j);
//                 writer.attachDofData(*moleFraction[j], oss.str(), isBox);
//             }
//
//         if (velocityOutput.enableOutput()) // check if velocity output is demanded
//         {
//             writer.attachDofData(*velocity,  "velocity", isBox, dim);
//         }
//
//         writer.attachCellData(*rank, "process rank");
//     }

    /*!
     * \brief Update the static data of all vertices in the grid.
     *
     * \param curGlobalSol The current global solution
     * \param oldGlobalSol The previous global solution
     */
//     void updateStaticData(SolutionVector &curGlobalSol,
//                           const SolutionVector &oldGlobalSol)
//     {
//         for (unsigned i = 0; i < this->staticDat_.size(); ++i)
//             this->staticDat_[i].visited = false;
//
//         for (const auto& element : elements(this->gridView_()))
//         {
//             FVElementGeometry fvGeometry;
//             fvGeometry.update(this->gridView_(), element);
//             for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
//             {
//                 auto dofIdxGlobal = this->dofMapper().subIndex(element, scvIdx, dofCodim);
//
//                 if (this->staticDat_[dofIdxGlobal].visited)
//                     continue;
//
//                 this->staticDat_[dofIdxGlobal].visited = true;
//                 VolumeVariables volVars;
//                 volVars.update(curGlobalSol[dofIdxGlobal],
//                                this->problem_(),
//                                element,
//                                fvGeometry,
//                                scvIdx,
//                                false);
//             }
//         }
//     }

       /*!
     * \brief Write the current solution to a restart file.
     *
     * \param outStream The output stream of one entity for the restart file
     * \param entity The entity, either a vertex or an element
     */
    template<class Entity>
    void serializeEntity(std::ostream &outStream, const Entity &entity)
    {
        // write primary variables
        ParentType::serializeEntity(outStream, entity);

        int dofIdxGlobal = this->dofMapper().index(entity);

        if (!outStream.good())
            DUNE_THROW(Dune::IOError, "Could not serialize entity " << dofIdxGlobal);

    }

    /*!
     * \brief Reads the current solution from a restart file.
     *
     * \param inStream The input stream of one entity from the restart file
     * \param entity The entity, either a vertex or an element
     */
    template<class Entity>
    void deserializeEntity(std::istream &inStream, const Entity &entity)
    {
        // read primary variables
        ParentType::deserializeEntity(inStream, entity);

        // read phase presence
        int dofIdxGlobal = this->dofMapper().index(entity);

        if (!inStream.good())
            DUNE_THROW(Dune::IOError, "Could not deserialize entity " << dofIdxGlobal);
    }


    Tensor perm_(Scalar perm) const
    {
        Tensor K(0.0);

        for(int i=0; i<dim; i++)
            K[i][i] = perm;

       return K;
    }

    const Tensor& perm_(const Tensor& perm) const
    {
       return perm;
    }
};

}

#include "propertydefaults.hh"

#endif
