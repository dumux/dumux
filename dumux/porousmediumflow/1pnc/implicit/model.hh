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
* \brief  Base class for all models which use the single-phase, n-component fully implicit model.
*         Adaption of the fully implicit model to the one-phase n-component flow model.
*/

#ifndef DUMUX_1PNC_MODEL_HH
#define DUMUX_1PNC_MODEL_HH

#include <dumux/porousmediumflow/implicit/velocityoutput.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/model.hh>

#include "properties.hh"
#include "indices.hh"
// #include "localresidual.hh"

namespace Dumux
{
/*!
 * \ingroup OnePNCModel
 * \brief Adaption of the fully implicit scheme to the
 *        one-phase n-component flow model.
 *
* This model implements a one-phase flow of a compressible fluid, that consists of n components,
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
 \phi\frac{\partial \varrho}{\partial t} - \text{div} \left\{
   \varrho \frac{\textbf K}{\mu}  \left(\textbf{grad}\, p - \varrho {\textbf g} \right)
 \right\} = q \;,
 \f]
 *
 * The transport of the components \f$\kappa \in \{ w, a, ... \}\f$ is described by the following equation:
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

template<class TypeTag>
class OnePNCModel: public GET_PROP_TYPE(TypeTag, BaseModel)
{
    using ParentType = typename GET_PROP_TYPE(TypeTag, BaseModel);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using NonIsothermalModel = Dumux::NonIsothermalModel<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
//    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

    static const int phaseIdx = Indices::phaseIdx;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };
//
//     enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
//
    enum {  numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
//
//      enum {
//             pressureIdx = Indices::pressureIdx,
//             firstMoleFracIdx = Indices::firstMoleFracIdx,
//     };

    typedef typename GridView::template Codim<dim>::Entity Vertex;
    using Element = typename GridView::template Codim<0>::Entity;

//     using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using CoordScalar = typename GridView::ctype;
    using Tensor = Dune::FieldMatrix<CoordScalar, dimWorld, dimWorld>;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:
    /*!
     * \brief Apply the initial conditions to the model.
     *
     * \param problem The object representing the problem which needs to
     *             be simulated.
     */
    void init(Problem &problem)
    {
        ParentType::init(problem);

        // register standardized vtk output fields

        auto& vtkOutputModule = problem.vtkOutputModule();
        vtkOutputModule.addSecondaryVariable("p", [](const VolumeVariables& v){ return v.pressure(phaseIdx); });
        vtkOutputModule.addSecondaryVariable("rho", [](const VolumeVariables& v){ return v.density(phaseIdx); });
        vtkOutputModule.addSecondaryVariable("porosity", [](const VolumeVariables& v){ return v.porosity(); });
        vtkOutputModule.addSecondaryVariable("temperature", [](const VolumeVariables& v){ return v.temperature(); });

        vtkOutputModule.addSecondaryVariable("Kxx",
                                             [this](const VolumeVariables& v){ return this->perm_(v.permeability())[0][0]; });
        if (dim >= 2)
            vtkOutputModule.addSecondaryVariable("Kyy",
                                                 [this](const VolumeVariables& v){ return this->perm_(v.permeability())[1][1]; });
        if (dim >= 3)
            vtkOutputModule.addSecondaryVariable("Kzz",
                                                 [this](const VolumeVariables& v){ return this->perm_(v.permeability())[2][2]; });

       for (int i = 0; i < numComponents; ++i)
           vtkOutputModule.addSecondaryVariable("x" + FluidSystem::componentName(i),
                                                [i](const VolumeVariables& v){ return v.moleFraction(phaseIdx, i); });

       for (int i = 0; i < numComponents; ++i)
           vtkOutputModule.addSecondaryVariable("m^w_" + FluidSystem::componentName(i),
                                                 [i](const VolumeVariables& v){ return v.molarity(phaseIdx,i); });

        NonIsothermalModel::maybeAddTemperature(vtkOutputModule);
    }

    /*!
     * \brief Compute the total storage of all conservation quantities
     *
     * \param storage Contains the storage of each component
     * \param phaseIdx The phase index
     */
//     void globalPhaseStorage(PrimaryVariables &storage, int phaseIdx)
//     {
//         storage = 0;
//         for (const auto& element : elements(this->gridView_()))
//         {
//             if(element.partitionType() == Dune::InteriorEntity)
//             {
//                 this->localResidual().evalPhaseStorage(element, phaseIdx);
//
//                 for (unsigned int i = 0; i < this->localResidual().storageTerm().size(); ++i)
//                     storage += this->localResidual().storageTerm()[i];
//             }
//         }
//         this->gridView_().comm().sum(storage);
//     }

    /*!
     * \brief Called by the update() method if applying the newton
     *         method was unsuccessful.
     */
    void updateFailed()
    {
        ParentType::updateFailed();
    }

    /*!
     * \brief Called by the problem if a time integration was
     *        successful, post processing of the solution is done and the
     *        result has been written to disk.
     *
     * This should prepare the model for the next time integration.
     */
//     void advanceTimeLevel()
//     {
//         ParentType::advanceTimeLevel();
//     }

    /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     *
     * \param sol The solution vector
     * \param writer The writer for multi-file VTK datasets
     */
//     template<class MultiWriter>
//     void addOutputVtkFields(const SolutionVector &sol,
//                             MultiWriter &writer)
//     {
//
//         typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;
//         typedef Dune::BlockVector<Dune::FieldVector<Scalar, dim> > VectorField;
//
//         // get the number of degrees of freedom
//         auto numDofs = this->numDofs();
//
//         ScalarField *pressure            = writer.allocateManagedBuffer (numDofs);
//         ScalarField *rho          = writer.allocateManagedBuffer (numDofs);
//         ScalarField *temperature   = writer.allocateManagedBuffer (numDofs);
//         ScalarField *poro          = writer.allocateManagedBuffer (numDofs);
//         VectorField *velocity = writer.template allocateManagedBuffer<double, dim>(numDofs);
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
//         ScalarField *moleFraction[numComponents];
// //         for (int i = 0; i < numPhases; ++i)
//             for (int j = 0; j < numComponents; ++j)
//                 moleFraction[j] = writer.allocateManagedBuffer(numDofs);
//
//         ScalarField *molarity[numComponents];
//         for (int j = 0; j < numComponents ; ++j)
//             molarity[j] = writer.allocateManagedBuffer(numDofs);
//
//         ScalarField *Perm[dim];
//         for (int j = 0; j < dim; ++j) //Permeability only in main directions xx and yy
//             Perm[j] = writer.allocateManagedBuffer(numDofs);
//
//         auto numElements = this->gridView_().size(0);
//         ScalarField *rank = writer.allocateManagedBuffer (numElements);
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
//                 GlobalPosition globalPos = fvGeometry.subContVol[scvIdx].global;
//                 (*pressure)[dofIdxGlobal]             = elemVolVars[scvIdx].pressure(phaseIdx);
//                 (*rho)[dofIdxGlobal]           = elemVolVars[scvIdx].density(phaseIdx);
//                 (*poro)[dofIdxGlobal]           = elemVolVars[scvIdx].porosity();
//                 (*temperature)[dofIdxGlobal]    = elemVolVars[scvIdx].temperature();
//
//                 for (int compIdx = 0; compIdx < numComponents; ++compIdx)
//                     (*moleFraction[compIdx])[dofIdxGlobal]= elemVolVars[scvIdx].moleFraction(compIdx);
//
//                 for (int compIdx = 0; compIdx < numComponents; ++compIdx)
//                     (*molarity[compIdx])[dofIdxGlobal] = (elemVolVars[scvIdx].molarity(compIdx));
//
//                 Tensor K = perm_(this->problem_().spatialParams().intrinsicPermeability(element, fvGeometry, scvIdx));
//
//                 for (int j = 0; j<dim; ++j){
//                     (*Perm[j])[dofIdxGlobal] = K[j][j];
//                 }
//             };
//
//             // velocity output
//             if(velocityOutput.enableOutput()){
//                 velocityOutput.calculateVelocity(*velocity, elemVolVars, fvGeometry, element, phaseIdx);
//             }
//
//         } // loop over element
//
//         writer.attachDofData(*pressure, "p", isBox);
//         writer.attachDofData(*rho, "rho", isBox);
//         writer.attachDofData(*poro, "porosity", isBox);
//         writer.attachDofData(*temperature, "temperature", isBox);
//         writer.attachDofData(*Perm[0], "Kxx", isBox);
//         if (dim >= 2)
//             writer.attachDofData(*Perm[1], "Kyy", isBox);
//         if (dim == 3)
//             writer.attachDofData(*Perm[2], "Kzz", isBox);
//
//             for (int j = 0; j < numComponents; ++j)
//             {
//                 std::ostringstream oss;
//                 oss << "x"
//                     << FluidSystem::componentName(j);
// //                     << FluidSystem::phaseName(i);
//                 writer.attachDofData(*moleFraction[j], oss.str(), isBox);
//             }
//
//         for (int j = 0; j < numComponents; ++j)
//         {
//             std::ostringstream oss;
//             oss << "m^w_"
//                 << FluidSystem::componentName(j);
//             writer.attachDofData(*molarity[j], oss.str().c_str(), isBox);
//         }
//
//         if (velocityOutput.enableOutput()) // check if velocity output is demanded
//         {
//             writer.attachDofData(*velocity,  "velocity", isBox, dim);
//         }
//
//         writer.attachCellData(*rank, "process rank");
//     }



private :
    Tensor perm_(Scalar perm)
    {
        Tensor K(0.0);

        for(int i=0; i<dim; i++)
            K[i][i] = perm;

       return K;
    }

    Tensor perm_(Tensor perm)
    {
       return perm;
    }

};

}

#include "propertydefaults.hh"

#endif
