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
 * \brief Base class for all models which use the one-phase two-component linear elasticity model.
 *        Adaption of the fully implicit scheme to the one-phase two-component linear elasticity model.
 */
#ifndef DUMUX_ELASTIC1P2C_MODEL_HH
#define DUMUX_ELASTIC1P2C_MODEL_HH

#include "properties.hh"
#include <dumux/common/eigenvalues.hh>

namespace Dumux {
/*!
 * \ingroup ElOnePTwoCBoxModel
 * \brief Adaption of the fully implicit scheme to the one-phase two-component linear elasticity model.
 *
 * This model implements a one-phase flow of an incompressible fluid, that consists of two components.
 * The deformation of the solid matrix is described with a quasi-stationary momentum balance equation.
 * The influence of the pore fluid is accounted for through the effective stress concept (Biot 1941 \cite Biot1941a).
 * The total stress acting on a rock is partially supported by the rock matrix and partially supported
 * by the pore fluid. The effective stress represents the share of the total stress which is supported
 * by the solid rock matrix and can be determined as a function of the strain according to Hooke's law.
 *
 * As an equation for the conservation of momentum within the fluid phase Darcy's approach is used:
 \f[
 v = - \frac{\textbf K}{\mu}
 \left(\textbf{grad}\, p - \varrho_w {\textbf g} \right)
 \f]
 *
 * Gravity can be enabled or disabled via the property system.
 * By inserting this into the volume balance of the solid-fluid mixture, one gets
 \f[
 \frac{\partial \text{div} \textbf{u}}{\partial t}
 - \text{div} \left\{
   \frac{\textbf K}{\mu}  \left(\textbf{grad}\, p - \varrho_w {\textbf g} \right)
   + \sum_\kappa D^\kappa_\text{pm} \frac{M^\kappa}{M_\alpha} \textbf{grad} x^\kappa
 \right\} = q \;,
 \f]
 *
 * The transport of the components \f$\kappa \in \{ w, a \}\f$ is described by the following equation:
 \f[
 \frac{ \partial \phi_{eff} X^\kappa}{\partial t}
 - \text{div} \left\lbrace
 X^\kappa \frac{{\textbf K}}{\mu} \left( \textbf{grad}\, p - \varrho_w {\textbf g} \right)
 + D^\kappa_\text{pm} \frac{M^\kappa}{M_\alpha} \textbf{grad} x^\kappa
 - \phi_{eff} X^\kappa \frac{\partial \boldsymbol{u}}{\partial t}
 \right\rbrace = q.
 \f]
 *
 * If the model encounters stability problems, a stabilization term can be switched on. The stabilization
 * term is defined in Aguilar et al. (2008) \cite Aguilar2008a.
 \f[
 \beta \text{div} \textbf{grad} \frac{\partial p}{\partial t}
 \f]
 with \f$\beta\f$:
 \f[
 \beta = h^2 / 4(\lambda + 2 \mu)
 \f]
 * where \f$h\f$ is the discretization length.
 *
 * The balance equations
 * with the stabilization term are given below:
 \f[
 \frac{\partial \text{div} \textbf{u}}{\partial t}
 - \text{div} \left\{
   \frac{\textbf K}{\mu}  \left(\textbf{grad}\, p - \varrho_w {\textbf g} \right)
   + \sum_\kappa D^\kappa_\text{pm} \frac{M^\kappa}{M_\alpha} \textbf{grad} x^\kappa
   + \varrho_w \beta \textbf{grad} \frac{\partial p}{\partial t}
 \right\} = q \;,
 \f]
 *
 * The transport of the components \f$\kappa \in \{ w, a \}\f$ is described by the following equation:
 *
 \f[
 \frac{ \partial \phi_{eff} X^\kappa}{\partial t}
 - \text{div} \left\lbrace
 X^\kappa \frac{{\textbf K}}{\mu} \left( \textbf{grad}\, p - \varrho_w {\textbf g} \right)
 + \varrho_w X^\kappa \beta \textbf{grad} \frac{\partial p}{\partial t}
 + D^\kappa_\text{pm} \frac{M^\kappa}{M_\alpha} \textbf{grad} x^\kappa
 - \phi_{eff} X^\kappa \frac{\partial \boldsymbol{u}}{\partial t}
 \right\rbrace = q.
 \f]
 *
 *
 * The quasi-stationary momentum balance equation is:
 \f[
 \text{div}\left( \boldsymbol{\sigma'}- p \boldsymbol{I} \right) + \left( \phi_{eff} \varrho_w + (1 - \phi_{eff}) * \varrho_s  \right)
  {\textbf g} = 0 \;,
 \f]
 * with the effective stress:
 \f[
  \boldsymbol{\sigma'} = 2\,G\,\boldsymbol{\epsilon} + \lambda \,\text{tr} (\boldsymbol{\epsilon}) \, \boldsymbol{I}.
 \f]
 *
 * and the strain tensor \f$\boldsymbol{\epsilon}\f$ as a function of the solid displacement gradient \f$\textbf{grad} \boldsymbol{u}\f$:
 \f[
  \boldsymbol{\epsilon} = \frac{1}{2} \, (\textbf{grad} \boldsymbol{u} + \textbf{grad}^T \boldsymbol{u}).
 \f]
 *
 * Here, the rock mechanics sign convention is switch off which means compressive stresses are < 0 and tensile stresses are > 0.
 * The rock mechanics sign convention can be switched on for the vtk output via the property system.
 *
 * The effective porosity is calculated as a function of the solid displacement:
 \f[
      \phi_{eff} = \frac{\phi_{init} + \text{div} \boldsymbol{u}}{1 + \text{div}}
 \f]
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial
 * and the implicit Euler method as time discretization.
 *
 * The primary variables are the pressure \f$p\f$ and the mole or mass fraction of dissolved component \f$x\f$ and the solid
 * displacement vector \f$\boldsymbol{u}\f$.
 */


template<class TypeTag>
class ElOnePTwoCModel: public GET_PROP_TYPE(TypeTag, BaseModel)
{
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum {
        dim = GridView::dimension
    };

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;


    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dim> DimVector;
    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;

public:
    /*!
     * \brief \copybrief ImplicitModel::addOutputVtkFields
     *
     * Specialization for the ElOnePTwoCBoxModel, add one-phase two-component
     * properties, solid displacement, stresses, effective properties and the
     * process rank to the VTK writer.
     */
    template<class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol, MultiWriter &writer) {

        // check whether compressive stresses are defined to be positive
        // (rockMechanicsSignConvention_ == true) or negative
        rockMechanicsSignConvention_ =  GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, RockMechanicsSignConvention);

        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, dim> > VectorField;

        // create the required scalar and vector fields
        unsigned numVertices = this->gridView_().size(dim);
        unsigned numElements = this->gridView_().size(0);

        // create the required fields for vertex data
        ScalarField &pressure = *writer.allocateManagedBuffer(numVertices);
        ScalarField &moleFraction0 = *writer.allocateManagedBuffer(numVertices);
        ScalarField &moleFraction1 = *writer.allocateManagedBuffer(numVertices);
        ScalarField &massFraction0 = *writer.allocateManagedBuffer(numVertices);
        ScalarField &massFraction1 = *writer.allocateManagedBuffer(numVertices);
        VectorField &displacement = *writer.template allocateManagedBuffer<Scalar, dim>(numVertices);
        ScalarField &density = *writer.allocateManagedBuffer(numVertices);
        ScalarField &viscosity = *writer.allocateManagedBuffer(numVertices);
        ScalarField &porosity = *writer.allocateManagedBuffer(numVertices);
        ScalarField &Kx = *writer.allocateManagedBuffer(numVertices);

        // create the required fields for element data
        // effective stresses
        VectorField &effStressX = *writer.template allocateManagedBuffer<Scalar,
                dim>(numElements);
        VectorField &effStressY = *writer.template allocateManagedBuffer<Scalar,
                dim>(numElements);
        VectorField &effStressZ = *writer.template allocateManagedBuffer<Scalar,
                dim>(numElements);
        // total stresses
        VectorField &totalStressX = *writer.template allocateManagedBuffer<
                Scalar, dim>(numElements);
        VectorField &totalStressY = *writer.template allocateManagedBuffer<
                Scalar, dim>(numElements);
        VectorField &totalStressZ = *writer.template allocateManagedBuffer<
                Scalar, dim>(numElements);

        // principal stresses
        ScalarField &principalStress1 = *writer.allocateManagedBuffer(
                numElements);
        ScalarField &principalStress2 = *writer.allocateManagedBuffer(
                numElements);
        ScalarField &principalStress3 = *writer.allocateManagedBuffer(
                numElements);

        ScalarField &effPorosity = *writer.allocateManagedBuffer(numElements);
        ScalarField &cellPorosity = *writer.allocateManagedBuffer(numElements);
        ScalarField &cellKx = *writer.allocateManagedBuffer(numElements);
        ScalarField &cellPressure = *writer.allocateManagedBuffer(numElements);

        // initialize cell stresses, cell-wise hydraulic parameters and cell pressure with zero
        for (unsigned int eIdx = 0; eIdx < numElements; ++eIdx) {
            effStressX[eIdx] = Scalar(0.0);
            if (dim >= 2)
                effStressY[eIdx] = Scalar(0.0);
            if (dim >= 3)
                effStressZ[eIdx] = Scalar(0.0);

            totalStressX[eIdx] = Scalar(0.0);
            if (dim >= 2)
                totalStressY[eIdx] = Scalar(0.0);
            if (dim >= 3)
                totalStressZ[eIdx] = Scalar(0.0);

            principalStress1[eIdx] = Scalar(0.0);
            if (dim >= 2)
                principalStress2[eIdx] = Scalar(0.0);
            if (dim >= 3)
                principalStress3[eIdx] = Scalar(0.0);

            effPorosity[eIdx] = Scalar(0.0);
            cellPorosity[eIdx] = Scalar(0.0);
            cellKx[eIdx] = Scalar(0.0);
            cellPressure[eIdx] = Scalar(0.0);
        }
        ScalarField &rank = *writer.allocateManagedBuffer(numElements);


        FVElementGeometry fvGeometry;
        ElementVolumeVariables elemVolVars;
        ElementBoundaryTypes elemBcTypes;

        // initialize start and end of element iterator
        // loop over all elements (cells)
        for (const auto& element : elements(this->gridView_(), Dune::Partitions::interior))
        {
            unsigned int eIdx = this->problem_().model().elementMapper().index(element);
            rank[eIdx] = this->gridView_().comm().rank();

            fvGeometry.update(this->gridView_(), element);
            elemBcTypes.update(this->problem_(), element, fvGeometry);
            elemVolVars.update(this->problem_(), element, fvGeometry, false);

            // loop over all local vertices of the cell
            int numScv = element.subEntities(dim);

            for (int scvIdx = 0; scvIdx < numScv; ++scvIdx)
            {
                unsigned int vIdxGlobal = this->dofMapper().subIndex(element, scvIdx, dim);

                pressure[vIdxGlobal] = elemVolVars[scvIdx].pressure();
                moleFraction0[vIdxGlobal] = elemVolVars[scvIdx].moleFraction(0);
                moleFraction1[vIdxGlobal] = elemVolVars[scvIdx].moleFraction(1);
                massFraction0[vIdxGlobal] = elemVolVars[scvIdx].massFraction(0);
                massFraction1[vIdxGlobal] = elemVolVars[scvIdx].massFraction(1);
                // in case of rock mechanics sign convention solid displacement is
                // defined to be negative if it points in positive coordinate direction
                if(rockMechanicsSignConvention_){
                    DimVector tmpDispl;
                    tmpDispl = Scalar(0);
                    tmpDispl -= elemVolVars[scvIdx].displacement();
                    displacement[vIdxGlobal] = tmpDispl;
                    }

                else
                    displacement[vIdxGlobal] = elemVolVars[scvIdx].displacement();

                density[vIdxGlobal] = elemVolVars[scvIdx].density();
                viscosity[vIdxGlobal] = elemVolVars[scvIdx].viscosity();
                porosity[vIdxGlobal] = elemVolVars[scvIdx].porosity();
                Kx[vIdxGlobal] =    this->problem_().spatialParams().intrinsicPermeability(
                                element, fvGeometry, scvIdx)[0][0];
                // calculate cell quantities by adding up scv quantities and dividing through numScv
                cellPorosity[eIdx] += elemVolVars[scvIdx].porosity()    / numScv;
                cellKx[eIdx] += this->problem_().spatialParams().intrinsicPermeability(
                                element, fvGeometry, scvIdx)[0][0] / numScv;
                cellPressure[eIdx] += elemVolVars[scvIdx].pressure()    / numScv;
            };

            // calculate cell quantities for variables which are defined at the integration point
            Scalar tmpEffPoro;
            DimMatrix tmpEffStress;
            tmpEffStress = Scalar(0);
            tmpEffPoro = Scalar(0);

            // loop over all scv-faces of the cell
            for (int fIdx = 0; fIdx < fvGeometry.numScvf; fIdx++) {

                //prepare the flux calculations (set up and prepare geometry, FE gradients)
                FluxVariables fluxVars;
                fluxVars.update(this->problem_(),
                                element, fvGeometry,
                                fIdx,
                                elemVolVars);

                // divide by number of scv-faces and sum up edge values
                tmpEffPoro = fluxVars.effPorosity() / fvGeometry.numScvf;
                tmpEffStress = fluxVars.sigma();
                tmpEffStress /= fvGeometry.numScvf;

                effPorosity[eIdx] += tmpEffPoro;

                // in case of rock mechanics sign convention compressive stresses
                // are defined to be positive
                if(rockMechanicsSignConvention_){
                    effStressX[eIdx] -= tmpEffStress[0];
                    if (dim >= 2) {
                        effStressY[eIdx] -= tmpEffStress[1];
                    }
                    if (dim >= 3) {
                        effStressZ[eIdx] -= tmpEffStress[2];
                    }
                }
                else{
                    effStressX[eIdx] += tmpEffStress[0];
                    if (dim >= 2) {
                        effStressY[eIdx] += tmpEffStress[1];
                    }
                    if (dim >= 3) {
                        effStressZ[eIdx] += tmpEffStress[2];
                    }
                }
            }

            // calculate total stresses
            // in case of rock mechanics sign convention compressive stresses
            // are defined to be positive and total stress is calculated by adding the pore pressure
            if(rockMechanicsSignConvention_){
                totalStressX[eIdx][0] = effStressX[eIdx][0]    + cellPressure[eIdx];
                totalStressX[eIdx][1] = effStressX[eIdx][1];
                totalStressX[eIdx][2] = effStressX[eIdx][2];
                if (dim >= 2) {
                    totalStressY[eIdx][0] = effStressY[eIdx][0];
                    totalStressY[eIdx][1] = effStressY[eIdx][1]    + cellPressure[eIdx];
                    totalStressY[eIdx][2] = effStressY[eIdx][2];
                }
                if (dim >= 3) {
                    totalStressZ[eIdx][0] = effStressZ[eIdx][0];
                    totalStressZ[eIdx][1] = effStressZ[eIdx][1];
                    totalStressZ[eIdx][2] = effStressZ[eIdx][2]    + cellPressure[eIdx];
                }
            }
            else{
                totalStressX[eIdx][0] = effStressX[eIdx][0]    - cellPressure[eIdx];
                totalStressX[eIdx][1] = effStressX[eIdx][1];
                totalStressX[eIdx][2] = effStressX[eIdx][2];
                if (dim >= 2) {
                    totalStressY[eIdx][0] = effStressY[eIdx][0];
                    totalStressY[eIdx][1] = effStressY[eIdx][1]    - cellPressure[eIdx];
                    totalStressY[eIdx][2] = effStressY[eIdx][2];
                }
                if (dim >= 3) {
                    totalStressZ[eIdx][0] = effStressZ[eIdx][0];
                    totalStressZ[eIdx][1] = effStressZ[eIdx][1];
                    totalStressZ[eIdx][2] = effStressZ[eIdx][2]    - cellPressure[eIdx];
                }
            }
        }

        // calculate principal stresses i.e. the eigenvalues of the total stress tensor
        Scalar a1, a2, a3;
        DimMatrix totalStress;
        DimVector eigenValues;
        a1=Scalar(0);
        a2=Scalar(0);
        a3=Scalar(0);

        for (unsigned int eIdx = 0; eIdx < numElements; eIdx++)
        {
            eigenValues = Scalar(0);
            totalStress = Scalar(0);

            totalStress[0] = totalStressX[eIdx];
            if (dim >= 2)
                totalStress[1] = totalStressY[eIdx];
            if (dim >= 3)
                totalStress[2] = totalStressZ[eIdx];

            calculateEigenValues<dim>(eigenValues, totalStress);


            for (int i = 0; i < dim; i++)
                {
                    if (std::isnan(eigenValues[i]))
                        eigenValues[i] = 0.0;
                }

            // sort principal stresses: principalStress1 >= principalStress2 >= principalStress3
            if (dim == 2) {
                a1 = eigenValues[0];
                a2 = eigenValues[1];

                if (a1 >= a2) {
                    principalStress1[eIdx] = a1;
                    principalStress2[eIdx] = a2;
                } else {
                    principalStress1[eIdx] = a2;
                    principalStress2[eIdx] = a1;
                }
            }

            if (dim == 3) {
                a1 = eigenValues[0];
                a2 = eigenValues[1];
                a3 = eigenValues[2];

                if (a1 >= a2) {
                    if (a1 >= a3) {
                        principalStress1[eIdx] = a1;
                        if (a2 >= a3) {
                            principalStress2[eIdx] = a2;
                            principalStress3[eIdx] = a3;
                        }
                        else //a3 > a2
                        {
                            principalStress2[eIdx] = a3;
                            principalStress3[eIdx] = a2;
                        }
                    }
                    else // a3 > a1
                    {
                        principalStress1[eIdx] = a3;
                        principalStress2[eIdx] = a1;
                        principalStress3[eIdx] = a2;
                    }
                } else // a2>a1
                {
                    if (a2 >= a3) {
                        principalStress1[eIdx] = a2;
                        if (a1 >= a3) {
                            principalStress2[eIdx] = a1;
                            principalStress3[eIdx] = a3;
                        }
                        else //a3>a1
                        {
                            principalStress2[eIdx] = a3;
                            principalStress3[eIdx] = a1;
                        }
                    }
                    else //a3>a2
                    {
                        principalStress1[eIdx] = a3;
                        principalStress2[eIdx] = a2;
                        principalStress3[eIdx] = a1;
                    }
                }
            }

        }

        writer.attachVertexData(pressure, "P");

        writer.attachVertexData(moleFraction0, "x_" + FluidSystem::componentName(0));
        writer.attachVertexData(moleFraction1, "x_" + FluidSystem::componentName(1));
        writer.attachVertexData(massFraction0, "X_" + FluidSystem::componentName(0));
        writer.attachVertexData(massFraction1, "X_" + FluidSystem::componentName(1));

        writer.attachVertexData(displacement, "u", dim);
        writer.attachVertexData(density, "rho");
        writer.attachVertexData(viscosity, "mu");
        writer.attachVertexData(porosity, "porosity");
        writer.attachVertexData(Kx, "Kx");
        writer.attachCellData(cellPorosity, "porosity");
        writer.attachCellData(cellKx, "Kx");
        writer.attachCellData(effPorosity, "effective porosity");

        writer.attachCellData(totalStressX, "total stresses X", dim);
        if (dim >= 2)
            writer.attachCellData(totalStressY, "total stresses Y", dim);
        if (dim >= 3)
            writer.attachCellData(totalStressZ, "total stresses Z", dim);

        writer.attachCellData(effStressX, "effective stress changes X", dim);
        if (dim >= 2)
            writer.attachCellData(effStressY, "effective stress changes Y",    dim);
        if (dim >= 3)
            writer.attachCellData(effStressZ, "effective stress changes Z",    dim);

        writer.attachCellData(principalStress1, "principal stress 1");
        if (dim >= 2)
            writer.attachCellData(principalStress2, "principal stress 2");
        if (dim >= 3)
            writer.attachCellData(principalStress3, "principal stress 3");

        writer.attachCellData(cellPressure, "P");

        writer.attachCellData(rank, "rank");

    }
private:
    bool rockMechanicsSignConvention_;

};
}
#include "propertydefaults.hh"
#endif
