/*****************************************************************************
 *   Copyright (C) 2009-2010 by Andreas Lauser                               *
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
#ifndef DUMUX_MPNC_MODEL_HH
#define DUMUX_MPNC_MODEL_HH

#include "MpNcproperties.hh"
#include "MpNcvtkwriter.hh"

#include <dumux/boxmodels/common/boxmodel.hh>
#include <tr1/array>

namespace Dumux
{
/*!
 * \ingroup MPNCModel
 * \brief Adaption of the box scheme to compositional twophase flows.
 *
 * This model implements a two-phase flow of a fluid mixture composed
 * of \f$N\f$ chemical species. The phases are denoted by \f$\alpha
 * \in \{ l, g \}\f$ for liquid and gas. The liquid and the gas phase
 * both are a mixture of of \f$N \geq 2\f$ species. The model assumes,
 * a fluid configuration of a solvent species and \f$N-1\f$ solute
 * species with very low solubility.
 *
 * The standard multiphase Darcy approach is used as the equation for
 * the conservation of momentum:
 * \f[
     v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \boldsymbol{K}
     \left(\text{grad} p_\alpha - \varrho_{\alpha} \boldsymbol{g} \right)
 * \f]
 *
 * By inserting this into the equations for the conservation of the
 * components, one gets one transport equation for each component \f$\kappa\f$
 * \f{eqnarray*}
    && \sum_{\kappa \in \alpha} \left(
    %
    \phi \frac{\partial \varrho_\alpha x_{\alpha\kappa} S_\alpha}{\partial t}
    -
    \nabla \cdot
    \left\{
       \frac{\varrho_\alpha}{\overline M_\alpha} x_{\alpha\kappa}
       \frac{k_{r\alpha}}{\mu_\alpha} \boldsymbol{K}
       (\nabla p_\alpha - \varrho_{\alpha} \boldsymbol{g})
    \right\}
    \left)
    \nonumber \\
    \nonumber \
    &-& \sum_{\kappa \in \ \nabla \cdot \left\{{\bf D_{pm}^\kappa} \frac{\varrho_{\alpha}}{\overline M_\alpha} {\bf \nabla} x^\kappa_{\alpha} \right\}
    - \sum_{\kappa \in \alpha} q_\alpha^\kappa = \quad 0 \qquad \kappa \in \{1, \dots, N\} \, ,
    \alpha \in \{l, g\}
    \f}
 * with \f$\overline M_\alpha\f$ being the average molar mass of the phase \f$\alpha\f$:
 * \f[ \overline M_\alpha = \sum_{\kappa = 1}^N M_\kappa \; x_{\alpha\kappa} \f]
 *
 * This is discretized in the model using the fully-coupled vertex
 * centered finite volume (box) scheme as spatial and
 * the implicit Euler method as temporal discretization.
 *
 * The model uses \f$x_{l1}, \dots, x_{lN}, S_g, p_g \f$ as primary
 * variables. \f$x_{g\kappa}\f$ is calculated using Henry's law (if
 * the componentis a solute), or the vapor pressure (if \f$\kappa\f$
 * is the solvent) as \f$\beta_\kappa\f$, the constant of
 * proportionality between the liquid mole fraction and the partial
 * pressure of a component:
 * \f[ x_{g\kappa} = \frac{\beta_kappa x_{l\kappa}}{\sum_{i=1}^N x_{l\kappa} \beta_kappa} \f]
 *
 * Additionally two auxiliary conditions are used to keep the solution physical.
 * \todo describe NCP approach
 */
template<class TypeTag>
class MPNCModel : public BoxModel<TypeTag>
{
    typedef MPNCModel<TypeTag> ThisType;
    typedef BoxModel<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVariables)) FluxVariables;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementBoundaryTypes)) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MPNCIndices)) Indices;
    typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

    typedef Dumux::MPNCVtkWriter<TypeTag> MPNCVtkWriter;


    enum {
        enableEnergy = GET_PROP_VALUE(TypeTag, PTAG(EnableEnergy)),
        enableDiffusion = GET_PROP_VALUE(TypeTag, PTAG(EnableDiffusion)),
        enableKinetic = GET_PROP_VALUE(TypeTag, PTAG(EnableKinetic)),
        enableKineticEnergy = GET_PROP_VALUE(TypeTag, PTAG(EnableKineticEnergy)),
        enableSmoothUpwinding = GET_PROP_VALUE(TypeTag, PTAG(EnableSmoothUpwinding)),
        enablePartialReassemble = GET_PROP_VALUE(TypeTag, PTAG(EnablePartialReassemble)),
        enableJacobianRecycling = GET_PROP_VALUE(TypeTag, PTAG(EnableJacobianRecycling)),
        numDiffMethod = GET_PROP_VALUE(TypeTag, PTAG(NumericDifferenceMethod)),
        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),
        numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)),
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
        dimWorld = GridView::dimensionworld,
        dim = GridView::dimension
    };

public:
    ~MPNCModel()
    { delete vtkWriter_; };

    void init(Problem &problem)
    {
        ParentType::init(problem);
        vtkWriter_ = new MPNCVtkWriter(problem);

        if (this->gridView_().comm().rank() == 0)
            std::cout
                << "Initializing M-phase N-component model: \n"
                << "    phases: " << numPhases << "\n"
                << "    components: " << numComponents << "\n"
                << "    equations: " << numEq << "\n"
                << "    kinetic mass transfer: " << enableKinetic<< "\n"
                << "    kinetic energy transfer: " << enableKineticEnergy<< "\n"
                << "    diffusion: " << enableDiffusion << "\n"
                << "    energy equation: " << enableEnergy << "\n"
                << "    smooth upwinding: " << enableSmoothUpwinding << "\n"
                << "    partial jacobian reassembly: " << enablePartialReassemble << "\n"
                << "    numeric differentiation method: " << numDiffMethod << " (-1: backward, 0: central, +1 forward)\n"
                << "    jacobian recycling: " << enableJacobianRecycling << "\n";
    }

    /*!
     * \brief Compute the total storage inside one phase of all
     *        conservation quantities.
     */
    void globalPhaseStorage(PrimaryVariables &dest, int phaseIdx)
    {
        dest = 0;

        ElementIterator elemIt = this->gridView_().template begin<0>();
        const ElementIterator elemEndIt = this->gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            this->localResidual().addPhaseStorage(dest, *elemIt, phaseIdx);
        };

        this->gridView_().comm().sum(dest);
    }

    /*!
     * \brief Add the result of the current timestep to the VTK output.
     */
    template <class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {
        vtkWriter_->addCurrentSolution(writer);
    }

    MPNCVtkWriter *vtkWriter_;
};

}

#include "MpNcpropertydefaults.hh"

#endif
