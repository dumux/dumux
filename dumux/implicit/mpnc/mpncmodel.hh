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
 * \brief A fully implicit model for MpNc flow using
 *        vertex centered finite volumes.
 *
 */

#ifndef DUMUX_MPNC_MODEL_HH
#define DUMUX_MPNC_MODEL_HH

#include "mpncproperties.hh"
#include "mpncvtkwriter.hh"

#include <dumux/implicit/common/implicitmodel.hh>

namespace Dumux
{
/*!
 * \ingroup MPNCModel
 * \brief A fully implicit model for MpNc flow using
 *        vertex centered finite volumes.
 *
 * This model implements a \f$M\f$-phase flow of a fluid mixture
 * composed of \f$N\f$ chemical species. The phases are denoted by
 * lower index \f$\alpha \in \{ 1, \dots, M \}\f$. All fluid phases
 * are mixtures of \f$N \geq M - 1\f$ chemical species which are
 * denoted by the upper index \f$\kappa \in \{ 1, \dots, N \} \f$.
 *
 * The momentum approximation can be selected via "BaseFluxVariables":
 * Darcy (ImplicitDarcyFluxVariables) and Forchheimer (ImplicitForchheimerFluxVariables)
 * relations are available for all Box models.
 *
 * By inserting this into the equations for the conservation of the
 * mass of each component, one gets one mass-continuity equation for
 * each component \f$\kappa\f$
 * \f[
 \sum_{\kappa} \left(
    \phi \frac{\partial \left(\varrho_\alpha x_\alpha^\kappa S_\alpha\right)}{\partial t}
    +
    \mathrm{div}\;
    \left\{
        v_\alpha
       \frac{\varrho_\alpha}{\overline M_\alpha} x_\alpha^\kappa
    \right\}
    \right)
    = q^\kappa
    \f]
 * with \f$\overline M_\alpha\f$ being the average molar mass of the
 * phase \f$\alpha\f$: \f[ \overline M_\alpha = \sum_\kappa M^\kappa
 * \; x_\alpha^\kappa \f]
 *
 * For the missing \f$M\f$ model assumptions, the model assumes that
 * if a fluid phase is not present, the sum of the mole fractions of
 * this fluid phase is smaller than \f$1\f$, i.e.
 * \f[
 * \forall \alpha: S_\alpha = 0 \implies \sum_\kappa x_\alpha^\kappa \leq 1
 * \f]
 *
 * Also, if a fluid phase may be present at a given spatial location
 * its saturation must be positive:
 * \f[ \forall \alpha: \sum_\kappa x_\alpha^\kappa = 1 \implies S_\alpha \geq 0 \f]
 *
 * Since at any given spatial location, a phase is always either
 * present or not present, one of the strict equalities on the
 * right hand side is always true, i.e.
 * \f[ \forall \alpha: S_\alpha \left( \sum_\kappa x_\alpha^\kappa - 1 \right) = 0 \f]
 * always holds.
 *
 * These three equations constitute a non-linear complementarity
 * problem, which can be solved using so-called non-linear
 * complementarity functions \f$\Phi(a, b)\f$ which have the property
 * \f[\Phi(a,b) = 0 \iff a \geq0 \land b \geq0  \land a \cdot b = 0 \f]
 *
 * Several non-linear complementarity functions have been suggested,
 * e.g. the Fischer-Burmeister function
 * \f[ \Phi(a,b) = a + b - \sqrt{a^2 + b^2} \;. \f]
 * This model uses
 * \f[ \Phi(a,b) = \min \{a,  b \}\;, \f]
 * because of its piecewise linearity.
 *
 * These equations are then discretized using a fully-implicit vertex
 * centered finite volume scheme (often known as 'box'-scheme) for
 * spatial discretization and the implicit Euler method as temporal
 * discretization.
 *
 * The model assumes local thermodynamic equilibrium and uses the
 * following primary variables:
 * - The component fugacities \f$f^1, \dots, f^{N}\f$
 * - The pressure of the first phase \f$p_1\f$
 * - The saturations of the first \f$M-1\f$ phases \f$S_1, \dots, S_{M-1}\f$
 * - Temperature \f$T\f$ if the energy equation is enabled
 */
template<class TypeTag>
class MPNCModel : public GET_PROP_TYPE(TypeTag, BaseModel)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseModel) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef Dumux::MPNCVtkWriter<TypeTag> MPNCVtkWriter;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum {enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy)};
    enum {enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion)};
    enum {enableKinetic = GET_PROP_VALUE(TypeTag, EnableKinetic)};
    enum {numEnergyEquations = GET_PROP_VALUE(TypeTag, NumEnergyEquations)};
    enum {numPhases = GET_PROP_VALUE(TypeTag, NumPhases)};
    enum {numComponents = GET_PROP_VALUE(TypeTag, NumComponents)};
    enum {numEq = GET_PROP_VALUE(TypeTag, NumEq)};
    enum { numEnergyEqs = Indices::numPrimaryEnergyVars};
    enum { dimWorld = GridView::dimensionworld};
    enum { dim = GridView::dimension};

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;


public:
    MPNCModel()
    {
        enableSmoothUpwinding_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Implicit, EnableSmoothUpwinding);
        enablePartialReassemble_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Implicit, EnablePartialReassemble);
        enableJacobianRecycling_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Implicit, EnableJacobianRecycling);
        numDiffMethod_ = GET_PARAM_FROM_GROUP(TypeTag, int, Implicit, NumericDifferenceMethod);
    }

    void init(Problem &problem)
    {
        ParentType::init(problem);
        vtkWriter_ = std::make_shared<MPNCVtkWriter>(problem);

        if (this->gridView_().comm().rank() == 0)
            std::cout
                << "Initializing M-phase N-component model: \n"
                << "    phases: " << numPhases << "\n"
                << "    components: " << numComponents << "\n"
                << "    equations: " << numEq << "\n"
                << "    kinetic mass transfer: " << enableKinetic<< "\n"
                << "    number of energy equations: " << numEnergyEquations<< "\n"
                << "    diffusion: " << enableDiffusion << "\n"
                << "    energy equation: " << enableEnergy << "\n"
                << "    smooth upwinding: " << enableSmoothUpwinding_ << "\n"
                << "    partial jacobian reassembly: " << enablePartialReassemble_ << "\n"
                << "    numeric differentiation method: " << numDiffMethod_ << " (-1: backward, 0: central, +1 forward)\n"
                << "    jacobian recycling: " << enableJacobianRecycling_ << "\n";
    }

    /*!
     * \brief Compute the total storage inside one phase of all
     *        conservation quantities.
     *
     *        \param phaseStorage The conserved quantity within the phase in the whole domain
     *        \param phaseIdx The local index of the phases
     */
    void globalPhaseStorage(PrimaryVariables &phaseStorage, const unsigned int phaseIdx)
    {
        phaseStorage = 0;

        ElementIterator eIt = this->gridView_().template begin<0>();
        const ElementIterator eEndIt = this->gridView_().template end<0>();
        for (; eIt != eEndIt; ++eIt) {
            this->localResidual().addPhaseStorage(phaseStorage, *eIt, phaseIdx);
        }

        if (this->gridView_().comm().size() > 1)
            phaseStorage = this->gridView_().comm().sum(phaseStorage);
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

    std::shared_ptr<MPNCVtkWriter> vtkWriter_;

private:
    bool enableSmoothUpwinding_;
    bool enablePartialReassemble_;
    bool enableJacobianRecycling_;
    int numDiffMethod_;
};

}

#include "mpncpropertydefaults.hh"

#endif
