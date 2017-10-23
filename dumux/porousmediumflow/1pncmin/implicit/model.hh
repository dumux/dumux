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

#include <dumux/porousmediumflow/1pnc/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/velocityoutput.hh>

namespace Dumux
{
/*!
 * \ingroup OnePNCMinModel
 * \brief Adaption of the fully implicit scheme to the
 *        one-phase n-component fully implicit model with additional solid/mineral phases.
 *
 * This model implements one-phase n-component flow of a compressible fluid composed of
 * the n components \f$\kappa \f$ in combination with mineral precipitation and dissolution
 * of the solid phases. The standard multiphase Darcy
 * approach is used as the equation for the conservation of momentum:
 * \f[
 v = - \frac{k_{r}}{\mu_} \mbox{\bf K}
 \left(\text{grad}\, p - \varrho_{f} \mbox{\bf g} \right)
 * \f]
 *
 * By inserting this into the equations for the conservation of the
 * components, one gets one transport equation for each component
 * \f{eqnarray}
 && \frac{\partial ( \varrho_f X^\kappa \phi  )}
 {\partial t} -  \text{div} \left\{ \varrho_f X^\kappa
 \frac{k_{r}}{\mu} \mbox{\bf K}
 (\text{grad}\, p - \varrho_{f}  \mbox{\bf g}) \right\}
 \nonumber \\ \nonumber \\
 &-& \text{div} \left\{{\bf D_{pm}^\kappa} \varrho_{f} \text{grad}\, X^\kappa \right\}
 -  q_\kappa = 0 \qquad \kappa \in \{w, a,\cdots \} \,
 \f}
 *
 * The solid or mineral phases are assumed to consist of a single component.
 * Their mass balance consist only of a storage and a source term:
 *  \f$\frac{\partial \varrho_\lambda \phi_\lambda )} {\partial t}
 *  = q_\lambda\f$
 *
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial and the implicit Euler method as time
 * discretization.
 *
 * The primary variables are the pressure \f$p\f$ and the mole fractions of the
 * dissolved components \f$x^k\f$. The primary variable of the solid phases is the volume
 * fraction \f$\phi_\lambda = \frac{V_\lambda}{V_{total}}\f$.
 *
 * The source an sink terms link the mass balances of the n-transported component to the
 * solid phases. The porosity \f$\phi\f$ is updated according to the reduction of the initial * (or solid-phase-free porous medium) porosity \f$\phi_0\f$ by the accumulated volume
 * fractions of the solid phases:
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

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

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

        vtkOutputModule.addSecondaryVariable("Kxx",
                                             [this](const VolumeVariables& v){ return this->perm_(v.permeability())[0][0]; });
        if (dim >= 2)
            vtkOutputModule.addSecondaryVariable("Kyy",
                                                 [this](const VolumeVariables& v){ return this->perm_(v.permeability())[1][1]; });
        if (dim >= 3)
            vtkOutputModule.addSecondaryVariable("Kzz",
                                                 [this](const VolumeVariables& v){ return this->perm_(v.permeability())[2][2]; });

        vtkOutputModule.addSecondaryVariable("permeabilityFactor", [](const VolumeVariables& v)
                                             { return v.permeabilityFactor(); });


        for (int sPhaseIdx = 0; sPhaseIdx < numSPhases; ++sPhaseIdx)
            vtkOutputModule.addSecondaryVariable("precipitateVolumeFraction_" + FluidSystem::phaseName(numPhases + sPhaseIdx),
                                                 [sPhaseIdx](const VolumeVariables& v)
                                                 { return v.precipitateVolumeFraction(numPhases + sPhaseIdx); });
    }


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

private:

    Tensor perm_(Scalar perm) const
    {
        Tensor K(0.0);

        for(int i=0; i<dimWorld; i++)
            K[i][i] = perm;

       return K;
    }

    const Tensor& perm_(const Tensor& perm) const
    {
       return perm;
    }
};

} // end namespace Dumux

#include "propertydefaults.hh"

#endif
