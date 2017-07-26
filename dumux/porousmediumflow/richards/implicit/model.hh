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
* \brief Adaption of the fully implicit scheme to the Richards model.
*/
#ifndef DUMUX_RICHARDS_MODEL_HH
#define DUMUX_RICHARDS_MODEL_HH

#include <dumux/implicit/model.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/model.hh>
#include "properties.hh"
#include "primaryvariableswitch.hh"

namespace Dumux
{
/*!
 * \ingroup RichardsModel
 *
 * \brief This model which implements a variant of the Richards'
 *        equation for quasi-twophase flow.
 *
 * In the unsaturated zone, Richards' equation
 \f[
 \frac{\partial\;\phi S_w \varrho_w}{\partial t}
 -
 \text{div} \left\lbrace
 \varrho_w \frac{k_{rw}}{\mu_w} \; \mathbf{K} \;
 \left( \text{\textbf{grad}}
 p_w - \varrho_w \textbf{g}
 \right)
 \right\rbrace
 =
 q_w,
 \f]
 * is frequently used to
 * approximate the water distribution above the groundwater level.
 *
 * It can be derived from the two-phase equations, i.e.
 \f[
 \phi\frac{\partial S_\alpha \varrho_\alpha}{\partial t}
 -
 \text{div} \left\lbrace
 \varrho_\alpha \frac{k_{r\alpha}}{\mu_\alpha}\; \mathbf{K} \;
 \left( \text{\textbf{grad}}
 p_\alpha - \varrho_\alpha \textbf{g}
 \right)
 \right\rbrace
 =
 q_\alpha,
 \f]
 * where \f$\alpha \in \{w, n\}\f$ is the fluid phase,
 * \f$\kappa \in \{ w, a \}\f$ are the components,
 * \f$\rho_\alpha\f$ is the fluid density, \f$S_\alpha\f$ is the fluid
 * saturation, \f$\phi\f$ is the porosity of the soil,
 * \f$k_{r\alpha}\f$ is the relative permeability for the fluid,
 * \f$\mu_\alpha\f$ is the fluid's dynamic viscosity, \f$\mathbf{K}\f$ is the
 * intrinsic permeability, \f$p_\alpha\f$ is the fluid pressure and
 * \f$g\f$ is the potential of the gravity field.
 *
 * In contrast to the full two-phase model, the Richards model assumes
 * gas as the non-wetting fluid and that it exhibits a much lower
 * viscosity than the (liquid) wetting phase. (For example at
 * atmospheric pressure and at room temperature, the viscosity of air
 * is only about \f$1\%\f$ of the viscosity of liquid water.) As a
 * consequence, the \f$\frac{k_{r\alpha}}{\mu_\alpha}\f$ term
 * typically is much larger for the gas phase than for the wetting
 * phase. For this reason, the Richards model assumes that
 * \f$\frac{k_{rn}}{\mu_n}\f$ is infinitly large. This implies that
 * the pressure of the gas phase is equivalent to the static pressure
 * distribution and that therefore, mass conservation only needs to be
 * considered for the wetting phase.
 *
 * The model thus choses the absolute pressure of the wetting phase
 * \f$p_w\f$ as its only primary variable. The wetting phase
 * saturation is calculated using the inverse of the capillary
 * pressure, i.e.
 \f[
 S_w = p_c^{-1}(p_n - p_w)
 \f]
 * holds, where \f$p_n\f$ is a given reference pressure. Nota bene,
 * that the last step is assumes that the capillary
 * pressure-saturation curve can be uniquely inverted, so it is not
 * possible to set the capillary pressure to zero when using the
 * Richards model!
 */
template<class TypeTag >
class RichardsModel : public GET_PROP_TYPE(TypeTag, BaseModel)
{
    // the parent class needs to access the variable switch
    friend typename GET_PROP_TYPE(TypeTag, BaseModel);
    using ParentType = typename GET_PROP_TYPE(TypeTag, BaseModel);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    using NonIsothermalModel = Dumux::NonIsothermalModel<TypeTag>;

    enum {
        nPhaseIdx = Indices::nPhaseIdx,
        wPhaseIdx = Indices::wPhaseIdx
    };

    static constexpr int dim = GridView::dimension;
    static constexpr bool isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox);
    enum { dofCodim = isBox ? dim : 0 };

    static constexpr bool enableWaterDiffusionInAir
        = GET_PROP_VALUE(TypeTag, EnableWaterDiffusionInAir);

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
        vtkOutputModule.addSecondaryVariable("Sw", [](const VolumeVariables& v){ return v.saturation(wPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("Sn", [](const VolumeVariables& v){ return v.saturation(nPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("pw", [](const VolumeVariables& v){ return v.pressure(wPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("pn", [](const VolumeVariables& v){ return v.pressure(nPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("pc", [](const VolumeVariables& v){ return v.capillaryPressure(); });
        vtkOutputModule.addSecondaryVariable("density", [](const VolumeVariables& v){ return v.density(wPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("mobility", [](const VolumeVariables& v){ return v.mobility(wPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("kr", [](const VolumeVariables& v){ return v.relativePermeability(wPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("porosity", [](const VolumeVariables& v){ return v.porosity(); });
        if(GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
            vtkOutputModule.addSecondaryVariable("pressure head", [](const VolumeVariables& v){ return v.pressureHead(wPhaseIdx); });
        if (enableWaterDiffusionInAir)
            vtkOutputModule.addSecondaryVariable("x^w_air", [](const VolumeVariables& v){ return v.moleFraction(1, 0); });
        vtkOutputModule.addSecondaryVariable("water content", [](const VolumeVariables& v){ return v.waterContent(wPhaseIdx); });

        NonIsothermalModel::maybeAddTemperature(vtkOutputModule);
    }

    /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    template<class VtkOutputModule>
    void addVtkOutputFields(VtkOutputModule& outputModule) const
    {
        auto& phasePresence = outputModule.createScalarField("phase presence", dofCodim);
        for (std::size_t i = 0; i < phasePresence.size(); ++i)
            phasePresence[i] = this->curSol()[i].state();
    }

    /*!
     * \brief One Newton iteration was finished.
     * \param uCurrent The solution after the current Newton iteration
     */
    template<typename T = TypeTag>
    typename std::enable_if<GET_PROP_VALUE(T, EnableGlobalVolumeVariablesCache), void>::type
    newtonEndStep()
    {
        // \todo resize volvars vector if grid was adapted

        // update the variable switch
        switchFlag_ = priVarSwitch_().update(this->problem_(), this->curSol());

        // update the secondary variables if global caching is enabled
        // \note we only updated if phase presence changed as the volume variables
        //       are already updated once by the switch
        for (const auto& element : elements(this->problem_().gridView()))
        {
            // make sure FVElementGeometry & vol vars are bound to the element
            auto fvGeometry = localView(this->globalFvGeometry());
            fvGeometry.bindElement(element);

            if (switchFlag_)
            {
                for (auto&& scv : scvs(fvGeometry))
                {
                    auto dofIdxGlobal = scv.dofIndex();
                    if (priVarSwitch_().wasSwitched(dofIdxGlobal))
                    {
                        const auto eIdx = this->problem_().elementMapper().index(element);
                        const auto elemSol = this->elementSolution(element, this->curSol());
                        this->nonConstCurGlobalVolVars().volVars(eIdx, scv.indexInElement()).update(elemSol,
                                                                                                    this->problem_(),
                                                                                                    element,
                                                                                                    scv);
                    }
                }
            }

            // handle the boundary volume variables for cell-centered models
            if(!isBox)
            {
                for (auto&& scvf : scvfs(fvGeometry))
                {
                    // if we are not on a boundary, skip the rest
                    if (!scvf.boundary())
                        continue;

                    // check if boundary is a pure dirichlet boundary
                    const auto bcTypes = this->problem_().boundaryTypes(element, scvf);
                    if (bcTypes.hasOnlyDirichlet())
                    {
                        const auto insideScvIdx = scvf.insideScvIdx();
                        const auto& insideScv = fvGeometry.scv(insideScvIdx);
                        const auto elemSol = ElementSolutionVector{this->problem_().dirichlet(element, scvf)};

                        this->nonConstCurGlobalVolVars().volVars(scvf.outsideScvIdx(), 0/*indexInElement*/).update(elemSol, this->problem_(), element, insideScv);
                    }
                }
            }
        }
    }

    /*!
     * \brief One Newton iteration was finished.
     * \param uCurrent The solution after the current Newton iteration
     */
    template<typename T = TypeTag>
    typename std::enable_if<!GET_PROP_VALUE(T, EnableGlobalVolumeVariablesCache), void>::type
    newtonEndStep()
    {
        // update the variable switch
        switchFlag_ = priVarSwitch_().update(this->problem_(), this->curSol());
    }

    /*!
     * \brief Called by the update() method if applying the Newton
     *        method was unsuccessful.
     */
    void updateFailed()
    {
        ParentType::updateFailed();
        // reset privar switch flag
        switchFlag_ = false;
    }

    /*!
     * \brief Called by the problem if a time integration was
     *        successful, post processing of the solution is done and the
     *        result has been written to disk.
     *
     * This should prepare the model for the next time integration.
     */
    void advanceTimeLevel()
    {
        ParentType::advanceTimeLevel();
        // reset privar switch flag
        switchFlag_ = false;
    }

    /*!
     * \brief Returns true if the primary variables were switched for
     *        at least one dof after the last timestep.
     */
    bool switched() const
    {
        return switchFlag_;
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

        outStream << this->curSol()[dofIdxGlobal].state() << " ";
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

        int phasePresence;
        inStream >> phasePresence;

        this->curSol()[dofIdxGlobal].setState(phasePresence);
        this->prevSol()[dofIdxGlobal].setState(phasePresence);
    }

    const Dumux::ExtendedRichardsPrimaryVariableSwitch<TypeTag>& priVarSwitch() const
    { return switch_; }

protected:

    Dumux::ExtendedRichardsPrimaryVariableSwitch<TypeTag>& priVarSwitch_()
    { return switch_; }

    /*!
     * \brief Applies the initial solution for all vertices of the grid.
     *
     * \todo the initial condition needs to be unique for
     *       each vertex. we should think about the API...
     */
    void applyInitialSolution_()
    {
        ParentType::applyInitialSolution_();

        // initialize the primary variable switch
        priVarSwitch_().init(this->problem_());
    }

    //! the class handling the primary variable switch
    ExtendedRichardsPrimaryVariableSwitch<TypeTag> switch_;
    bool switchFlag_;
};

} // end namespace Dumux

#include "propertydefaults.hh"

#endif
