// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FacetTests
 * \brief The problem for the bulk domain of the tracer facet coupling test.
 */

#ifndef DUMUX_TEST_TPFAFACETCOUPLING_TRACER_BULK_PROBLEM_HH
#define DUMUX_TEST_TPFAFACETCOUPLING_TRACER_BULK_PROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup FacetTests
 * \brief The problem for the bulk domain of the tracer facet coupling test.
 */
template <class TypeTag>
class TracerBulkProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    //! property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using typename ParentType::SpatialParams;

    TracerBulkProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                      std::shared_ptr<SpatialParams> spatialParams,
                      std::shared_ptr<CouplingManager> couplingManager,
                      const std::string& paramGroup = "")
    : ParentType(gridGeometry, spatialParams, paramGroup)
    , couplingManagerPtr_(couplingManager)
    , initialMassFraction_(getParamFromGroup<Scalar>(paramGroup, "Problem.ContaminationMassFraction"))
    {
        // stating in the console whether mole or mass fractions are used
        const auto problemName = getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
        std::cout<< "problem " << problemName << " uses " << (useMoles ? "mole" : "mass") << " fractions" << '\n';
        problemName_  =  getParamFromGroup<std::string>(this->paramGroup(), "Vtk.OutputName") + "_" + problemName;
    }

    //! The problem name.
    const std::string& name() const
    { return problemName_; }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    //! Specifies the type of interior boundary condition at a given position.
    BoundaryTypes interiorBoundaryTypes(const Element& element, const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        if (isInContaminatedRegion_(globalPos))
            return PrimaryVariables(initialMassFraction_);
        return PrimaryVariables(0.0);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param elemFluxVarsCache Flux variables caches for all faces in stencil
     * \param scvf The sub control volume face
     *
     * Negative values mean influx.
     * E.g. for the mass balance that would the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    template<class ElementVolumeVariables, class ElementFluxVarsCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVarsCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        // get the volume flux on this segment
        const auto flux = this->spatialParams().volumeFlux(element, fvGeometry, elemVolVars, scvf);
        if (flux > 0.0)
        {
            const auto& insideVolVars = elemVolVars[fvGeometry.scv(scvf.insideScvIdx())];
            const auto tracerFlux = insideVolVars.massFraction(/*phaseIdx*/0, /*compIdx*/0)*flux;
            return NumEqVector(tracerFlux);
        }

        return NumEqVector(0.0);
    }

    //! Returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

private:
    bool isInContaminatedRegion_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] > 0.75 && globalPos[0] < 1.3
               && globalPos[1] > 5.3 && globalPos[1] < 5.8;
    }

    std::shared_ptr<CouplingManager> couplingManagerPtr_;
    Scalar initialMassFraction_;
    std::string problemName_;
};

} // end namespace Dumux

#endif
