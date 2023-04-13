// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/**
 * \file
 * \ingroup TracerTests
 * \brief Definition of a problem for the tracer problem:
 * Two tracers are diluted by diffusion and constant two-phase saturation field.
 */

#ifndef DUMUX_TRACER_MULTIPHASE_TEST_PROBLEM_HH
#define DUMUX_TRACER_MULTIPHASE_TEST_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup TracerTests
 * \brief A problem, where two tracers are diluted by diffusion and constant two-phase saturation field.
 */
template <class TypeTag>
class TracerTest : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;

    //! property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    TracerTest(std::shared_ptr<const GridGeometry> gg)
    : ParentType(gg)
    {
        // stating in the console whether mole or mass fractions are used
        if(useMoles)
            std::cout<<"problem uses mole fractions" << '\n';
        else
            std::cout<<"problem uses mass fractions" << '\n';
    }

    /*!
     * \brief The boundary types
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();

        // tracer is only present in water
        if (this->spatialParams().isWater(globalPos) && globalPos[0] < eps_)
            values.setAllDirichlet();

        return values;
    }

    /*!
     * \brief The Dirichlet values
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        if (globalPos[0] < eps_)
            return PrimaryVariables({1e-8, 1e-8});
        else
            return PrimaryVariables({0.0, 0.0});
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * This is the method for the case where the Neumann condition is
     * potentially solution dependent
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
    template<class ElemFluxVarsCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElemFluxVarsCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);
        const auto& globalPos = scvf.ipGlobal();
        if (globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_)
        {
            const auto& volVars = elemVolVars[scvf.insideScvIdx()];
            values = NumEqVector({volVars.massFraction(0,0), volVars.massFraction(0,1)});
            values *= volVars.density()*(this->spatialParams().velocity(scvf)*scvf.unitOuterNormal());
        }
        return values;
    }

    /*!
     * \brief The initial values
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        return PrimaryVariables({0.0, 0.0});
    }
private:
    static constexpr Scalar eps_ = 1e-7;
};

} // end namespace Dumux

#endif
