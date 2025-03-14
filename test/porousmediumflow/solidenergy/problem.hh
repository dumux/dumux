// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/**
 * \file
 * \ingroup SolidEnergyTests
 * \brief A piece of stone heated by a network of channels
 */
#ifndef DUMUX_TEST_SOLIDENERGY_PROBLEM_HH
#define DUMUX_TEST_SOLIDENERGY_PROBLEM_HH

#include <string>
#include <cmath>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup SolidEnergyTests
 * \brief A piece of stone heated by a network of channels
 */
template <class TypeTag>
class SolidEnergyProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NeumannValues = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;

public:
    SolidEnergyProblem(std::shared_ptr<const GridGeometry> gridGeometry, const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    {
        name_ = getParam<std::string>("Problem.Name");
        temperatureHigh_ = getParam<Scalar>("Problem.TemperatureHigh");
        temperatureInit_ = getParam<Scalar>("Problem.TemperatureInit");
        lambdaSolid_ = getParam<Scalar>("Component.SolidThermalConductivity");
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    {
        return name_;
    }
    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        bcTypes.setAllNeumann();
        return bcTypes;
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
    NeumannValues neumann(const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const ElementVolumeVariables& elemVolVars,
                          const ElementFluxVariablesCache& elemFluxVarsCache,
                          const SubControlVolumeFace& scvf) const
    {
        const auto& volVars = elemVolVars[scvf.insideScvIdx()];
        const Scalar value = 0.02*(volVars.temperatureSolid() - temperatureInit_)/0.002;
        return NeumannValues(value);
    }



    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        return PrimaryVariables(temperatureInit_);
    }

    /*!
     * \brief Applies a vector of point sources. The point sources
     *        are possibly solution dependent.
     *
     * \param pointSources A vector of PointSource s that contain
              source values for all phases and space positions.
     *
     * For this method, the values method of the point source
     * has to return the absolute rate values in units
     * \f$ [ \textnormal{unit of conserved quantity} / s ] \f$.
     * Positive values mean that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / s ] \f$.
     */
    void addPointSources(std::vector<PointSource>& pointSources) const
    {
        // add point sources along a sinus
        const auto& gg = this->gridGeometry();
        const auto ampl = (gg.bBoxMax()[1] - gg.bBoxMin()[1])/3.0;
        const auto len = (gg.bBoxMax()[0] - gg.bBoxMin()[0]);
        const auto freq = 1.0/len;

        for (int i = 0; i < 100; ++i)
        {
            const auto posx = double(i)/99.0*len;
            const auto posy = ampl*std::sin(2*M_PI*freq*posx) + ampl*3.0/2.0;
            pointSources.emplace_back(GlobalPosition({posx, posy}), typename PointSource::Values({len/100.0}));
        }
    }

    /*!
     * \brief Evaluate the point sources (added by addPointSources)
     *        for all phases within a given sub-control-volume.
     *
     * This is the method for the case where the point source is
     * solution dependent
     *
     * \param source A single point source
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The sub control volume
     *
     * For this method, the values() method of the point sources returns
     * the absolute conserved quantity rate generated or annihilate in
     * units \f$ [ \textnormal{unit of conserved quantity} / s ] \f$.
     * Positive values mean that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / s ] \f$.
     */
    void pointSource(PointSource& source,
                     const Element &element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolume &scv) const
    {
        const auto& volVars = elemVolVars[scv];
        const Scalar charLen = 0.1;
        static const Scalar nusselt = getParam<Scalar>("Problem.NusseltNumber");
        const Scalar interfacialArea = 2*M_PI*0.001;
        const Scalar lambdaSolidFluidTransfer = lambdaSolid_;
        source *= -nusselt*interfacialArea*lambdaSolidFluidTransfer*(volVars.temperatureSolid() - temperatureHigh_)/charLen;
    }

    // \}

private:
    Scalar temperatureHigh_, temperatureInit_, lambdaSolid_;
    static constexpr Scalar eps_ = 1e-6;
    std::string name_;
};

} //end namespace Dumux

#endif
