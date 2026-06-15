// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BedloadTests
 * \copydoc Dumux::BumpTestProblemBedload
 */
#ifndef DUMUX_BUMP_TEST_PROBLEM_BEDLOAD_HH
#define DUMUX_BUMP_TEST_PROBLEM_BEDLOAD_HH

#include "spatialparams_bedload.hh"
#include <dumux/common/parameters.hh>

// needed to calculate the analytical solution
#include <vector>
#include <cmath>
#include <functional>
#include <limits>

#include <dumux/freeflow/bedload/problem.hh>

namespace Dumux {

template <class TypeTag>
class BumpTestProblemBedload;

/*!
 * \ingroup BedloadTests
 * \brief The problem class for the bump test (bedload part)
 *
 * This test case involves a weak interaction between flow and mobile bed
 * and has an approximate analytical solution. The test case and the analytical
 * solution can be found in "Numerical Techniqes for Morphodynamic Modelling"
 * Hudson 2001, page 82 ff.
 */
template <class TypeTag>
class BumpTestProblemBedload : public BedloadProblem<TypeTag>
{
    using ParentType = BedloadProblem<TypeTag>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using NeumannFluxes = NumEqVector;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using VolumeVariables = typename ElementVolumeVariables::VolumeVariables;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

public:
    BumpTestProblemBedload(std::shared_ptr<const GridGeometry> gridGeometry,
                           std::shared_ptr<GetPropType<TypeTag, Properties::SpatialParams>> spatialParams)
    : ParentType(gridGeometry, spatialParams)
    {
        using std::sin;
        using std::pow;

        name_ = getParam<std::string>("Problem.Name");
        nGrainClasses_ = Dumux::getParam<int>("Sediment.NumberGrainClasses");
        if (nGrainClasses_ != TypeTag::nGrainClasses) {
            DUNE_THROW(Dune::InvalidStateException, "Please set 'Sediment.NumberGrainClasses = 1'.");
        }
        flowRate_ = Dumux::getParam<Scalar>("ShallowWater.FlowRate");
        freeSurfaceOutflow_ = Dumux::getParam<Scalar>("ShallowWater.FreeSurfaceOutflow");

        // prepare the data needed to initialise the layer model
        xBoundaries_ = getParam<std::array<Scalar, 2>>("Grid.Positions0");
        const Scalar initialBumpHight = peakInitialBump_ - initialBedLevel_;
        int nElems = gridGeometry->gridView().size(0);
        initialBedSurface_.resize(nElems);
        for (const auto& element : elements(gridGeometry->gridView()))
        {
            // Calculate initial bed surface
            auto elemId = this->gridGeometry().elementMapper().index(element);
            if (element.geometry().center()[0] < 300.0) {
                initialBedSurface_[elemId] = initialBedLevel_;
            }
            else if (element.geometry().center()[0] < 500.0) {
                initialBedSurface_[elemId] = initialBedLevel_ + (initialBumpHight-initialBedLevel_)*pow(sin(M_PI*(element.geometry().center()[0]-300)/ 200), 2);
            }
            else {
                initialBedSurface_[elemId] = initialBedLevel_;
            }
        }
    }

    //! Get the analytical solution
    const std::vector<Scalar>& getAnalyticalSolution()
    {
        return analyticalBedSurface_;
    }

    //! Get the initial bed surface
    const std::vector<Scalar>& getInitialBedSurface()
    {
        return initialBedSurface_;
    }

    /*!
     * \brief Update the analytical solution
     *
     * More details how to calculate the analytical solutions can be found in "Numerical
     * Techniqes for Morphodynamic Modelling" Hudson 2001, page 82 ff.
     *
     * \param time the time at which the analytical solution shall be calculated
     */
    void updateAnalyticalSolution(const Scalar time)
    {
        using std::pow;
        using std::sin;

        // collect x value of all elements
        std::vector<Scalar> xCenters;
        xCenters.reserve(this->gridGeometry().gridView().size(0));
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            const auto center = element.geometry().center();
            xCenters.push_back(center[0]);
        }

        int grassExponent = 3;  // grassExponent = 3, because dumux uses the quadratic version of the Grass bedlaod discharge formula
        Scalar temp = this->spatialParams().grassAlpha() * grassExponent * pow(flowRate_, grassExponent) / (1.0 - this->spatialParams().porosity());

        auto f = [&](Scalar x0, Scalar x, Scalar t) {
            Scalar val = x0 + temp * t * pow(freeSurfaceOutflow_ - pow(sin(M_PI * (x0 - 300.0) / 200.0), 2), -grassExponent - 1) - x;
            return val * val;
        };

        auto xnull = [&](Scalar x, Scalar t) {
            if (x < 300.0 + temp * t * pow(freeSurfaceOutflow_, -grassExponent - 1) || x > 500.0 + temp * t * pow(freeSurfaceOutflow_, -grassExponent - 1)) {
                return x - temp * t * pow(freeSurfaceOutflow_, -grassExponent - 1);
            } else {
                auto g = [&](Scalar x0) { return f(x0, x, t); };
                return minimize_scalar(g, 300.0, 500.0);
            }
        };

        auto z = [&](Scalar x, Scalar t) {
            Scalar x0 = xnull(x, t);
            if (x0 < 300.0 || x0 > 500.0) {
                return 0.0;
            } else {
                return pow(sin(M_PI * (x0 - 300.0) / 200.0), 2);
            }
        };

        std::vector<Scalar> sol(xCenters.size());
        for (size_t j = 0; j < xCenters.size(); j++) {
            sol[j] = z(xCenters[j], time);
        }
        analyticalBedSurface_ = sol;
    }
    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    {
        return name_;
    }

     /*!
     * \brief Evaluate the source term for all balance equations within a given
     *        sub-control-volume.
     *
     * There are no sources in this test. In general, a source could be artificial
     * sediment feeding or suspended load.
     *
     * For this method, the \a values parameter stores the conserved quantity rate
     * generated or annihilate per volume unit. Positive values mean
     * that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the sediment mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The sub control volume
     */
     NumEqVector source(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const SubControlVolume &scv) const
    {
        NumEqVector values (0.0);

        return values;
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
     * \param globalPos The position for which the boundary type is set
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        bcTypes.setAllNeumann();
        return bcTypes;
    }

    /*!
     * \brief Specifies the Neumann boundary
     *
     * Specify the Neumann boundary for the bedload flux.
     *
     * For the bedload flux we implement a euqilibrium boundary condition. This means that as much sediment enters
     * or leaves the domain as is tranported in the corresponding boundary cell.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scvf The sub control volume face
     */
    NeumannFluxes neumann(const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const ElementVolumeVariables& elemVolVars,
                          const ElementFluxVariablesCache& elemFluxVarsCache,
                          const SubControlVolumeFace& scvf) const
    {
        NeumannFluxes values(0.0);

        // inflow bedload discharge
        if (scvf.center()[0] < xBoundaries_[0] + eps_)
        {
            values[0] = -bedloadBoundaryDischarge_;
        }
        // outflow bedload discharge
        else if (xBoundaries_[1] - eps_ < scvf.center()[0])
        {
            values[0] = bedloadBoundaryDischarge_;
        }
        return values;
    }

    /*!
     * \brief Get coupled variables
     *
     * Get the bed surface for exchange with the shallow water model.
     *
     * \param curSol The current solution
     * \param gridVariables The grid variables
     */
    template<class SolutionVector, class GridVariables>
    std::map<std::string, std::vector<Scalar>> getCoupledVariables(const SolutionVector& curSol,
                                                                   const GridVariables& gridVariables)
    {
        int size = this->gridGeometry().gridView().size(0);
        // make map with a string and a vector to store the variable names and their values
        std::map<std::string, std::vector<Scalar>> coupledVariables;
        // add the variable names
        coupledVariables["bedSurface"] = std::vector<Scalar>(size, 0.0);
        // add the values
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVariables.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, curSol);

            // get local index of the element
            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            // loop over the gost entities, too
            for (auto&& scv : scvs(fvGeometry))
            {
                coupledVariables["bedSurface"][eIdx] = elemVolVars[scv].bedSurface();
            }
        }
        return coupledVariables;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial values for a control volume.
     *
     * For this method, the \a values parameter stores primary variables.
     *
     * \param element The finite element
     */
    PrimaryVariables initial(const Element& element) const
    {
        PrimaryVariables values(0.0);
        auto elemId = this->gridGeometry().elementMapper().index(element);
        values[0] = (initialBedSurface_[elemId] - this->spatialParams().fixedGroundLevel(element)) * element.geometry().volume()
                    * (1 - this->spatialParams().porosity()) * this->spatialParams().grainDensity(0);
        return values;
    }

    // \}

    /*!
     * \brief Return the bottom of the active layer of an element.
     *
     * \param element The finite element
     */
    const Scalar bottomActiveLayer(const Element& element,
                                   const SubControlVolume& scv) const
    {
        // The bottom of the active layer is the fixed ground level because no layer model is used
        return this->spatialParams().fixedGroundLevel(element);
    }

private:
    // simple Brent search for minimum of a function in [a,b]
    Scalar minimize_scalar(const std::function<Scalar(Scalar)>& f, Scalar a, Scalar b, Scalar tol = 1e-6, int max_iter = 100)
    {
        const Scalar phi = (1 + std::sqrt(5.0)) / 2.0;
        Scalar x1 = b - (b - a) / phi;
        Scalar x2 = a + (b - a) / phi;
        Scalar f1 = f(x1);
        Scalar f2 = f(x2);

        int iter = 0;
        while ((b - a) > tol && iter < max_iter) {
            if (f1 < f2) {
                b = x2;
                x2 = x1;
                f2 = f1;
                x1 = b - (b - a) / phi;
                f1 = f(x1);
            } else {
                a = x1;
                x1 = x2;
                f1 = f2;
                x2 = a + (b - a) / phi;
                f2 = f(x2);
            }
            iter++;
        }
        return (f1 < f2) ? x1 : x2;
    }

    std::vector<Scalar> analyticalBedSurface_; // [m]
    Scalar bedloadBoundaryDischarge_ = 2.65;   // [kg/s]
    static constexpr Scalar eps_ = 1.0e-6;     // [-]
    Scalar flowRate_;                          // [m^2/s]
    Scalar freeSurfaceOutflow_;                // [m]
    Scalar initialBedLevel_ = 0.0;             // [m]
    std::vector<Scalar> initialBedSurface_;    // [m]
    std::string name_;
    int nGrainClasses_;                        // [-]
    Scalar peakInitialBump_ = 1.0;             // [m]
    std::array<Scalar, 2> xBoundaries_;
};

} //end namespace Dumux

#endif
