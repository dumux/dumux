// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test problems for the complex-valued Helmholtz model
 */
#ifndef DUMUX_TEST_COMPLEX_HELMHOLTZ_PROBLEM_HH
#define DUMUX_TEST_COMPLEX_HELMHOLTZ_PROBLEM_HH

#include <cmath>
#include <complex>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblem.hh>

namespace Dumux {

/*!
 * \brief Base problem holding the complex wave number
 */
template<class TypeTag>
class ComplexHelmholtzProblemBase : public FVProblem<TypeTag>
{
    using ParentType = FVProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using Complex = std::complex<Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using GlobalPosition = typename GridGeometry::LocalView::Element::Geometry::GlobalCoordinate;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;

    ComplexHelmholtzProblemBase(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , waveNumberSquared_(0.0)
    {}

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    void setWaveNumberSquared(const Complex& waveNumberSquared)
    { waveNumberSquared_ = waveNumberSquared; }

    const Complex& waveNumberSquared() const
    { return waveNumberSquared_; }

private:
    Complex waveNumberSquared_;
};

/*!
 * \brief Manufactured solution \f$ u = a \prod_i \sin(\pi x_i) \f$ with complex amplitude \f$ a \f$
 *        on the unit cube with homogeneous Dirichlet conditions,
 *        so that \f$ f = (d \pi^2 - k^2) u \f$.
 */
template<class TypeTag>
class ComplexHelmholtzManufacturedProblem : public ComplexHelmholtzProblemBase<TypeTag>
{
    using ParentType = ComplexHelmholtzProblemBase<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr int dimWorld = GridGeometry::GridView::dimensionworld;
public:
    using typename ParentType::Complex;
    using typename ParentType::PrimaryVariables;
    using typename ParentType::NumEqVector;
    using typename ParentType::GlobalPosition;

    ComplexHelmholtzManufacturedProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        amplitude_ = Complex(1.0, 0.5);
        this->setWaveNumberSquared(Complex(
            getParam<Scalar>("Problem.WaveNumberSquaredReal", 10.0),
            getParam<Scalar>("Problem.WaveNumberSquaredImag", -3.0)
        ));
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(exactSolution(globalPos)); }

    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector source(0.0);
        source[0] = (dimWorld*M_PI*M_PI - this->waveNumberSquared())*exactSolution(globalPos);
        return source;
    }

    Complex exactSolution(const GlobalPosition& globalPos) const
    {
        double product = 1.0;
        for (int i = 0; i < dimWorld; ++i)
            product *= std::sin(M_PI*globalPos[i]);
        return amplitude_*product;
    }

private:
    Complex amplitude_;
};

/*!
 * \brief Compact source in a damped cavity, used for the resonance sweep
 */
template<class TypeTag>
class ComplexHelmholtzSweepProblem : public ComplexHelmholtzProblemBase<TypeTag>
{
    using ParentType = ComplexHelmholtzProblemBase<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
public:
    using typename ParentType::PrimaryVariables;
    using typename ParentType::NumEqVector;
    using typename ParentType::GlobalPosition;
    using ParentType::ParentType;

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector source(0.0);
        const double R = 0.1;
        if (std::hypot(globalPos[0]-0.37, globalPos[1]-0.43) < R)
            source[0] = 1.0/(M_PI*R*R);
        return source;
    }
};

} // end namespace Dumux

#endif
