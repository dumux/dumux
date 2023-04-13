// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePNCTests
 * \brief Definition of the spatial parameters for the 1pnc problems.
 */

#ifndef DUMUX_1PNC_TEST_SPATIAL_PARAMS_HH
#define DUMUX_1PNC_TEST_SPATIAL_PARAMS_HH

#include <dune/common/exceptions.hh>
#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup OnePNCTests
 * \brief Definition of the spatial parameters for the 1pnc test problems.
 */
template<class GridGeometry, class Scalar>
class OnePNCTestSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                             OnePNCTestSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                                           OnePNCTestSpatialParams<GridGeometry, Scalar>>;

    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Dune::FieldVector<Scalar, dimWorld>;

public:
    // export permeability type
    using PermeabilityType = Scalar;
    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    OnePNCTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        permeability_ = 1e-10;
        porosity_ = 0.4;
        alphaL_ = getParam<Scalar>("Problem.AlphaL");
        alphaT_ = getParam<Scalar>("Problem.AlphaT");
        dispersionTensorCoefficients_ = getParam<std::vector<Scalar>>("Problem.DispersionTensor");

        if (dispersionTensorCoefficients_.size() > 1)
        {
            if (dispersionTensorCoefficients_.size() != (dimWorld*dimWorld))
                DUNE_THROW(Dune::InvalidStateException, "For anisotropic dispersion tensors, please list all entries (dim x dim).");

            int k = 0;
            for (int i = 0; i < dimWorld; i++)
            {
                for (int j = 0; j < dimWorld; j++)
                {
                    dispersionTensor_[i][j] = dispersionTensorCoefficients_[k];
                    k++;
                }
            }
        }
        else
        {
            if (dispersionTensorCoefficients_.size() != 1)
                DUNE_THROW(Dune::InvalidStateException, "For isotropic dispersion tensors, please one scalar value.");

            for (int i = 0; i < dimWorld; i++)
                dispersionTensor_[i][i] = dispersionTensorCoefficients_[0];
        }
    }

    /*!
     * \brief Defines the intrinsic permeability \f$\mathrm{[m^2]}\f$.
     *
     * \param globalPos The global position
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return permeability_; }

    /*!
     * \brief Defines the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }

    /*!
     * \brief Returns the temperature within the domain [K].
     *
     * \param globalPos The global position
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    { return 273.15 + 20; } // in [K]

    /*!
     * \brief Defines the dispersion tensor \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    std::array<Scalar, 2> dispersionAlphas(const GlobalPosition& globalPos, int phaseIdx = 0, int compIdx = 0) const
    { return { alphaL_, alphaT_ }; }

    /*!
     * \brief Defines the dispersion tensor \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    const DimWorldMatrix &dispersionTensor(const GlobalPosition& globalPos, int phaseIdx = 0, int compIdx = 0) const
    { return dispersionTensor_; }

private:
    Scalar permeability_;
    Scalar porosity_;
    Scalar alphaL_;
    Scalar alphaT_;
    std::vector<Scalar> dispersionTensorCoefficients_;
    DimWorldMatrix dispersionTensor_;
};

} // end namespace Dumux

#endif
