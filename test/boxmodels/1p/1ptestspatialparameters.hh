// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Bernd Flemisch                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
/*!
 * \file
 *
 * \brief The spatial parameters class for the test problem using the
 *        1p box model
 */
#ifndef DUMUX_1P_TEST_SPATIALPARAMETERS_HH
#define DUMUX_1P_TEST_SPATIALPARAMETERS_HH

#include <dumux/material/spatialparameters/boxspatialparameters1p.hh>

namespace Dumux
{

/*!
 * \ingroup OnePBoxModel
 * \ingroup BoxTestProblems
 *
 * \brief The spatial parameters class for the test problem using the
 *        1p box model
 */
template<class TypeTag>
class OnePTestSpatialParameters : public BoxSpatialParametersOneP<TypeTag>
{
    typedef BoxSpatialParametersOneP<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;
    typedef typename GridView::template Codim<0>::Entity Element;

public:
    OnePTestSpatialParameters(const GridView& gridView)
        : ParentType(gridView)
    {
        try
        {
            lensLowerLeft_[0] = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParameters.lensLowerLeftX);
            lensLowerLeft_[1] = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParameters.lensLowerLeftY);
            lensUpperRight_[0] = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParameters.lensUpperRightX);
            lensUpperRight_[1] = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParameters.lensUpperRightY);

            permeability_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParameters.permeability);
            permeabilityLens_=GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParameters.permeabilityLens);
        }
        catch (Dumux::ParameterException &e) {
            std::cerr << e << ". Abort!\n";
            exit(1) ;
        }
        catch (...) {
            std::cerr << "Unknown exception thrown!\n";
            exit(1);
        }
    }

    /*!
     * \brief Return the intrinsic permeability for the current sub-control volume.
     *
     * \param element The current finite element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index sub-control volume face where the
     *                      intrinsic velocity ought to be calculated.
     */
    Scalar intrinsicPermeability(const Element &element,
                                 const FVElementGeometry &fvElemGeom,
                                 int scvIdx) const
    {
        const GlobalPosition &globalPos = fvElemGeom.subContVol[scvIdx].global;

        if (isInLens_(globalPos))
            return permeabilityLens_;
        else
            return permeability_;
    }

    /*! \brief Define the porosity.
   *
   * \param element The finite element
   * \param fvElemGeom The finite volume geometry
   * \param scvIdx The local index of the sub-control volume where
   */
    Scalar porosity(const Element &element,
                    const FVElementGeometry &fvElemGeom,
                    int scvIdx) const
    { return 0.4; }

private:
    bool isInLens_(const GlobalPosition &pos) const
    {
        for (int i = 0; i < dim; ++i) {
            if (pos[i] < lensLowerLeft_[i] || pos[i] > lensUpperRight_[i])
                return false;
        }
        return true;
    }

    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;

    Scalar permeability_, permeabilityLens_;
};

} // end namespace
#endif

