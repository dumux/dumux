// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
 * \brief Definition of the spatial parameters for the 1p2c
 *        outlfow problem.
 */
#ifndef DUMUX_1P2C_OUTFLOW_SPATIAL_PARAMETERS_HH
#define DUMUX_1P2C_OUTFLOW_SPATIAL_PARAMETERS_HH

#include <dumux/material/spatialparameters/boxspatialparameters1p.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{

/*!
 * \ingroup OnePTwoCBoxModel
 * \ingroup BoxTestProblems
 *
 * \brief Definition of the spatial parameters for the 1p2c
 *        outflow problem.
 */
template<class TypeTag>
class OnePTwoCOutflowSpatialParameters : public BoxSpatialParametersOneP<TypeTag>
{
    typedef BoxSpatialParametersOneP<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;
    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;


    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;


    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GridView::template Codim<0>::Entity Element;

    //typedef LinearMaterial<Scalar> EffMaterialLaw;
public:
    OnePTwoCOutflowSpatialParameters(const GridView &gv)
        : ParentType(gv)
    {
        permeability_ = 1e-10;
        porosity_ = 0.4;
        tortuosity_ = 0.28;
    }

    ~OnePTwoCOutflowSpatialParameters()
    {}


    /*!
     * \brief Update the spatial parameters with the flow solution
     *        after a timestep.
     *
     * \param globalSolution the global solution vector
     */
    void update(const SolutionVector &globalSolution)
    {
    };

    /*!
     * \brief Define the intrinsic permeability \f$\mathrm{[m^2]}\f$.
     *
     * \param element The current finite element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    const Scalar intrinsicPermeability(const Element &element,
                                       const FVElementGeometry &fvElemGeom,
                                       int scvIdx) const
    {
            return permeability_;
    }

    /*!
     * \brief Define the porosity \f$\mathrm{[-]}\f$.
     *
     * \param element The finite element
     * \param fvElemGeom The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     */
    double porosity(const Element &element,
                    const FVElementGeometry &fvElemGeom,
                    int scvIdx) const
    {
            return porosity_;
    }

    /*!
     * \brief Define the tortuosity \f$\mathrm{[-]}\f$.
     *
     * \param element The finite element
     * \param fvElemGeom The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     */
    double tortuosity(const Element &element,
                    const FVElementGeometry &fvElemGeom,
                    int scvIdx) const
    {
            return tortuosity_;
    }

    /*!
     * \brief Define the dispersivity.
     *
     * \param element The finite element
     * \param fvElemGeom The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     */
    double dispersivity(const Element &element,
                    const FVElementGeometry &fvElemGeom,
                    int scvIdx) const
    {
        return 0;
    }

    bool useTwoPointGradient(const Element &elem,
                             int vertexI,
                             int vertexJ) const
    {
        return false;
    }

private:
    Scalar permeability_;
    Scalar porosity_;
    Scalar tortuosity_;
};

}

#endif
