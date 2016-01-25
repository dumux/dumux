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
#ifndef DUMUX_INJECTION_SPATIAL_PARAMETERS_HH
#define DUMUX_INJECTION_SPATIAL_PARAMETERS_HH

#include <dumux/porousmediumflow/2pncmin/implicit/indices.hh>
#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{
//forward declaration
template<class TypeTag>
class DissolutionSpatialparams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(DissolutionSpatialparams);

// Set the spatial parameters
SET_TYPE_PROP(DissolutionSpatialparams, SpatialParams, Dumux::DissolutionSpatialparams<TypeTag>);

// Set the material Law
SET_PROP(DissolutionSpatialparams, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef RegularizedBrooksCorey<Scalar> EffMaterialLaw;
public:
    // define the material law parameterized by absolute saturations
    typedef EffToAbsLaw<EffMaterialLaw> type;
};
}

/**
 * \brief Definition of the spatial parameters for the brine-co2 problem
 *
 */
template<class TypeTag>
class DissolutionSpatialparams : public ImplicitSpatialParams<TypeTag>
{
    typedef ImplicitSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename Grid::ctype CoordScalar;
    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld,
    };

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        wPhaseIdx = FluidSystem::wPhaseIdx,
        nPhaseIdx = FluidSystem::nPhaseIdx,
    };

    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<CoordScalar, dimWorld, dimWorld> Tensor;

    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GridView::template Codim<0>::Entity Element;

    typedef std::vector<Scalar> PermeabilityType;
    typedef std::vector<MaterialLawParams> MaterialLawParamsVector;

    typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;

public:
    DissolutionSpatialparams(const GridView &gridView)
        : ParentType(gridView),
          K_(0)
    {
        //set main diagonal entries of the permeability tensor to a value
        //setting to one value means: isotropic, homogeneous
        for (int i = 0; i < dim; i++)
            K_[i][i] = 2.23e-14;

        // residual saturations
        materialParams_.setSwr(0.2);
        materialParams_.setSnr(1E-3);

        //parameters of Brooks & Corey Law
        materialParams_.setPe(500);
        materialParams_.setLambda(2);
    }

    ~DissolutionSpatialparams()
    {}
    /*!
     * \brief Update the spatial parameters with the flow solution
     *        after a timestep.
     *
     * \param globalSolution The global solution vector
     */
    void update(const SolutionVector &globalSolution)
    { };

    /*! Intrinsic permeability tensor K \f$[m^2]\f$ depending
     *  on the position in the domain
     *
     *  \param element The finite volume element
     *  \param fvGeometry The finite-volume geometry in the box scheme
     *  \param scvIdx The local vertex index
     *
     *  Alternatively, the function intrinsicPermeabilityAtPos(const GlobalPosition& globalPos)
     *  could be defined, where globalPos is the vector including the global coordinates
     *  of the finite volume.
     */
    const Dune::FieldMatrix<Scalar, dim, dim> &intrinsicPermeability(const Element &element,
                                                                     const FVElementGeometry &fvGeometry,
                                                                     const int scvIdx) const
    {
        return K_;
    }

    /*!
     * \brief Define the minimum porosity \f$[-]\f$ after salt precipitation
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    double porosityMin(const Element &element,
                       const FVElementGeometry &fvGeometry,
                       int scvIdx) const
     {
        return 1e-5;
     }

    /*!
     * \brief Define the minimum porosity \f$[-]\f$ after clogging caused by mineralization
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    double porosity(const Element &element,
                     const FVElementGeometry &fvGeometry,
                     int scvIdx) const
     {
        return 0.11;
     }


    double solidity(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    int scvIdx) const
    {

        return 1 - 0.11;
    }

    double SolubilityLimit() const
    {
        return 0.26;
    }

    double theta(const Element &element,
            const FVElementGeometry &fvGeometry,
            int scvIdx) const
    {
        return 10.0;
    }


    // return the brooks-corey context depending on the position
    const MaterialLawParams& materialLawParams(const Element &element,
                                                const FVElementGeometry &fvGeometry,
                                                int scvIdx) const
    {
        return materialParams_;
    }

private:
    Dune::FieldMatrix<Scalar, dim, dim> K_;
    MaterialLawParams materialParams_;
};

}

#endif
