// $Id$
/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff                                      *
 *   Copyright (C) 2007-2008 by Klaus Mosthaf                                *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
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
 * \brief spatial parameters for the RichardsLensProblem
 */
#ifndef DUMUX_RICHARDS_LENS_SPATIAL_PARAMETERS_HH
#define DUMUX_RICHARDS_LENS_SPATIAL_PARAMETERS_HH

#include <dumux/material/spatialparameters/boxspatialparameters.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{

/*!
 * \ingroup RichardsModel
 * \brief The spatial parameters for the RichardsLensProblem
 */
template<class TypeTag>
class RichardsLensSpatialParameters : public BoxSpatialParameters<TypeTag>
{
    typedef BoxSpatialParameters<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename Grid::ctype CoordScalar;

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld,
    };

    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    // define the material law which is parameterized by effective
    // saturations
    typedef RegularizedVanGenuchten<Scalar> EffectiveLaw;

public:
    /*!
     * \brief The material law to be used.
     *
     * This problem uses the RegularizedVanGenuchten material law.
     */
    typedef EffToAbsLaw<EffectiveLaw> MaterialLaw;
    //! The parameters of the material law to be used
    typedef typename MaterialLaw::Params MaterialLawParams;

    /*!
     * \brief Constructor
     *
     * \param gridView The DUNE GridView representing the spatial
     *                 domain of the problem.
     */
    RichardsLensSpatialParameters(const GridView& gridView)
        : ParentType(gridView)
    {
        // residual saturations
        lensMaterialParams_.setSwr(0.18);
        lensMaterialParams_.setSnr(0.0);
        outerMaterialParams_.setSwr(0.05);
        outerMaterialParams_.setSnr(0.0);

        // parameters for the Van Genuchten law
        // alpha and n
        lensMaterialParams_.setVgAlpha(0.00045);
        lensMaterialParams_.setVgN(7.3);
        outerMaterialParams_.setVgAlpha(0.0037);
        outerMaterialParams_.setVgN(4.7);

        // parameters for the linear law
        // minimum and maximum pressures
//        lensMaterialParams_.setEntryPC(0);
//        outerMaterialParams_.setEntryPC(0);
//        lensMaterialParams_.setMaxPC(0);
//        outerMaterialParams_.setMaxPC(0);

        lensK_ = 1e-12;
        outerK_ = 5e-12;
    }

    /*!
     * \brief Returns the intrinsic permeability tensor [m^2] at a given location
     *
     * \param element An arbitrary DUNE Codim<0> entity of the grid view
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    Scalar intrinsicPermeability(const Element &element,
                                 const FVElementGeometry &fvElemGeom,
                                 int scvIdx) const
    {
        const GlobalPosition &globalPos = fvElemGeom.subContVol[scvIdx].global;
        if (isInLens_(globalPos))
            return lensK_;
        return outerK_;
    }

    /*!
     * \brief Returns the porosity [] at a given location
     *
     * \param element An arbitrary DUNE Codim<0> entity of the grid view
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    Scalar porosity(const Element &element,
                    const FVElementGeometry &fvElemGeom,
                    int scvIdx) const
    { return 0.4; }

    /*!
     * \brief Returns the parameters for the material law at a given location
     *
     * \param element An arbitrary DUNE Codim<0> entity of the grid view
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    const MaterialLawParams& materialLawParams(const Element &element,
                                                const FVElementGeometry &fvElemGeom,
                                                int scvIdx) const
    {
        return materialLawParams(fvElemGeom.subContVol[scvIdx].global);
    }

    /*!
     * \brief Returns the parameters for the material law at a given location
     *
     * This method is not actually required by the Richards model, but provided
     * for the convenience of the RichardsLensProblem
     *
     * \param globalPos A global coordinate vector
     */
    const MaterialLawParams& materialLawParams(const GlobalPosition &globalPos) const
    {

        if (isInLens_(globalPos))
            return lensMaterialParams_;
        return outerMaterialParams_;
    }

    /*!
     * \brief Set the bounding box of the low-permeability lens
     *
     * This method is not actually required by the Richards model, but provided
     * for the convenience of the RichardsLensProblem
     *
     * \param lensLowerLeft the lower left corner of the lens
     * \param lensUpperRight the upper right corner of the lens
     */
    void setLensCoords(const GlobalPosition& lensLowerLeft,
                       const GlobalPosition& lensUpperRight)
    {
        lensLowerLeft_ = lensLowerLeft;
        lensUpperRight_ = lensUpperRight;
    }

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

    Scalar lensK_;
    Scalar outerK_;
    MaterialLawParams lensMaterialParams_;
    MaterialLawParams outerMaterialParams_;
};

} // end namespace

#endif

