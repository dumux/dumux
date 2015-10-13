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
/*!
 * \file
 * \brief Volume variables gathered on an element
 */
#ifndef DUMUX_BOX_EL1P2C_ELEMENT_VOLUME_VARIABLES_HH
#define DUMUX_BOX_EL1P2C_ELEMENT_VOLUME_VARIABLES_HH

#include <dumux/implicit/box/boxproperties.hh>
#include <dumux/implicit/box/boxelementvolumevariables.hh>

namespace Dumux
{

/*!
 * \ingroup ElOnePTwoCBoxModel
 *
 * \brief This class stores an array of VolumeVariables objects, one
 *        volume variables object for each of the element's vertices
 */
template<class TypeTag>
class ElOnePTwoCElementVolumeVariables : public BoxElementVolumeVariables<TypeTag>
{
    typedef BoxElementVolumeVariables<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    enum { dim = GridView::dimension };

    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;

public:
    /*!
     * \brief The constructor.
     */
    ElOnePTwoCElementVolumeVariables()
    { }

    /*!
     * \brief Construct the volume variables for all vertices of an element.
     *
     * \param problem The problem which needs to be simulated.
     * \param element The DUNE Codim<0> entity for which the volume variables ought to be calculated
     * \param fvGeometry The finite volume geometry of the element
     * \param oldSol Tells whether the model's previous or current solution should be used.
     *
     * This class is required for the update of the effective porosity values at the
     * vertices since it is a function of the divergence of the solid displacement
     * at the integration points
     */
    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                bool oldSol)
    {
        ParentType::update(problem, element, fvGeometry, oldSol);
        this->updateEffPorosity(problem, element, fvGeometry);

    };

    /*!
     * \brief Update the effective porosities and the volumetric strain divU for all vertices of an element.
     *
     * \param problem The problem which needs to be simulated.
     * \param element The DUNE Codim<0> entity for which the volume variables ought to be calculated
     * \param fvGeometry The finite volume geometry of the element
     *
     * This function is required for the update of the effective porosity / divU values at the
     * vertices.
     *
     * During the partial derivative calculation, changes of the solid displacement
     * at vertex i can affect effective porosities / divU of all element vertices.
     * To correctly update the effective porosities / divU of all element vertices
     * an iteration over all scv faces is required.
     * The remaining volvars are only updated for the vertex whose primary variable
     * is changed for the derivative calculation.
     */
    void updateEffPorosity(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry)
    {
        // we assert that the i-th shape function is
        // associated to the i-th vert of the element.
        int numScv = element.subEntities(dim);

        // number of faces which contribute to the porosity value in the sub-control volume
        std::vector<double>  numContributingFaces;
        numContributingFaces.resize(numScv);

        for (int scvIdx = 0; scvIdx < numScv; scvIdx++) {
            (*this)[scvIdx].effPorosity = 0.0;
            (*this)[scvIdx].divU = 0.0;
            numContributingFaces[scvIdx] = 0.0;
        }
        for (int fIdx = 0; fIdx < fvGeometry.numScvf; fIdx++)
        {
            // evaluate the gradients at the IPs for each subcontrol volume face
            FluxVariables fluxVars(problem,
                               element,
                               fvGeometry,
                               fIdx,
                               *this);

            numContributingFaces[fluxVars.face().i] += 1;
            numContributingFaces[fluxVars.face().j] += 1;

            // average value for the effective porosity
            (*this)[fluxVars.face().i].effPorosity += fluxVars.effPorosity();
            (*this)[fluxVars.face().j].effPorosity += fluxVars.effPorosity();
            // average value for the volumetric strain
            (*this)[fluxVars.face().i].divU += fluxVars.divU();
            (*this)[fluxVars.face().j].divU += fluxVars.divU();

        }
        for (int scvIdx = 0; scvIdx < numScv; scvIdx++) {
            (*this)[scvIdx].effPorosity /= numContributingFaces[scvIdx];
            (*this)[scvIdx].divU /= numContributingFaces[scvIdx];
        }
    };


};

} // namespace Dumux

#endif
