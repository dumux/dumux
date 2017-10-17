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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 * \ingroup BoundaryCoupling
 * \brief @copybrief Dumux::CouplingMapper
 */

#ifndef DUMUX_COUPLINGMAPPER_HH
#define DUMUX_COUPLINGMAPPER_HH

#include <dumux/common/propertysystem.hh>

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(StokesProblemTypeTag);
NEW_PROP_TAG(DarcyProblemTypeTag);
NEW_PROP_TAG(DarcyProblem);
NEW_PROP_TAG(StokesProblem);
NEW_PROP_TAG(CouplingManager);
NEW_PROP_TAG(DarcyToStokesMapValue);
}

/*!
 * \ingroup BoundaryCoupling
 * \brief Mapper class to create coupling maps for both models.
 */
template<typename TypeTag>
class CouplingMapperStokesDarcy
{
    // extract some types from the actual problem
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);

    // obtain the type tags of the sub problems
    using StokesProblemTypeTag = typename GET_PROP_TYPE(TypeTag, StokesProblemTypeTag);
    using DarcyProblemTypeTag= typename GET_PROP_TYPE(TypeTag, DarcyProblemTypeTag);

    using StokesProblem = typename GET_PROP_TYPE(StokesProblemTypeTag, Problem);
    using DarcyProblem = typename GET_PROP_TYPE(DarcyProblemTypeTag, Problem);

    using StokesElementVolumeVariables = typename GET_PROP_TYPE(StokesProblemTypeTag, ElementVolumeVariables);
    using StokesFVElementGeometry = typename GET_PROP_TYPE(StokesProblemTypeTag, FVElementGeometry);
    using StokesSubControlVolume = typename GET_PROP_TYPE(StokesProblemTypeTag, SubControlVolume);

    using DarcyFVElementGeometry = typename GET_PROP_TYPE(DarcyProblemTypeTag, FVElementGeometry);

    using StokesGridView = typename GET_PROP_TYPE(StokesProblemTypeTag, GridView);
    using DarcyGridView = typename GET_PROP_TYPE(DarcyProblemTypeTag, GridView);

    using StokesElement = typename StokesGridView::template Codim<0>::Entity;
    using DarcyElement = typename DarcyGridView::template Codim<0>::Entity;

    using CoordScalar = typename StokesGridView::ctype;

    enum {
        stokesDim = StokesGridView::dimension,
        darcyDim = DarcyGridView::dimension,
        dimWorld = StokesGridView::dimensionworld
    };

    enum{maxNC = (stokesDim < 3 ? 4 : 8)};

    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using Vertex = typename StokesGridView::template Codim<stokesDim>::Entity;
    using DarcyVertex = typename DarcyGridView::template Codim<darcyDim>::Entity;

    // aliases for the stokes grid ansatz function shape for darcyToStokesMapValue
    using ReferenceElement = typename Dune::ReferenceElement<CoordScalar, stokesDim>;
    using ReferenceElements = typename Dune::ReferenceElements<CoordScalar, stokesDim>;

    using DarcyToStokesMapValue = typename GET_PROP_TYPE(TypeTag, DarcyToStokesMapValue);
    using DarcyToStokesMap = std::map<unsigned int, DarcyToStokesMapValue>; // key: global darcy dof index (vIdx)

    // aliases for the stokesToDarcy maps
    struct StokesToDarcyMapValue
    {
        unsigned int darcyDofIdx;
        unsigned int darcyElementIdx;
    };
    //    using StokesCCToDarcyMap = std::map<unsigned int, std::vector<StokesToDarcyMapValue>>; // key: global scv index (eIdx, scvIdx)
    //    using StokesFaceToDarcyMap = std::map<unsigned int, std::vector<StokesToDarcyMapValue>>; // key: global scv index (eIdx, scvIdx)

    // TODO only one Darcy element mapped to each Stokes element --> no vector needed?! (changes in couplingmanager necessary!)
    using StokesCCToDarcyMap = std::map<unsigned int, StokesToDarcyMapValue>; // key: global scv index (eIdx, scvIdx)
    using StokesFaceToDarcyMap = std::map<unsigned int, StokesToDarcyMapValue>; // key: global scv index (eIdx, scvIdx)

public:
    CouplingMapperStokesDarcy(StokesProblem &stokesProblem, DarcyProblem &darcyProblem, CouplingManager &couplingManager)
: stokesGridView_(stokesProblem.gridView()),
  darcyGridView_(darcyProblem.gridView()),
  stokesProblem_(stokesProblem),
  darcyProblem_(darcyProblem),
  couplingManager_(couplingManager)
{}

    /*!
     * \brief Computes the coupling maps
     */
    void computeCouplingMaps()
    {
        // get the Stokes boundingBoxTree
        const auto &stokesTree = stokesProblem_.boundingBoxTree();

        //iterate over all Darcy elements
        for(const auto &darcyElement : elements(darcyGridView_))
        {
            DarcyFVElementGeometry darcyFVElementGeometry = localView(darcyProblem_.model().globalFvGeometry());
            darcyFVElementGeometry.bind(darcyElement);

            const auto darcyElementIdx = darcyGridView_.indexSet().index(darcyElement);

            // iterate over the Darcy element's scvs
            for(auto&& darcyScv : scvs(darcyFVElementGeometry))
            {
                // create a unique global vertex (=DOF) index for the Darcy model vertex
                const unsigned int darcyDofIdxGlobal = darcyScv.dofIndex();

                // only consider Darcy DOFS on the Darcy domain's boundary
                if(!darcyProblem_.model().onBoundary(darcyDofIdxGlobal))
                    continue;

                // if the Darcy scv is at the coupling interface determine the respective scvf
                auto darcyPos = darcyScv.center();

                for(auto& darcyScvf : scvfs(darcyFVElementGeometry))
                {
                    if(darcyProblem_.onCouplingInterface(darcyScvf.center()))
                    {
                        darcyPos = darcyScvf.center();
                        continue; // assumption: coupling interface is a straight line, no corners
                    }
                }

                // determine the Stokes element that is coupled to the Darcy DOFs
                std::vector<unsigned int> stokesElementIndices = stokesTree.computeEntityCollisions(darcyPos);

                // do nothing if there are no coupled elements
                if(stokesElementIndices.empty())
                    continue;

                const int stokesElementIndex = stokesElementIndices[0]; // only one Stokes element coupled to one Darcy element

                // set the Darcy eIdx
                darcyToStokesMap_[darcyDofIdxGlobal].setDarcyElementIndex(darcyElementIdx);
                darcyToStokesMap_[darcyDofIdxGlobal].setDarcyScvIdx(darcyScv.dofIndex());

                const auto& stokesElement = stokesTree.entity(stokesElementIndex);
                StokesFVElementGeometry stokesFvGeometry = localView(stokesProblem_.model().globalFvGeometry());
                stokesFvGeometry.bind(stokesElement);

                // loop over all Stokes sub control volumes and check if the Darcy vertex is inside
                for(auto&& stokesScv : scvs(stokesFvGeometry))
                {
                    // Darcy vertex is inside Stokes scv
                    // create a unique index for the Stokes scv
                    const unsigned int stokesCCDofIdxGlobal = stokesScv.dofIndex();

                    StokesToDarcyMapValue value;
                    value.darcyDofIdx = darcyDofIdxGlobal;
                    value.darcyElementIdx = darcyElementIdx;
                    stokesCCToDarcyMap_[stokesCCDofIdxGlobal] = value;

                    darcyToStokesMap_[darcyDofIdxGlobal].addStokesElementIndex(stokesElementIndex);
                    darcyToStokesMap_[darcyDofIdxGlobal].addStokesCCDofIndex(stokesCCDofIdxGlobal);

                    for(auto&& stokesScvf : scvfs(stokesFvGeometry))
                    {
                        // make sure that the faces and the the Darcy vertices lie within the same
                        if(stokesScvf.boundary())
                        {
                            const auto delta = stokesScvf.center() - darcyPos;
                            if(delta[1] < 1e-8)
                            {
                                stokesFaceToDarcyMap_[stokesScvf.dofIndex()] = value;
                                darcyToStokesMap_[darcyDofIdxGlobal].addstokesFaceDofIndex(stokesScvf.dofIndex());
                            }
                        }
                    }
                }
            }
        }
    }

    /*!
     * \brief Returns a reference to the darcyToStokesMap
     */
    const DarcyToStokesMap &darcyToStokesMap() const
    { return darcyToStokesMap_; }

    /*!
     * \brief Returns a reference to the stokesCCToDarcyMap
     */
    const StokesCCToDarcyMap &stokesCCToDarcyMap() const
    { return stokesCCToDarcyMap_; }

    /*!
     * \brief Returns a reference to the stokesFaceToDarcyMap
     */
    const StokesFaceToDarcyMap &stokesFaceToDarcyMap() const
    { return stokesFaceToDarcyMap_; }

private:

    const StokesGridView stokesGridView_;
    const DarcyGridView darcyGridView_;

    StokesProblem &stokesProblem_;
    DarcyProblem &darcyProblem_;

    DarcyToStokesMap darcyToStokesMap_;

    StokesCCToDarcyMap stokesCCToDarcyMap_;
    StokesFaceToDarcyMap stokesFaceToDarcyMap_;

    CouplingManager &couplingManager_;
};

} // end namespace

#endif
