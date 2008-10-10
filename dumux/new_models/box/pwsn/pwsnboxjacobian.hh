/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUMUX_PWSN_BOX_JACOBIAN_IMP_HH
#define DUMUX_PWSN_BOX_JACOBIAN_IMP_HH

#include <dumux/new_models/box/boxjacobian.hh>
#include <dumux/fvgeometry/fvelementgeometry.hh>

#include <boost/format.hpp>

namespace Dune
{
    /*!
     * \brief pw-Sn formulation specific details needed to calculate
     *        the jacobian in the BOX scheme numerically.
     *
     * This class is used to fill the gaps for BoxJacobian, not
     * standalone. Please use PwSnBoxJacobian to get it nicely
     * wrapped...
     */
    template<class ProblemT>
    class _PwSnBoxJacobianImp
    {
        typedef ProblemT                                Problem;
        typedef typename Problem::DomainTraits          DomTraits;
        typedef typename Problem::BoxTraits             BoxTraits;
        
        typedef typename DomTraits::Scalar             Scalar;
        typedef typename DomTraits::CoordScalar        CoordScalar;
        typedef typename DomTraits::Grid               Grid;
        typedef typename DomTraits::Cell               Cell;
        typedef typename DomTraits::LocalCoord         LocalCoord;

        enum {
            GridDim     = DomTraits::GridDim,
            WorldDim    = DomTraits::WorldDim,

            NumUnknowns = BoxTraits::NumUnknowns,
            PwIndex     = BoxTraits::PwIndex,
            SnIndex     = BoxTraits::SnIndex
        };

        typedef typename BoxTraits::UnknownsVector     UnknownsVector;
        typedef typename BoxTraits::FVElementGeometry  FVElementGeometry;
        typedef typename BoxTraits::CachedCellData       CachedCellData;
        typedef typename BoxTraits::CachedSubContVolData CachedSubContVolData;
        typedef typename BoxTraits::LocalFunction        LocalFunction;

    public:
        /*!
         * \brief Evaluate the rate of change of mass within a sub
         *        control volume of a dual cell in the pw-Sn
         *        formulation.
         */
        static void evalMassBalance(UnknownsVector &result,
                                    const Problem &problem,
                                    const Cell &cell,
                                    const FVElementGeometry &dualCell,
                                    const LocalFunction &localSol,
                                    const int scvId)
            {
                Scalar porosity = problem.porosity(cell);
                // derivative of the wetting phase mass regarding time
                result[PwIndex] = -problem.densityW()
                                   *porosity
                                   *localSol.atSubContVol[scvId][SnIndex];
                // derivative of the non-wetting phase mass regarding time
                result[SnIndex] = problem.densityN()
                                  *porosity
                                  *localSol.atSubContVol[scvId][SnIndex];
            }

        
        /*!
         * \brief Evaluates the mass flux over a face of a subcontrol
         *        volume.
         */
        static void evalMassFlux(UnknownsVector &flux,
                                 const Problem &problem,
                                 const Cell& cell,
                                 const FVElementGeometry &dualCell,
                                 const LocalFunction &localSol,
                                 const CachedCellData &cachedData,
                                 const int faceId)
            {
                Api::require<Api::BasicDomainTraits, typename ProblemT::DomainTraits>();
                assert(NumUnknowns == 2);

                const typename FVElementGeometry::SubControlVolumeFace
                    &face = dualCell.subContVolFace[faceId];
                const int i = face.i;
                const int j = face.j;

                // permeability in edge direction
                LocalCoord Kij(0);

                // Kij += K*normal
                problem.permeability(cell).umv(face.normal, Kij);

                for (int phase = 0; phase < NumUnknowns; phase++) {
                    // calculate FE gradient
                    LocalCoord pGrad(0);
                    for (int k = 0; k < dualCell.nNodes; k++) {
                        LocalCoord grad(face.grad[k]);
                        if (phase == SnIndex)
                            grad *= cachedData.atSubContVol[k].pN;
                        else
                            grad *= localSol.atSubContVol[k][PwIndex];

                        pGrad += grad;
                    }

                    // adjust pressure gradient by gravity force
                    Scalar phaseDensity = problem.density(phase);
                    LocalCoord gravity = problem.gravity();
                    gravity *= phaseDensity;
                    pGrad   -= gravity;
                    
                    // calculate the flux using upwind
                    Scalar outward = pGrad*Kij;
                    if (outward < 0)
                        flux[phase] = phaseDensity*cachedData.atSubContVol[i].mobility[phase]*outward;
                    else
                        flux[phase] = phaseDensity*cachedData.atSubContVol[j].mobility[phase]*outward;
                }
            }

        /*!
         * \brief Pre-compute the cell cache data.
         *
         * This method is called by BoxJacobian (which in turn is
         * called by the operator assembler) every time the current
         * cell changes.
         */
        static void updateCellCache(CachedCellData &dest,
                                    const Problem &problem,
                                    const Cell& cell,
                                    const FVElementGeometry &dualCell,
                                    const LocalFunction &sol)
            {
                assert(NumUnknowns == 2);

                int cellIndex   = problem.cellIndex(cell);
                for (int i = 0; i < dualCell.nNodes; i++) {
                    int iGlobal = problem.vertexIndex(cell, i);
                    _partialUpdateCellCache(dest,
                                            problem,
                                            cell,
                                            cellIndex,
                                            dualCell,
                                            sol,
                                            i,        // index of sub volume to update,
                                            iGlobal); // global vertex index of the sub volume's grid vertex
                }
            }

        /*!
         * \brief Update all pre-computed data which is dependend of
         *        the solution at an individual vertex.
         *
         * This method is called by BoxJacobian (which in turn is
         * called by the operator assembler) every time the current
         * cell changes.
         */
        static void partialUpdateCellCache(CachedCellData &dest,
                                           const CachedCellData &source,
                                           const Problem &problem,
                                           const Cell& cell,
                                           const FVElementGeometry &dualCell,
                                           const LocalFunction &sol,
                                           int subVolIndex) // ID of the subvolume/grid vertex
            {
                assert(NumUnknowns == 2);

                int cellIndex   = problem.cellIndex(cell);
                int globalVertexIndex = problem.vertexIndex(cell, subVolIndex);

                // copy source cache into the destination cache
                dest = source;
                // evaluate cached data at subvolume
                _partialUpdateCellCache(dest,
                                        problem,
                                        cell,
                                        cellIndex,
                                        dualCell,
                                        sol,
                                        subVolIndex,
                                        globalVertexIndex);
            }

    private:
        static void _partialUpdateCellCache(CachedCellData &dest,
                                            const Problem &problem,
                                            const Cell& cell,
                                            int         cellIndex,
                                            const FVElementGeometry &dualCell,
                                            const LocalFunction &sol,
                                            int i, // index of the subvolume/grid vertex
                                            int iGlobal) // global index of the sub-volume's vertex
            {
                // Current cache at sub-controlvolume
                CachedSubContVolData &curSCVCache = dest.atSubContVol[i];
                // Current solution for sub-controlvolume
                const UnknownsVector &curSCVSol = sol.atSubContVol[i];

                curSCVCache.Sw = 1.0 - curSCVSol[SnIndex];
                curSCVCache.pC = problem.pC(cell,
                                            cellIndex,
                                            i,
                                            iGlobal,
                                            curSCVCache.Sw);
                curSCVCache.pN = curSCVSol[PwIndex] + curSCVCache.pC;
                curSCVCache.mobility[PwIndex] = problem.mobilityW(cell,
                                                                  cellIndex,
                                                                  i,
                                                                  iGlobal,
                                                                  curSCVCache.Sw);
                curSCVCache.mobility[SnIndex] = problem.mobilityN(cell,
                                                                  cellIndex,
                                                                  i,
                                                                  iGlobal,
                                                                  curSCVSol[SnIndex]);

/*                if (HACKY_HACK) {
                    std::cout << boost::format("_partialUpdateCellCache\n");
                    std::cout << boost::format("cellIndex: %d\n")%cellIndex;
                    std::cout << boost::format("i: %i\n")%i;
                    std::cout << boost::format("iGlobal: %i\n")%iGlobal;
                    std::cout << boost::format("---> pW: %g\n")%curSCVSol[PwIndex];
                    std::cout << boost::format("---> Sn: %g\n")%curSCVSol[SnIndex];
                    std::cout << boost::format("Sw: %g\n")%curSCVCache.Sw;
                    std::cout << boost::format("pC: %g\n")%curSCVCache.pC;
                    std::cout << boost::format("pN: %g\n")%curSCVCache.pN;
                    std::cout << boost::format("mobility[PwIndex]: %g\n")%curSCVCache.mobility[PwIndex];
                    std::cout << boost::format("mobility[SnIndex]: %g\n")%curSCVCache.mobility[SnIndex];
                    std::cout << boost::format("end _partialUpdateCellCache\n");
                    }
*/
            }

    };


    /*!
     * \brief The local jacobian for the BOX model in the Pw-Sn formulation.
     */
    template<class ProblemT, class BoxTraitsT>
    class PwSnBoxJacobian : public BoxJacobian<ProblemT, BoxTraitsT, _PwSnBoxJacobianImp<ProblemT> >
    {
    private:
        typedef BoxJacobian<ProblemT, BoxTraitsT, _PwSnBoxJacobianImp<ProblemT> >   ParentType;

    public:
        PwSnBoxJacobian(ProblemT &problem,
                        bool levelBoundaryAsDirichlet,
                        bool procBoundaryAsDirichlet=true) 
            : ParentType(problem,
                         levelBoundaryAsDirichlet,
                         procBoundaryAsDirichlet)
            {};
    };
}

#endif
