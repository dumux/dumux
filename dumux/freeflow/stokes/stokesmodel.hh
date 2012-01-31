// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
 *   Copyright (C) 2010 by Katherina Baber, Klaus Mosthaf                    *
 *   Copyright (C) 2008-2009 by Bernd Flemisch, Andreas Lauser               *
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
#ifndef DUMUX_STOKES_MODEL_HH
#define DUMUX_STOKES_MODEL_HH

#include "stokeslocalresidual.hh"
#include "stokesnewtoncontroller.hh"
#include "stokeslocaljacobian.hh"
#include "stokesproblem.hh"
#include "stokesproperties.hh"

namespace Dumux
{
/*!
 * \ingroup BoxStokesModel
 * \brief Adaption of the BOX scheme to the stokes flow model.
 *
 * This model implements laminar Stokes flow of a single fluid, solving a momentum balance:
 * \f[
\frac{\partial \left(\varrho_g {\boldsymbol{v}}_g\right)}{\partial t}
+ \boldsymbol{\nabla} \boldsymbol{\cdot} \left(p_g {\bf {I}}
- \mu_g \left(\boldsymbol{\nabla} \boldsymbol{v}_g
+ \boldsymbol{\nabla} \boldsymbol{v}_g^T\right)\right)
- \varrho_g {\bf g} = 0,
 * \f]
 *
 * and the mass balance equation:
 * \f[
\frac{\partial \varrho_g}{\partial t} + \boldsymbol{\nabla}\boldsymbol{\cdot}\left(\varrho_g {\boldsymbol{v}}_g\right) - q_g = 0
 * \f]
 *
 * This is discretized by a fully-coupled vertex- centered finite volume
 * (box) scheme as spatial and the implicit Euler method
 * as temporal discretization.
 */
template<class TypeTag >
class StokesModel : public BoxModel<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

public:
    /*!
     * \brief Calculate the fluxes across a certain layer in the domain.
     * The layer is situated perpendicular to the coordinate axis "coord" and cuts
     * the axis at the value "coordValue"
     *
     */
    void calculateFluxAcrossLayer (const SolutionVector &globalSol, Dune::FieldVector<Scalar, 2> &flux, int coord, Scalar coordVal)
    {
        GlobalPosition globalI, globalJ;
        PrimaryVariables tmpFlux(0.0);
        int sign;

        FVElementGeometry fvElemGeom;
        ElementVolumeVariables elemVolVars;

        // Loop over elements
        ElementIterator elemIt = this->problem_.gridView().template begin<0>();
        ElementIterator endit = this->problem_.gridView().template end<0>();
        for (; elemIt != endit; ++elemIt)
        {
            if (elemIt->partitionType() != Dune::InteriorEntity)
                continue;

            fvElemGeom.update(this->gridView_(), *elemIt);
            elemVolVars.update(this->problem_(), *elemIt, fvElemGeom);
            this->localResidual().evalFluxes(*elemIt, elemVolVars);

            bool hasLeft = false;
            bool hasRight = false;
            for (int i = 0; i < fvElemGeom.numVertices; i++) {
                const GlobalPosition &global = fvElemGeom.subContVol[i].global;
                if (globalI[coord] < coordVal)
                    hasLeft = true;
                else if (globalI[coord] >= coordVal)
                    hasRight = true;
            }
            if (!hasLeft || !hasRight)
                continue;

            for (int i = 0; i < fvElemGeom.numVertices; i++) {
                const int &globalIdx =
                    this->vertexMapper().map(*elemIt, i, dim);
                const GlobalPosition &global = fvElemGeom.subContVol[i].global;
                if (globalI[coord] < coordVal)
                    flux += this->localResidual().residual(i);
            }
        }

        flux = this->problem_.gridView().comm().sum(flux);
    }

    /*!
     * \brief Calculate mass in the whole model domain
     *         and get minimum and maximum values of primary variables
     *
     */
    void calculateMass(const SolutionVector &sol, Dune::FieldVector<Scalar, 2> &mass)
    {
        //         const DofMapper &dofMapper = this->dofEntityMapper();
        VolumeVariables tmp;
        Scalar vol, poro, rhoN, rhoW, satN, satW, pW;//, Te;
        Scalar massNPhase(0.), massWPhase(0.);

        mass = 0;
        Scalar minSat = 1e100;
        Scalar maxSat = -1e100;
        Scalar minP = 1e100;
        Scalar maxP = -1e100;
        //         Scalar minTe = 1e100;
        //         Scalar maxTe = -1e100;

        FVElementGeometry fvElemGeom;
        VolumeVariables volVars;
        ElementBoundaryTypes elemBcTypes;

        // Loop over elements
        ElementIterator elemIt = this->problem_.gridView().template begin<0>();
        ElementIterator endit = this->problem_.gridView().template end<0>();
        for (; elemIt != endit; ++elemIt)
        {
            if (elemIt->partitionType() != Dune::InteriorEntity)
                continue;

            fvElemGeom.update(this->gridView_(), *elemIt);
            elemBcTypes.update(this->problem_(), *elemIt, fvElemGeom);

            int numVerts = elemIt->template count<dim>();
            for (int i = 0; i < numVerts; ++i)
            {
                int globalIdx = this->vertexMapper().map(*elemIt, i, dim);
                volVars.update(sol[globalIdx],
                               this->problem_(),
                               *elemIt,
                               fvElemGeom,
                               i,
                               false);

                //                 int globalIdx = dofMapper.map(*elemIt, i, dim);
                vol = fvElemGeom.subContVol[i].volume;
                poro = volVars.porosity;
                rhoN = volVars.density;
                pW = volVars.pressure;
                //                 Te = asImp_()->temperature((*sol)[globalIdx]);


                massNPhase = vol * poro * rhoN;

                // get minimum and maximum values of primary variables
                minP = std::min(minP, pW);
                maxP = std::max(maxP, pW);
                //                 minTe = std::min(minTe, Te);
                //                 maxTe = std::max(maxTe, Te);

                // calculate total mass
                mass[0] += massNPhase; // mass nonwetting phase
            }
        }

        // IF PARALLEL: calculate total mass including all processors
        // also works for sequential calculation
        mass = this->problem_.gridView().comm().sum(mass);

        if(this->problem_.gridView().comm().rank() == 0) // IF PARALLEL: only print by processor with rank() == 0
        {
            // print minimum and maximum values
            std::cout << "nonwetting phase saturation: min = "<< minSat
                      << ", max = "<< maxSat << std::endl;
            std::cout << "wetting phase pressure: min = "<< minP
                      << ", max = "<< maxP << std::endl;
            //             std::cout << "temperature: min = "<< minTe
            //             << ", max = "<< maxTe << std::endl;
        }
    }


    /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     */
    template <class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol, MultiWriter &writer)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, dim> > VelocityField;

        // create the required scalar fields
        unsigned numVertices = this->gridView_().size(dim);
        ScalarField &pN = *writer.allocateManagedBuffer(numVertices);
        ScalarField &delP = *writer.allocateManagedBuffer(numVertices);
        ScalarField &rho = *writer.allocateManagedBuffer(numVertices);
        ScalarField &mu = *writer.allocateManagedBuffer(numVertices);
        VelocityField &velocity = *writer.template allocateManagedBuffer<Scalar, dim> (numVertices);

        unsigned numElements = this->gridView_().size(0);
        ScalarField &rank = *writer.allocateManagedBuffer(numElements);

        FVElementGeometry fvElemGeom;
        VolumeVariables volVars;
        ElementBoundaryTypes elemBcTypes;

        ElementIterator elemIt = this->gridView_().template begin<0>();
        ElementIterator endit = this->gridView_().template end<0>();
        for (; elemIt != endit; ++elemIt)
        {
            int idx = this->elementMapper().map(*elemIt);
            rank[idx] = this->gridView_().comm().rank();

            fvElemGeom.update(this->gridView_(), *elemIt);
            elemBcTypes.update(this->problem_(), *elemIt, fvElemGeom);

            int numLocalVerts = elemIt->template count<dim>();
            for (int i = 0; i < numLocalVerts; ++i)
            {
                int globalIdx = this->vertexMapper().map(*elemIt, i, dim);
                volVars.update(sol[globalIdx],
                               this->problem_(),
                               *elemIt,
                               fvElemGeom,
                               i,
                               false);

                pN[globalIdx] = volVars.pressure();
                delP[globalIdx] = volVars.pressure() - 1e5;
                rho[globalIdx] = volVars.density();
                mu[globalIdx] = volVars.viscosity();
                velocity[globalIdx] = volVars.velocity();
            };
        }
        writer.attachVertexData(pN, "P");
        writer.attachVertexData(delP, "delP");
        writer.attachVertexData(rho, "rho");
        writer.attachVertexData(mu, "mu");
        writer.attachVertexData(velocity, "v", dim);
    }
};
}

#include "stokespropertydefaults.hh"

#endif
