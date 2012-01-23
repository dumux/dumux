// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff                                      *
 *   Copyright (C) 2007-2010 by Yufei Cao                                    *
 *   Institute of Applied Analysis and Numerical Simulation                  *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@mathematik.uni-stuttgart.de                   *
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
#ifndef DUMUX_FVMPFAOPRESSURE2P_HH
#define DUMUX_FVMPFAOPRESSURE2P_HH

// dumux environment
#include <dumux/decoupled/common/fv/fvpressure.hh>

/**
 * @file
 * @brief  Finite Volume MPFA O-method discretization of a pressure equation.
 * @brief  Remark1: only for 2-D quadrilateral grid.
 * @brief  Remark2: can use UGGrid or SGrid (YaspGrid).
 * @brief  Remark3: without capillary pressure and gravity!
 * @author Yufei Cao
 */

namespace Dumux
{
/*! \ingroup FV2p
 *
 * \brief MPFA-O method for the pressure equation
 *
 * An interface for defining a numerical diffusion model for the
 *  solution of equations of the form
 * \f$ - \text{div}\, (\lambda K \text{grad}\, p ) = 0, \f$,
 * \f$p = g\f$ on \f$\Gamma_1\f$, and
 * \f$-\lambda K \text{grad}\, p \cdot \mathbf{n} = J\f$
 * on \f$\Gamma_2\f$. Here,
 * \f$p\f$ denotes the pressure, \f$K\f$ the absolute permeability,
 * and \f$\lambda\f$ the total mobility, possibly depending on the
 * saturation.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class FVMPFAOPressure2P: public FVPressure<TypeTag>
{
    typedef FVPressure<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) SpatialParameters;
    typedef typename SpatialParameters::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, TwoPIndices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;

    typedef typename SolutionTypes::ScalarSolution ScalarSolutionType;

    typedef typename GET_PROP_TYPE(TypeTag, GridTypeIndices) GridTypeIndices;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum
    {
        Sw = Indices::saturationW,
        Sn = Indices::saturationNW,
        pGlobal = Indices::pressureGlobal
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        pressEqIdx = Indices::pressEqIdx,
        satEqIdx = Indices::satEqIdx,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> FieldMatrix;

    typedef typename GET_PROP_TYPE(TypeTag, PressureCoefficientMatrix) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, PressureRHSVector) Vector;
    typedef Dune::FieldVector<Scalar, dim> FieldVector;

    //initializes the matrix to store the system of equations
    void initializeMatrix();

    //function which assembles the system of equations to be solved
    void assemble();

protected:
    Problem& problem()
    {
        return problem_;
    }

    const Problem& problem() const
    {
        return problem_;
    }

public:

    //constitutive functions are initialized and stored in the variables object
    //! updates and stores constitutive relations
    void updateMaterialLaws();

    void initialize(bool solveTwice = true)
    {
        initializeMatrix();
        updateMaterialLaws();

        assemble();
        this->solve();

        storePressureSolution();

        return;
    }

    //! updates the pressure field (analog to update function in Dumux::IMPET)
    void update()
    {
        assemble();
        this->solve();

        storePressureSolution();

        return;
    }

    void storePressureSolution()
    {
        int size = problem_.gridView().size(0);
        for (int i = 0; i < size; i++)
        {
            CellData& cellData = problem_.variables().cellData(i);
            storePressureSolution(i, cellData);
            cellData.fluxData().resetVelocity();
        }
    }

    void storePressureSolution(int globalIdx, CellData& cellData)
    {
            cellData.setGlobalPressure(this->pressure()[globalIdx]);
    }

    //! \brief Write data files
     /*  \param name file name */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        ScalarSolutionType *pressure = writer.allocateManagedBuffer (problem_.gridView().size(0));

        *pressure = this->pressure();

        writer.attachCellData(*pressure, "global pressure");

        // output  phase-dependent stuff
        ScalarSolutionType *pC = writer.allocateManagedBuffer (problem_.gridView().size(0));
        int size = problem_.gridView().size(0);
        for (int i = 0; i < size; i++)
        {
            CellData& cellData = problem_.variables().cellData(i);
            (*pC)[i] = cellData.capillaryPressure();
        }
        writer.attachCellData(*pC, "capillary pressure");

        return;
    }

    //! Constructs a FVMPFAOPressure2P object
    /**
     * \param problem a problem class object
     */
    FVMPFAOPressure2P(Problem& problem)
    : ParentType(problem), problem_(problem)
    {
        if (pressureType_ != pGlobal)
        {
            DUNE_THROW(Dune::NotImplemented, "Pressure type not supported!");
        }

        const Element& element = *(problem_.gridView().template begin<0> ());
        FluidState fluidState;
        fluidState.setPressure(wPhaseIdx, problem_.referencePressure(element));
        fluidState.setPressure(nPhaseIdx, problem_.referencePressure(element));
        fluidState.setTemperature(problem_.temperature(element));
        fluidState.setSaturation(wPhaseIdx, 1.);
        fluidState.setSaturation(nPhaseIdx, 0.);
        density_[wPhaseIdx] = FluidSystem::density(fluidState, wPhaseIdx);
        density_[nPhaseIdx] = FluidSystem::density(fluidState, nPhaseIdx);
        viscosity_[wPhaseIdx] = FluidSystem::viscosity(fluidState, wPhaseIdx);
        viscosity_[nPhaseIdx] = FluidSystem::viscosity(fluidState, nPhaseIdx);
    }

private:
    Problem& problem_;
    Scalar density_[numPhases];
    Scalar viscosity_[numPhases];

protected:
    static const int pressureType_ = GET_PROP_VALUE(TypeTag, PressureFormulation); //!< gives kind of pressure used (\f$p_w\f$, \f$p_n\f$, \f$p_{global}\f$)
    static const int saturationType_ = GET_PROP_VALUE(TypeTag, SaturationFormulation); //!< gives kind of saturation used (\f$S_w\f$, \f$S_n\f$)
};

template<class TypeTag>
void FVMPFAOPressure2P<TypeTag>::initializeMatrix()
{
    // determine matrix row sizes
    ElementIterator eItBegin = problem_.gridView().template begin<0>();
    ElementIterator eItEnd = problem_.gridView().template end<0>();
    for (ElementIterator eIt = eItBegin; eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = problem_.variables().index(*eIt);

        // initialize row size
        int rowSize = 1;

        // run through all intersections with neighbors
        IntersectionIterator isItBegin = problem_.gridView().ibegin(*eIt);
        IntersectionIterator isItEnd = problem_.gridView().iend(*eIt);
        for (IntersectionIterator isIt = isItBegin; isIt!=isItEnd; ++isIt)
        {
            IntersectionIterator tempisIt = isIt;
            IntersectionIterator tempisItBegin = isItBegin;

            // 'nextisIt' iterates over next codimension 1 intersection neighboring with 'isIt'
            IntersectionIterator nextisIt = ++tempisIt;

            // get 'nextisIt'
            switch (GET_PROP_VALUE(TypeTag, GridImplementation))
            {
                // for SGrid
                case GridTypeIndices::sGrid:
                {
                    if (nextisIt == isItEnd)
                    nextisIt = isItBegin;
                    else
                    {
                        nextisIt = ++tempisIt;

                        if (nextisIt == isItEnd)
                        {
                            nextisIt = ++tempisItBegin;
                        }
                    }

                    break;
                }
                // for YaspGrid
                case GridTypeIndices::yaspGrid:
                {
                    if (nextisIt == isItEnd)
                    {
                        nextisIt = isItBegin;
                    }
                    else
                    {
                        nextisIt = ++tempisIt;

                        if (nextisIt == isItEnd)
                        {
                            nextisIt = ++tempisItBegin;
                        }
                    }

                    break;
                }
                // for UGGrid
                case GridTypeIndices::ugGrid:
                {
                    if (nextisIt == isItEnd)
                    nextisIt = isItBegin;

                    break;
                }
                default:
                {
                    DUNE_THROW(Dune::NotImplemented, "GridType can not be used with MPFAO implementation!");
                    break;
                }
            }

            if (isIt->neighbor())
            rowSize++;

            if (isIt->neighbor() && nextisIt->neighbor())
            rowSize++;
        } // end of 'for' IntersectionIterator

        // set number of indices in row globalIdxI to rowSize
        this->A_.setrowsize(globalIdxI, rowSize);

    } // end of 'for' ElementIterator

    // indicate that size of all rows is defined
    this->A_.endrowsizes();

    // determine position of matrix entries
    for (ElementIterator eIt = eItBegin; eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = problem_.variables().index(*eIt);

        // add diagonal index
        this->A_.addindex(globalIdxI, globalIdxI);

        // run through all intersections with neighbors
        IntersectionIterator isItBegin = problem_.gridView().ibegin(*eIt);
        IntersectionIterator isItEnd = problem_.gridView().iend(*eIt);
        for (IntersectionIterator isIt = isItBegin; isIt!=isItEnd; ++isIt)
        {
            IntersectionIterator tempisIt = isIt;
            IntersectionIterator tempisItBegin = isItBegin;

            // 'nextisIt' iterates over next codimension 1 intersection neighboring with 'isIt'
            // sequence of next is anticlockwise of 'isIt'
            IntersectionIterator nextisIt = ++tempisIt;

            // get 'nextisIt'
            switch (GET_PROP_VALUE(TypeTag, GridImplementation))
            {
                // for SGrid
                case GridTypeIndices::sGrid:
                {
                    if (nextisIt == isItEnd)
                    {
                        nextisIt = isItBegin;
                    }
                    else
                    {
                        nextisIt = ++tempisIt;

                        if (nextisIt == isItEnd)
                        {
                            nextisIt = ++tempisItBegin;
                        }
                    }

                    break;
                }
                // for YaspGrid
                case GridTypeIndices::yaspGrid:
                {
                    if (nextisIt == isItEnd)
                    {
                        nextisIt = isItBegin;
                    }
                    else
                    {
                        nextisIt = ++tempisIt;

                        if (nextisIt == isItEnd)
                        {
                            nextisIt = ++tempisItBegin;
                        }
                    }

                    break;
                }
                // for UGGrid
                case GridTypeIndices::ugGrid:
                {
                    if (nextisIt == isItEnd)
                    nextisIt = isItBegin;

                    break;
                }
                default:
                {
                    DUNE_THROW(Dune::NotImplemented, "GridType can not be used with MPFAO implementation!");
                    break;
                }
            }

            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer outside = isIt->outside();
                int globalIdxJ = problem_.variables().index(*outside);

                // add off diagonal index
                // add index (row,col) to the matrix
                this->A_.addindex(globalIdxI, globalIdxJ);
            }

            if (isIt->neighbor() && nextisIt->neighbor())
            {
                // access the common neighbor of isIt's and nextisIt's outside
                ElementPointer outside = isIt->outside();
                ElementPointer nextisItoutside = nextisIt->outside();

                IntersectionIterator innerisItEnd = problem_.gridView().iend(*outside);
                IntersectionIterator innernextisItEnd = problem_.gridView().iend(*nextisItoutside);

                for (IntersectionIterator innerisIt = problem_.gridView().ibegin(*outside);
                        innerisIt!=innerisItEnd; ++innerisIt )
                for (IntersectionIterator innernextisIt = problem_.gridView().ibegin(*nextisItoutside);
                        innernextisIt!=innernextisItEnd; ++innernextisIt)
                {
                    if (innerisIt->neighbor() && innernextisIt->neighbor())
                    {
                        ElementPointer innerisItoutside = innerisIt->outside();
                        ElementPointer innernextisItoutside = innernextisIt->outside();

                        if (innerisItoutside == innernextisItoutside && innerisItoutside != isIt->inside())
                        {
                            int globalIdxJ = problem_.variables().index(*innerisItoutside);

                            this->A_.addindex(globalIdxI, globalIdxJ);
                        }
                    }
                }
            }
        } // end of 'for' IntersectionIterator
    } // end of 'for' ElementIterator

    // indicate that all indices are defined, check consistency
    this->A_.endindices();

    return;
}

// only for 2-D general quadrilateral
template<class TypeTag>
void FVMPFAOPressure2P<TypeTag>::assemble()
{
    // introduce matrix R for vector rotation and R is initialized as zero matrix
    FieldMatrix R(0);

    // evaluate matrix R
    if (dim==2)
    for (int i=0; i<dim; ++i)
    {
        R[0][1] = 1;
        R[1][0] = -1;
    }

    // run through all elements
    ElementIterator eItEnd = problem_.gridView().template end<0>();
    for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // get common geometry information for the following computation

        // cell 1 geometry type
        //Dune::GeometryType gt1 = eIt->geometry().type();

        // get global coordinate of cell 1 center
        GlobalPosition globalPos1 = eIt->geometry().center();

        // cell 1 volume
        Scalar volume1 = eIt->geometry().volume();

        // cell 1 index
        int globalIdx1 = problem_.variables().index(*eIt);

        CellData& cellData1 = problem_.variables().cellData(globalIdx1);

        // evaluate right hand side
        PrimaryVariables source(0.0);
        problem_.source(source, *eIt);

        this->f_[globalIdx1] = volume1*(source[wPhaseIdx]/density_[wPhaseIdx] + source[nPhaseIdx]/density_[nPhaseIdx]);

        // get absolute permeability of cell 1
        FieldMatrix K1(problem_.spatialParameters().intrinsicPermeability(*eIt));

        //compute total mobility of cell 1
        Scalar lambda1 = 0;
        lambda1 = cellData1.mobility(wPhaseIdx) + cellData1.mobility(nPhaseIdx);

        // if K1 is zero, no flux through cell1
        // for 2-D
        if (K1[0][0] == 0 && K1[0][1] == 0 && K1[1][0] == 0 && K1[1][1] == 0)
        {
            this->A_[globalIdx1][globalIdx1] += 1.0;
            continue;
        }

        IntersectionIterator isItBegin = problem_.gridView().ibegin(*eIt);
        IntersectionIterator isItEnd = problem_.gridView().iend(*eIt);
        for (IntersectionIterator isIt = isItBegin; isIt!=isItEnd; ++isIt)
        {
            // intersection iterator 'nextisIt' is used to get geometry information
            IntersectionIterator tempisIt = isIt;
            IntersectionIterator tempisItBegin = isItBegin;

            IntersectionIterator nextisIt = ++tempisIt;

            //get nextisIt
            switch (GET_PROP_VALUE(TypeTag, GridImplementation))
            {
                // for SGrid
                case GridTypeIndices::sGrid:
                {
                    if (nextisIt == isItEnd)
                    {
                        nextisIt = isItBegin;
                    }
                    else
                    {
                        nextisIt = ++tempisIt;

                        if (nextisIt == isItEnd)
                        {
                            nextisIt = ++tempisItBegin;
                        }
                    }

                    break;
                }
                // for YaspGrid
                case GridTypeIndices::yaspGrid:
                {
                    if (nextisIt == isItEnd)
                    {
                        nextisIt = isItBegin;
                    }
                    else
                    {
                        nextisIt = ++tempisIt;

                        if (nextisIt == isItEnd)
                        {
                            nextisIt = ++tempisItBegin;
                        }
                    }

                    break;
                }
                // for UGGrid
                case GridTypeIndices::ugGrid:
                {
                    if (nextisIt == isItEnd)
                    nextisIt = isItBegin;

                    break;
                }
                default:
                {
                    DUNE_THROW(Dune::NotImplemented, "GridType can not be used with MPFAO implementation!");
                    break;
                }
            }

            // get geometry type of face 'isIt', i.e., the face between cell1 and cell2 (locally numbered)
            //Dune::GeometryType gtf12 = isIt->geometryInInside().type();

            // center of face in global coordinates, i.e., the midpoint of edge 'isIt'
            GlobalPosition
            globalPosFace12 = isIt->geometry().center();

            // get face volume
            Scalar face12vol = isIt->geometry().volume();

            // get outer normal vector scaled with half volume of face 'isIt'
            Dune::FieldVector<Scalar,dimWorld> integrationOuterNormaln1
            = isIt->centerUnitOuterNormal();
            integrationOuterNormaln1
            *= face12vol/2.0;

            // get geometry type of 'nextisIt', i.e., face between cell1 and cell3 (locally numbered)
            //Dune::GeometryType gtf13 = nextisIt->geometryInInside().type();

            // center of face in global coordinates, i.e., the midpoint of edge 'nextisIt'
            GlobalPosition globalPosFace13
            = nextisIt->geometry().center();

            // get face volume
            Scalar face13vol = nextisIt->geometry().volume();

            // get outer normal vector scaled with half volume of face 'nextisIt'
            Dune::FieldVector<Scalar,dimWorld> integrationOuterNormaln3
            = nextisIt->centerUnitOuterNormal();
            integrationOuterNormaln3
            *= face13vol/2.0;

            // get the intersection node /bar^{x_3} between 'isIt' and 'nextisIt', denoted as 'corner1234'
            // initialization of corner1234
            GlobalPosition corner1234(0);

            // get the global coordinate of corner1234
            for (int i=0; i<isIt->geometry().corners(); ++i)
            {
                GlobalPosition isItcorner = isIt->geometry().corner(i);

                for (int j=0; j<nextisIt->geometry().corners(); ++j)
                {
                    GlobalPosition nextisItcorner = nextisIt->geometry().corner(j);

                    if (nextisItcorner == isItcorner)
                    {
                        corner1234 = isItcorner;
                        continue;
                    }
                }
            }

            // handle interior face
            if (isIt->neighbor())
            {
                // access neighbor cell 2 of 'isIt'
                ElementPointer outside = isIt->outside();
                int globalIdx2 = problem_.variables().index(*outside);

                CellData& cellData2 = problem_.variables().cellData(globalIdx2);

                // neighbor cell 2 geometry type
                //Dune::GeometryType gt2 = outside->geometry().type();

                // get global coordinate of neighbor cell 2 center
                GlobalPosition
                globalPos2 = outside->geometry().center();

                // get absolute permeability of neighbor cell 2
                FieldMatrix K2(problem_.spatialParameters().intrinsicPermeability(*outside));

                // get total mobility of neighbor cell 2
                Scalar lambda2 = 0;
                lambda2 = cellData2.mobility(wPhaseIdx) + cellData2.mobility(nPhaseIdx);

                // 'nextisIt' is an interior face
                if (nextisIt->neighbor())
                {
                    // get basic information of cell 1,2's neighbor cell 3,4
                    // neighbor cell 3
                    // access neighbor cell 3
                    ElementPointer nextisItoutside = nextisIt->outside();
                    int globalIdx3 = problem_.variables().index(*nextisItoutside);

                    CellData& cellData3 = problem_.variables().cellData(globalIdx3);

                    // neighbor cell 3 geometry type
                    //Dune::GeometryType gt3 = nextisItoutside->geometry().type();

                    // get global coordinate of neighbor cell 3 center
                    GlobalPosition
                    globalPos3 = nextisItoutside->geometry().center();

                    // get absolute permeability of neighbor cell 3
                    FieldMatrix K3(problem_.spatialParameters().intrinsicPermeability(*nextisItoutside));

                    // get total mobility of neighbor cell 3
                    Scalar lambda3 = 0;
                    lambda3 = cellData3.mobility(wPhaseIdx) + cellData3.mobility(nPhaseIdx);

                    // neighbor cell 4
                    GlobalPosition globalPos4(0);
                    FieldMatrix K4(0);
                    Scalar lambda4 = 0;
                    int globalIdx4 = 0;

                    IntersectionIterator innerisItEnd = problem_.gridView().iend(*outside);
                    IntersectionIterator innernextisItEnd = problem_.gridView().iend(*nextisItoutside);
                    for (IntersectionIterator innerisIt = problem_.gridView().ibegin(*outside);
                            innerisIt!=innerisItEnd; ++innerisIt )
                    for (IntersectionIterator innernextisIt = problem_.gridView().ibegin(*nextisItoutside);
                            innernextisIt!=innernextisItEnd; ++innernextisIt)
                    {
                        if (innerisIt->neighbor() && innernextisIt->neighbor())
                        {
                            ElementPointer innerisItoutside = innerisIt->outside();
                            ElementPointer innernextisItoutside = innernextisIt->outside();

                            // find the common neighbor cell between cell 2 and cell 3, except cell 1
                            if (innerisItoutside == innernextisItoutside && innerisItoutside != isIt->inside())
                            {
                                // access neighbor cell 4
                                globalIdx4 = problem_.variables().index(*innerisItoutside);

                                CellData& cellData4 = problem_.variables().cellData(globalIdx4);

                                // get global coordinate of neighbor cell 4 center
                                globalPos4 = innerisItoutside->geometry().center();

                                // get absolute permeability of neighbor cell 4
                                K4 += problem_.spatialParameters().intrinsicPermeability(*innerisItoutside);

                                // get total mobility of neighbor cell 4
                                lambda4 = cellData4.mobility(wPhaseIdx) + cellData4.mobility(nPhaseIdx);
                            }
                        }
                    }

                    // computation of flux through the first half edge of 'isIt' and the flux
                    // through the second half edge of 'nextisIt'

                    // get the information of the face 'isIt24' between cell2 and cell4 (locally numbered)
                    IntersectionIterator isIt24 = problem_.gridView().ibegin(*outside);

                    for (IntersectionIterator innerisIt = problem_.gridView().ibegin(*outside);
                            innerisIt != innerisItEnd; ++innerisIt)
                    {
                        if (innerisIt->neighbor())
                        {
                            if (innerisIt->outside() != isIt->inside())
                            {
                                for (int i=0; i<innerisIt->geometry().corners(); ++i)
                                {
                                    GlobalPosition innerisItcorner = innerisIt->geometry().corner(i);

                                    if (innerisItcorner == corner1234)
                                    {
                                        isIt24 = innerisIt;
                                        continue;
                                    }
                                }
                            }
                        }
                    }

                    // get geometry type of face 'isIt24'
                    //Dune::GeometryType gtf24 = isIt24->geometryInInside().type();

                    // center of face in global coordinates, i.e., the midpoint of edge 'isIt24'
                    GlobalPosition
                    globalPosFace24 = isIt24->geometry().center();

                    // get face volume
                    Scalar face24vol = isIt24->geometry().volume();

                    // get outer normal vector scaled with half volume of face 'isIt24'
                    Dune::FieldVector<Scalar,dimWorld> integrationOuterNormaln4
                    = isIt24->centerUnitOuterNormal();
                    integrationOuterNormaln4
                    *= face24vol/2.0;

                    // get the information of the face 'isIt34' between cell3 and cell4 (locally numbered)
                    IntersectionIterator isIt34 = problem_.gridView().ibegin(*nextisItoutside);

                    for (IntersectionIterator innerisIt = problem_.gridView().ibegin(*nextisItoutside);
                            innerisIt != innernextisItEnd; ++innerisIt)
                    {
                        if (innerisIt->neighbor())
                        {
                            if (innerisIt->outside() != isIt->inside())
                            {
                                for (int i=0; i<innerisIt->geometry().corners(); ++i)
                                {
                                    GlobalPosition innerisItcorner = innerisIt->geometry().corner(i);

                                    if (innerisItcorner == corner1234)
                                    {
                                        isIt34 = innerisIt;
                                        continue;
                                    }
                                }
                            }
                        }
                    }

                    // get geometry type of face 'isIt34'
                    //Dune::GeometryType gtf34 = isIt34->geometryInInside().type();

                    // center of face in global coordinates, i.e., the midpoint of edge 'isIt34'
                    GlobalPosition
                    globalPosFace34 = isIt34->geometry().center();

                    // get face volume
                    Scalar face34vol = isIt34->geometry().volume();

                    // get outer normal vector scaled with half volume of face 'isIt34'
                    Dune::FieldVector<Scalar,dimWorld> integrationOuterNormaln2
                    = isIt34->centerUnitOuterNormal();
                    integrationOuterNormaln2
                    *= face34vol/2.0;

                    // compute normal vectors nu11,nu21; nu12, nu22; nu13, nu23; nu14, nu24;
                    FieldVector nu11(0);
                    R.umv(globalPosFace13-globalPos1 ,nu11);

                    FieldVector nu21(0);
                    R.umv(globalPos1-globalPosFace12, nu21);

                    FieldVector nu12(0);
                    R.umv(globalPosFace24-globalPos2, nu12);

                    FieldVector nu22(0);
                    R.umv(globalPosFace12-globalPos2, nu22);

                    FieldVector nu13(0);
                    R.umv(globalPos3-globalPosFace13, nu13);

                    FieldVector nu23(0);
                    R.umv(globalPos3-globalPosFace34, nu23);

                    FieldVector nu14(0);
                    R.umv(globalPos4-globalPosFace24, nu14);

                    FieldVector nu24(0);
                    R.umv(globalPosFace34-globalPos4, nu24);

                    // compute dF1, dF2, dF3, dF4 i.e., the area of quadrilateral made by normal vectors 'nu'
                    FieldVector Rnu21(0);
                    R.umv(nu21, Rnu21);
                    Scalar dF1 = fabs(nu11 * Rnu21);

                    FieldVector Rnu22(0);
                    R.umv(nu22, Rnu22);
                    Scalar dF2 = fabs(nu12 * Rnu22);

                    FieldVector Rnu23(0);
                    R.umv(nu23, Rnu23);
                    Scalar dF3 = fabs(nu13 * Rnu23);

                    FieldVector Rnu24(0);
                    R.umv(nu24, Rnu24);
                    Scalar dF4 = fabs(nu14 * Rnu24);

                    // compute components needed for flux calculation, denoted as 'g'
                    FieldVector K1nu11(0);
                    K1.umv(nu11, K1nu11);
                    FieldVector K1nu21(0);
                    K1.umv(nu21, K1nu21);
                    FieldVector K2nu12(0);
                    K2.umv(nu12, K2nu12);
                    FieldVector K2nu22(0);
                    K2.umv(nu22, K2nu22);
                    FieldVector K3nu13(0);
                    K3.umv(nu13, K3nu13);
                    FieldVector K3nu23(0);
                    K3.umv(nu23, K3nu23);
                    FieldVector K4nu14(0);
                    K4.umv(nu14, K4nu14);
                    FieldVector K4nu24(0);
                    K4.umv(nu24, K4nu24);
                    Scalar g111 = lambda1 * (integrationOuterNormaln1 * K1nu11)/dF1;
                    Scalar g121 = lambda1 * (integrationOuterNormaln1 * K1nu21)/dF1;
                    Scalar g211 = lambda1 * (integrationOuterNormaln3 * K1nu11)/dF1;
                    Scalar g221 = lambda1 * (integrationOuterNormaln3 * K1nu21)/dF1;
                    Scalar g112 = lambda2 * (integrationOuterNormaln1 * K2nu12)/dF2;
                    Scalar g122 = lambda2 * (integrationOuterNormaln1 * K2nu22)/dF2;
                    Scalar g212 = lambda2 * (integrationOuterNormaln4 * K2nu12)/dF2;
                    Scalar g222 = lambda2 * (integrationOuterNormaln4 * K2nu22)/dF2;
                    Scalar g113 = lambda3 * (integrationOuterNormaln2 * K3nu13)/dF3;
                    Scalar g123 = lambda3 * (integrationOuterNormaln2 * K3nu23)/dF3;
                    Scalar g213 = lambda3 * (integrationOuterNormaln3 * K3nu13)/dF3;
                    Scalar g223 = lambda3 * (integrationOuterNormaln3 * K3nu23)/dF3;
                    Scalar g114 = lambda4 * (integrationOuterNormaln2 * K4nu14)/dF4;
                    Scalar g124 = lambda4 * (integrationOuterNormaln2 * K4nu24)/dF4;
                    Scalar g214 = lambda4 * (integrationOuterNormaln4 * K4nu14)/dF4;
                    Scalar g224 = lambda4 * (integrationOuterNormaln4 * K4nu24)/dF4;

                    // compute transmissibility matrix T = CA^{-1}B+F
                    Dune::FieldMatrix<Scalar,2*dim,2*dim> C(0), F(0), A(0), B(0);

                    // evaluate matrix C, F, A, B
                    C[0][0] = -g111;
                    C[0][2] = -g121;
                    C[1][1] = g114;
                    C[1][3] = g124;
                    C[2][1] = -g213;
                    C[2][2] = g223;
                    C[3][0] = g212;
                    C[3][3] = -g222;

                    F[0][0] = g111 + g121;
                    F[1][3] = -g114 - g124;
                    F[2][2] = g213 - g223;
                    F[3][1] = -g212 + g222;

                    A[0][0] = g111 + g112;
                    A[0][2] = g121;
                    A[0][3] = -g122;
                    A[1][1] = g114 + g113;
                    A[1][2] = -g123;
                    A[1][3] = g124;
                    A[2][0] = g211;
                    A[2][1] = -g213;
                    A[2][2] = g223 + g221;
                    A[3][0] = -g212;
                    A[3][1] = g214;
                    A[3][3] = g222 + g224;

                    B[0][0] = g111 + g121;
                    B[0][1] = g112 - g122;
                    B[1][2] = g113 - g123;
                    B[1][3] = g114 + g124;
                    B[2][0] = g211 + g221;
                    B[2][2] = -g213 + g223;
                    B[3][1] = -g212 + g222;
                    B[3][3] = g214 + g224;

                    // compute T
                    A.invert();
                    F += B.leftmultiply(C.rightmultiply(A));
                    Dune::FieldMatrix<Scalar,2*dim,2*dim> T(F);

                    // assemble the global matrix this->A_ and right hand side f
                    this->A_[globalIdx1][globalIdx1] += T[0][0] + T[2][0];
                    this->A_[globalIdx1][globalIdx2] += T[0][1] + T[2][1];
                    this->A_[globalIdx1][globalIdx3] += T[0][2] + T[2][2];
                    this->A_[globalIdx1][globalIdx4] += T[0][3] + T[2][3];

                }
                // 'nextisIt' is on the boundary

                else
                {
                    // computation of flux through the first half edge of 'isIt' and the flux
                    // through the second half edge of 'nextisIt'

                    // get common geometry information for the following computation
                    // get the information of the face 'isIt24' between cell2 and cell4 (locally numbered)
                    IntersectionIterator isIt24 = problem_.gridView().ibegin(*outside);
                    IntersectionIterator innerisItEnd = problem_.gridView().iend(*outside);
                    for (IntersectionIterator innerisIt = problem_.gridView().ibegin(*outside);
                            innerisIt != innerisItEnd; ++innerisIt)
                    {
                        if (innerisIt->boundary())
                        {
                            for (int i=0; i<innerisIt->geometry().corners(); ++i)
                            {
                                GlobalPosition innerisItcorner = innerisIt->geometry().corner(i);

                                if (innerisItcorner == corner1234)
                                {
                                    isIt24 = innerisIt;
                                    continue;
                                }
                            }
                        }
                    }

                    // get geometry type of face 'isIt24'
                    //Dune::GeometryType gtf24 = isIt24->geometryInInside().type();

                    // center of face in global coordinates, i.e., the midpoint of edge 'isIt24'
                    GlobalPosition
                    globalPosFace24 = isIt24->geometry().center();

                    // get face volume
                    Scalar face24vol = isIt24->geometry().volume();

                    // get outer normal vector scaled with half volume of face 'isIt24'
                    Dune::FieldVector<Scalar,dimWorld> integrationOuterNormaln4
                    = isIt24->centerUnitOuterNormal();
                    integrationOuterNormaln4
                    *= face24vol/2.0;

                    BoundaryTypes nextIsItBcType;
                    problem_.boundaryTypes(nextIsItBcType, *nextisIt);

                    // get boundary condition for boundary face (isIt24) center
                    BoundaryTypes isIt24BcType;
                    problem_.boundaryTypes(isIt24BcType, *isIt24);

                    PrimaryVariables boundValues(0.0);

                    // 'nextisIt': Neumann boundary
                    if (nextIsItBcType.isNeumann(pressEqIdx))
                    {
                        // get Neumann boundary value of 'nextisIt'
                        problem_.neumann(boundValues, *nextisIt);
                        Scalar J3 = (boundValues[wPhaseIdx]/density_[wPhaseIdx]+boundValues[nPhaseIdx]/density_[nPhaseIdx]);

                        // 'isIt24': Neumann boundary
                        if (isIt24BcType.isNeumann(pressEqIdx))
                        {
                            // get neumann boundary value of 'isIt24'
                            problem_.neumann(boundValues, *isIt24);
                            Scalar J4 = (boundValues[wPhaseIdx]/density_[wPhaseIdx]+boundValues[nPhaseIdx]/density_[nPhaseIdx]);

                            // compute normal vectors nu11,nu21; nu12, nu22;
                            FieldVector nu11(0);
                            R.umv(globalPosFace13-globalPos1 ,nu11);

                            FieldVector nu21(0);
                            R.umv(globalPos1-globalPosFace12, nu21);

                            FieldVector nu12(0);
                            R.umv(globalPosFace24-globalPos2, nu12);

                            FieldVector nu22(0);
                            R.umv(globalPosFace12-globalPos2, nu22);

                            // compute dF1, dF2 i.e., the area of quadrilateral made by normal vectors 'nu'
                            FieldVector Rnu21(0);
                            R.umv(nu21, Rnu21);
                            Scalar dF1 = fabs(nu11 * Rnu21);

                            FieldVector Rnu22(0);
                            R.umv(nu22, Rnu22);
                            Scalar dF2 = fabs(nu12 * Rnu22);

                            // compute components needed for flux calculation, denoted as 'g'
                            FieldVector K1nu11(0);
                            K1.umv(nu11, K1nu11);
                            FieldVector K1nu21(0);
                            K1.umv(nu21, K1nu21);
                            FieldVector K2nu12(0);
                            K2.umv(nu12, K2nu12);
                            FieldVector K2nu22(0);
                            K2.umv(nu22, K2nu22);
                            Scalar g111 = lambda1 * (integrationOuterNormaln1 * K1nu11)/dF1;
                            Scalar g121 = lambda1 * (integrationOuterNormaln1 * K1nu21)/dF1;
                            Scalar g211 = lambda1 * (integrationOuterNormaln3 * K1nu11)/dF1;
                            Scalar g221 = lambda1 * (integrationOuterNormaln3 * K1nu21)/dF1;
                            Scalar g112 = lambda2 * (integrationOuterNormaln1 * K2nu12)/dF2;
                            Scalar g122 = lambda2 * (integrationOuterNormaln1 * K2nu22)/dF2;
                            Scalar g212 = lambda2 * (integrationOuterNormaln4 * K2nu12)/dF2;
                            Scalar g222 = lambda2 * (integrationOuterNormaln4 * K2nu22)/dF2;

                            // compute the matrix T & vector r in v = A^{-1}(Bu + r1) = Tu + r
                            Dune::FieldMatrix<Scalar,2*dim-1,2*dim-1> A(0);
                            Dune::FieldMatrix<Scalar,2*dim-1,dim> B(0);
                            Dune::FieldVector<Scalar,2*dim-1> r1(0), r(0);

                            // evaluate matrix A, B
                            A[0][0] = g111 + g112;
                            A[0][1] = g121;
                            A[0][2] = -g122;
                            A[1][0] = g211;
                            A[1][1] = g221;
                            A[2][0] = -g212;
                            A[2][2] = g222;

                            B[0][0] = g111 + g121;
                            B[0][1] = g112 - g122;
                            B[1][0] = g211 + g221;
                            B[2][1] = g222 - g212;

                            // evaluate vector r1
                            r1[1] = -J3 * nextisIt->geometry().volume()/2.0;
                            r1[2] = -J4 * isIt24->geometry().volume()/2.0;

                            // compute T and r
                            A.invert();
                            B.leftmultiply(A);
                            Dune::FieldMatrix<Scalar,2*dim-1,dim> T(B);
                            A.umv(r1, r);

                            // assemble the global matrix this->A_ and right hand side f
                            this->A_[globalIdx1][globalIdx1] += g111 + g121 - g111 * T[0][0] - g121 * T[1][0];
                            this->A_[globalIdx1][globalIdx2] += -g111 * T[0][1] - g121 * T[1][1];
                            this->f_[globalIdx1] += g111 * r[0] + g121 * r[1];

                        }
                        // 'isIt24': Dirichlet boundary

                        else if (isIt24BcType.isDirichlet(pressEqIdx))
                        {
                            // get Dirichlet boundary value on 'isIt24'
                            problem_.dirichlet(boundValues, *isIt24);
                            Scalar g4 = boundValues[pressureIdx];

                            // compute total mobility for Dirichlet boundary 'isIt24'
                            //determine lambda at the boundary -> if no saturation is known directly at the boundary use the cell saturation
                            Scalar alambda2 = 0;
                            if (isIt24BcType.isDirichlet(satEqIdx))
                            {
                                Scalar satBound = boundValues[saturationIdx];

                                //determine phase saturations from primary saturation variable
                                Scalar satW = 0;
                                switch (saturationType_)
                                {
                                case Sw:
                                {
                                    satW = satBound;
                                    break;
                                }
                                case Sn:
                                {
                                    satW = 1 - satBound;
                                    break;
                                }
                                }

                                Scalar lambdaWBound = 0;
                                Scalar lambdaNWBound = 0;

                                lambdaWBound = MaterialLaw::krw(
                                        problem_.spatialParameters().materialLawParams(*eIt), satW)
                                        / viscosity_[wPhaseIdx];
                                lambdaNWBound = MaterialLaw::krn(
                                        problem_.spatialParameters().materialLawParams(*eIt), satW)
                                        / viscosity_[nPhaseIdx];
                                alambda2 = lambdaWBound + lambdaNWBound;
                            }
                            else
                            {
                                alambda2 = lambda2;
                            }

                            // compute normal vectors nu11,nu21; nu12, nu22;
                            FieldVector nu11(0);
                            R.umv(globalPosFace13-globalPos1 ,nu11);

                            FieldVector nu21(0);
                            R.umv(globalPos1-globalPosFace12, nu21);

                            FieldVector nu12(0);
                            R.umv(globalPosFace24-globalPos2, nu12);

                            FieldVector nu22(0);
                            R.umv(globalPosFace12-globalPos2, nu22);

                            // compute dF1, dF2 i.e., the area of quadrilateral made by normal vectors 'nu'
                            FieldVector Rnu21(0);
                            R.umv(nu21, Rnu21);
                            Scalar dF1 = fabs(nu11 * Rnu21);

                            FieldVector Rnu22(0);
                            R.umv(nu22, Rnu22);
                            Scalar dF2 = fabs(nu12 * Rnu22);

                            // compute components needed for flux calculation, denoted as 'g'
                            FieldVector K1nu11(0);
                            K1.umv(nu11, K1nu11);
                            FieldVector K1nu21(0);
                            K1.umv(nu21, K1nu21);
                            FieldVector K2nu12(0);
                            K2.umv(nu12, K2nu12);
                            FieldVector K2nu22(0);
                            K2.umv(nu22, K2nu22);
                            Scalar g111 = lambda1 * (integrationOuterNormaln1 * K1nu11)/dF1;
                            Scalar g121 = lambda1 * (integrationOuterNormaln1 * K1nu21)/dF1;
                            Scalar g211 = lambda1 * (integrationOuterNormaln3 * K1nu11)/dF1;
                            Scalar g221 = lambda1 * (integrationOuterNormaln3 * K1nu21)/dF1;
                            Scalar g112 = alambda2 * (integrationOuterNormaln1 * K2nu12)/dF2;
                            Scalar g122 = alambda2 * (integrationOuterNormaln1 * K2nu22)/dF2;

                            // compute the matrix T & vector r in v = A^{-1}(Bu + r1) = Tu + r
                            FieldMatrix A(0), B(0);
                            FieldVector r1(0), r(0);

                            // evaluate matrix A, B
                            A[0][0] = g111 + g112;
                            A[0][1] = g121;
                            A[1][0] = g211;
                            A[1][1] = g221;

                            B[0][0] = g111 + g121;
                            B[0][1] = g112 - g122;
                            B[1][0] = g211 + g221;

                            // evaluate vector r1
                            r1[0] = g122 * g4;
                            r1[1] = -J3 * nextisIt->geometry().volume()/2.0;

                            // compute T and r
                            A.invert();
                            B.leftmultiply(A);
                            FieldMatrix T(B);
                            A.umv(r1, r);

                            // assemble the global matrix this->A_ and right hand side f
                            this->A_[globalIdx1][globalIdx1] += g111 + g121 - g111 * T[0][0] - g121 * T[1][0];
                            this->A_[globalIdx1][globalIdx2] += -g111 * T[0][1] - g121 * T[1][1];
                            this->f_[globalIdx1] += g111 * r[0] + g121 * r[1];

                        }
                    }
                    // 'nextisIt': Dirichlet boundary

                    else if (nextIsItBcType.isDirichlet(pressEqIdx))
                    {
                        // get Dirichlet boundary value of 'nextisIt'
                        problem_.dirichlet(boundValues, *nextisIt);
                        Scalar g3 = boundValues[pressureIdx];

                        // compute total mobility for Dirichlet boundary 'nextisIt'
                        //determine lambda at the boundary -> if no saturation is known directly at the boundary use the cell saturation
                        Scalar alambda1 = 0;
                        if (nextIsItBcType.isDirichlet(satEqIdx))
                        {
                            Scalar satBound = boundValues[saturationIdx];

                            //determine phase saturations from primary saturation variable
                            Scalar satW = 0;
                            switch (saturationType_)
                            {
                            case Sw:
                            {
                                satW = satBound;
                                break;
                            }
                            case Sn:
                            {
                                satW = 1 - satBound;
                                break;
                            }
                            }

                            Scalar lambdaWBound = 0;
                            Scalar lambdaNWBound = 0;

                            lambdaWBound = MaterialLaw::krw(
                                    problem_.spatialParameters().materialLawParams(*eIt), satW)
                                    / viscosity_[wPhaseIdx];
                            lambdaNWBound = MaterialLaw::krn(
                                    problem_.spatialParameters().materialLawParams(*eIt), satW)
                                    / viscosity_[nPhaseIdx];
                            alambda1 = lambdaWBound + lambdaNWBound;
                        }
                        else
                        {
                            alambda1 = lambda1;
                        }

                        // 'isIt24': Neumann boundary
                        if (isIt24BcType.isNeumann(pressEqIdx))
                        {
                            // get Neumann boundary value of 'isIt24'
                            problem_.neumann(boundValues, *isIt24);
                            Scalar J4 = (boundValues[wPhaseIdx]/density_[wPhaseIdx]+boundValues[nPhaseIdx]/density_[nPhaseIdx]);

                            // compute normal vectors nu11,nu21; nu12, nu22;
                            FieldVector nu11(0);
                            R.umv(globalPosFace13-globalPos1 ,nu11);

                            FieldVector nu21(0);
                            R.umv(globalPos1-globalPosFace12, nu21);

                            FieldVector nu12(0);
                            R.umv(globalPosFace24-globalPos2, nu12);

                            FieldVector nu22(0);
                            R.umv(globalPosFace12-globalPos2, nu22);

                            // compute dF1, dF2 i.e., the area of quadrilateral made by normal vectors 'nu'
                            FieldVector Rnu21(0);
                            R.umv(nu21, Rnu21);
                            Scalar dF1 = fabs(nu11 * Rnu21);

                            FieldVector Rnu22(0);
                            R.umv(nu22, Rnu22);
                            Scalar dF2 = fabs(nu12 * Rnu22);

                            // compute components needed for flux calculation, denoted as 'g'
                            FieldVector K1nu11(0);
                            K1.umv(nu11, K1nu11);
                            FieldVector K1nu21(0);
                            K1.umv(nu21, K1nu21);
                            FieldVector K2nu12(0);
                            K2.umv(nu12, K2nu12);
                            FieldVector K2nu22(0);
                            K2.umv(nu22, K2nu22);
                            Scalar g111 = alambda1 * (integrationOuterNormaln1 * K1nu11)/dF1;
                            Scalar g121 = alambda1 * (integrationOuterNormaln1 * K1nu21)/dF1;
                            Scalar g211 = alambda1 * (integrationOuterNormaln3 * K1nu11)/dF1;
                            Scalar g221 = alambda1 * (integrationOuterNormaln3 * K1nu21)/dF1;
                            Scalar g112 = lambda2 * (integrationOuterNormaln1 * K2nu12)/dF2;
                            Scalar g122 = lambda2 * (integrationOuterNormaln1 * K2nu22)/dF2;
                            Scalar g212 = lambda2 * (integrationOuterNormaln4 * K2nu12)/dF2;
                            Scalar g222 = lambda2 * (integrationOuterNormaln4 * K2nu22)/dF2;

                            // compute the matrix T & vector r in v = A^{-1}(Bu + r1) = Tu + r
                            FieldMatrix A(0), B(0);
                            FieldVector r1(0), r(0);

                            // evaluate matrix A, B
                            A[0][0] = g111 + g112;
                            A[0][1] = -g122;
                            A[1][0] = -g212;
                            A[1][1] = g222;

                            B[0][0] = g111 + g121;
                            B[0][1] = g112 - g122;
                            B[1][1] = g222 - g212;

                            // evaluate vector r1
                            r1[0] = -g121 * g3;
                            r1[1] = -J4 * isIt24->geometry().volume()/2.0;

                            // compute T and r
                            A.invert();
                            B.leftmultiply(A);
                            FieldMatrix T(B);
                            A.umv(r1, r);

                            // assemble the global matrix this->A_ and right hand side this->f_
                            this->A_[globalIdx1][globalIdx1] += (g111 + g121 - g111 * T[0][0]) + (g211 + g221 - g211 * T[0][0]);
                            this->A_[globalIdx1][globalIdx2] += -g111 * T[0][1] - g211 * T[0][1];
                            this->f_[globalIdx1] += (g121 + g221) * g3 + (g111 + g211) * r[0];

                        }
                        // 'isIt24': Dirichlet boundary

                        else if (isIt24BcType.isDirichlet(pressEqIdx))
                        {
                            // get Dirichlet boundary value on 'isIt24'
                            problem_.dirichlet(boundValues, *isIt24);
                            Scalar g4 = boundValues[pressureIdx];

                            // compute total mobility for Dirichlet boundary 'isIt24'
                            //determine lambda at the boundary -> if no saturation is known directly at the boundary use the cell saturation
                            Scalar alambda2 = 0;
                            if (isIt24BcType.isDirichlet(satEqIdx))
                            {
                                Scalar satBound = boundValues[saturationIdx];

                                //determine phase saturations from primary saturation variable
                                Scalar satW = 0;
                                switch (saturationType_)
                                {
                                case Sw:
                                {
                                    satW = satBound;
                                    break;
                                }
                                case Sn:
                                {
                                    satW = 1 - satBound;
                                    break;
                                }
                                }

                                Scalar lambdaWBound = 0;
                                Scalar lambdaNWBound = 0;

                                lambdaWBound = MaterialLaw::krw(
                                        problem_.spatialParameters().materialLawParams(*eIt), satW)
                                        / viscosity_[wPhaseIdx];
                                lambdaNWBound = MaterialLaw::krn(
                                        problem_.spatialParameters().materialLawParams(*eIt), satW)
                                        / viscosity_[nPhaseIdx];
                                alambda2 = lambdaWBound + lambdaNWBound;
                            }
                            else
                            {
                                alambda2 = lambda2;
                            }

                            // compute normal vectors nu11,nu21; nu12, nu22;
                            FieldVector nu11(0);
                            R.umv(globalPosFace13-globalPos1 ,nu11);

                            FieldVector nu21(0);
                            R.umv(globalPos1-globalPosFace12, nu21);

                            FieldVector nu12(0);
                            R.umv(globalPosFace24-globalPos2, nu12);

                            FieldVector nu22(0);
                            R.umv(globalPosFace12-globalPos2, nu22);

                            // compute dF1, dF2 i.e., the area of quadrilateral made by normal vectors 'nu'
                            FieldVector Rnu21(0);
                            R.umv(nu21, Rnu21);
                            Scalar dF1 = fabs(nu11 * Rnu21);

                            FieldVector Rnu22(0);
                            R.umv(nu22, Rnu22);
                            Scalar dF2 = fabs(nu12 * Rnu22);

                            // compute components needed for flux calculation, denoted as 'g'
                            FieldVector K1nu11(0);
                            K1.umv(nu11, K1nu11);
                            FieldVector K1nu21(0);
                            K1.umv(nu21, K1nu21);
                            FieldVector K2nu12(0);
                            K2.umv(nu12, K2nu12);
                            FieldVector K2nu22(0);
                            K2.umv(nu22, K2nu22);
                            Scalar g111 = alambda1 * (integrationOuterNormaln1 * K1nu11)/dF1;
                            Scalar g121 = alambda1 * (integrationOuterNormaln1 * K1nu21)/dF1;
                            Scalar g211 = alambda1 * (integrationOuterNormaln3 * K1nu11)/dF1;
                            Scalar g221 = alambda1 * (integrationOuterNormaln3 * K1nu21)/dF1;
                            Scalar g112 = alambda2 * (integrationOuterNormaln1 * K2nu12)/dF2;
                            Scalar g122 = alambda2 * (integrationOuterNormaln1 * K2nu22)/dF2;

                            // compute the matrix T & vector r
                            FieldMatrix T(0);
                            FieldVector r(0);

                            Scalar coe = g111 + g112;

                            // evaluate matrix T
                            T[0][0] = g112 * (g111 + g121)/coe;
                            T[0][1] = -g111 * (g112 - g122)/coe;
                            T[1][0] = g221 + g211 * (g112 - g121)/coe;
                            T[1][1] = -g211 * (g112 - g122)/coe;

                            // evaluate vector r
                            r[0] = -(g4 * g122 * g111 + g3 * g112 * g121)/coe;
                            r[1] = -g221 * g3 + (g3 * g211 * g121 - g4 * g211 * g122)/coe;
                            // assemble the global matrix this->A_ and right hand side this->f_
                            this->A_[globalIdx1][globalIdx1] += T[0][0] + T[1][0];
                            this->A_[globalIdx1][globalIdx2] += T[0][1] + T[1][1];
                            this->f_[globalIdx1] -= r[0] + r[1];

                        }
                    }
                }
            }
            // handle boundary face 'isIt'

            else
            {
                BoundaryTypes isItBcType;
                problem_.boundaryTypes(isItBcType, *isIt);

                PrimaryVariables boundValues(0.0);

                // 'isIt' is on Neumann boundary
                if (isItBcType.isNeumann(pressEqIdx))
                {
                    // get Neumann boundary value
                    problem_.neumann(boundValues, *isIt);
                    Scalar J1 = (boundValues[wPhaseIdx]/density_[wPhaseIdx]+boundValues[nPhaseIdx]/density_[nPhaseIdx]);

                    // evaluate right hand side
                    this->f_[globalIdx1] -= face12vol*J1;

                    // 'nextisIt' is on boundary
                    if (nextisIt->boundary())
                    {
                        // get boundary condition for boundary face center of 'nextisIt'
                        BoundaryTypes nextIsItBcType;
                        problem_.boundaryTypes(nextIsItBcType, *nextisIt);

                        if (nextIsItBcType.isDirichlet(pressEqIdx))
                        {
                            problem_.dirichlet(boundValues, *nextisIt);

                            // compute total mobility for Dirichlet boundary 'nextisIt'
                            //determine lambda at the boundary -> if no saturation is known directly at the boundary use the cell saturation
                            Scalar alambda1 = 0;
                            if (nextIsItBcType.isDirichlet(satEqIdx))
                            {
                                Scalar satBound = boundValues[saturationIdx];

                                //determine phase saturations from primary saturation variable
                                Scalar satW = 0;
                                switch (saturationType_)
                                {
                                case Sw:
                                {
                                    satW = satBound;
                                    break;
                                }
                                case Sn:
                                {
                                    satW = 1 - satBound;
                                    break;
                                }
                                }

                                Scalar lambdaWBound = 0;
                                Scalar lambdaNWBound = 0;

                                 lambdaWBound = MaterialLaw::krw(
                                        problem_.spatialParameters().materialLawParams(*eIt), satW)
                                        / viscosity_[wPhaseIdx];
                                lambdaNWBound = MaterialLaw::krn(
                                        problem_.spatialParameters().materialLawParams(*eIt), satW)
                                        / viscosity_[nPhaseIdx];
                                alambda1 = lambdaWBound + lambdaNWBound;
                            }
                            else
                            {
                                alambda1 = lambda1;
                            }

                            // get Dirichlet boundary value
                            Scalar g3 = boundValues[pressureIdx];

                            // compute normal vectors nu11,nu21;
                            FieldVector nu11(0);
                            R.umv(globalPosFace13-globalPos1 ,nu11);

                            FieldVector nu21(0);
                            R.umv(globalPos1-globalPosFace12, nu21);

                            // compute dF1, dF2 i.e., the area of quadrilateral made by normal vectors 'nu'
                            FieldVector Rnu21(0);
                            R.umv(nu21, Rnu21);
                            Scalar dF1 = fabs(nu11 * Rnu21);

                            // compute components needed for flux calculation, denoted as 'g'
                            FieldVector K1nu11(0);
                            K1.umv(nu11, K1nu11);
                            FieldVector K1nu21(0);
                            K1.umv(nu21, K1nu21);
                            Scalar g111 = alambda1 * (integrationOuterNormaln1 * K1nu11)/dF1;
                            Scalar g121 = alambda1 * (integrationOuterNormaln1 * K1nu21)/dF1;
                            Scalar g211 = alambda1 * (integrationOuterNormaln3 * K1nu11)/dF1;
                            Scalar g221 = alambda1 * (integrationOuterNormaln3 * K1nu21)/dF1;

                            // assemble the global matrix this->A_ and right hand side this->f_
                            this->A_[globalIdx1][globalIdx1] += g221 - g211 * g121/g111;
                            this->f_[globalIdx1] -= (g211 * g121/g111 - g221) * g3 - (g211 * (-J1) * face12vol)/(2.0 * g111);

                        }
                    }
                    // 'nextisIt' is inside

                    else
                    {
                        // neighbor cell 3
                        // access neighbor cell 3
                        ElementPointer nextisItoutside = nextisIt->outside();
                        int globalIdx3 = problem_.variables().index(*nextisItoutside);

                        CellData& cellData3 = problem_.variables().cellData(globalIdx3);

                        // neighbor cell 3 geometry type
                        //Dune::GeometryType gt3 = nextisItoutside->geometry().type();

                        // get global coordinate of neighbor cell 3 center
                        GlobalPosition
                        globalPos3 = nextisItoutside->geometry().center();

                        // get absolute permeability of neighbor cell 3
                        FieldMatrix K3(problem_.spatialParameters().intrinsicPermeability(*nextisItoutside));

                        // get total mobility of neighbor cell 3
                        Scalar lambda3 = 0;
                        lambda3 = cellData3.mobility(wPhaseIdx) + cellData3.mobility(nPhaseIdx);

                        // get the information of the face 'isIt34' between cell3 and cell4 (locally numbered)
                        IntersectionIterator isIt34 = problem_.gridView().ibegin(*nextisItoutside);
                        IntersectionIterator innernextisItEnd = problem_.gridView().iend(*nextisItoutside);
                        for (IntersectionIterator innerisIt = problem_.gridView().ibegin(*nextisItoutside);
                                innerisIt != innernextisItEnd; ++innerisIt)
                        {
                            if (innerisIt->boundary())
                            {
                                for (int i=0; i<innerisIt->geometry().corners(); ++i)
                                {
                                    GlobalPosition innerisItcorner = innerisIt->geometry().corner(i);

                                    if (innerisItcorner == corner1234)
                                    {
                                        isIt34 = innerisIt;
                                        continue;
                                    }
                                }
                            }
                        }

                        // get geometry type of face 'isIt34'
                        //Dune::GeometryType gtf34 = isIt34->geometryInInside().type();

                        // center of face in global coordinates, i.e., the midpoint of edge 'isIt34'
                        GlobalPosition
                        globalPosFace34 = isIt34->geometry().center();

                        // get face volume
                        Scalar face34vol = isIt34->geometry().volume();

                        // get outer normal vector scaled with half volume of face 'isIt34'
                        Dune::FieldVector<Scalar,dimWorld> integrationOuterNormaln2
                        = isIt34->centerUnitOuterNormal();
                        integrationOuterNormaln2
                        *= face34vol/2.0;

                        // get boundary condition for boundary face center of 'isIt34'
                        BoundaryTypes isIt34BcType;
                        problem_.boundaryTypes(isIt34BcType, *isIt34);

                        // 'isIt34': Neumann boundary
                        if (isIt34BcType.isNeumann(pressEqIdx))
                        {
                            // get Neumann boundary value
                            problem_.neumann(boundValues, *isIt34);
                            Scalar J2 = (boundValues[wPhaseIdx]/density_[wPhaseIdx]+boundValues[nPhaseIdx]/density_[nPhaseIdx]);

                            // compute normal vectors nu11,nu21; nu13, nu23;
                            FieldVector nu11(0);
                            R.umv(globalPosFace13-globalPos1 ,nu11);

                            FieldVector nu21(0);
                            R.umv(globalPos1-globalPosFace12, nu21);

                            FieldVector nu13(0);
                            R.umv(globalPos3-globalPosFace13, nu13);

                            FieldVector nu23(0);
                            R.umv(globalPos3-globalPosFace34, nu23);

                            // compute dF1, dF3 i.e., the area of quadrilateral made by normal vectors 'nu'
                            FieldVector Rnu21(0);
                            R.umv(nu21, Rnu21);
                            Scalar dF1 = fabs(nu11 * Rnu21);

                            FieldVector Rnu23(0);
                            R.umv(nu23, Rnu23);
                            Scalar dF3 = fabs(nu13 * Rnu23);

                            // compute components needed for flux calculation, denoted as 'g'
                            FieldVector K1nu11(0);
                            K1.umv(nu11, K1nu11);
                            FieldVector K1nu21(0);
                            K1.umv(nu21, K1nu21);
                            FieldVector K3nu13(0);
                            K3.umv(nu13, K3nu13);
                            FieldVector K3nu23(0);
                            K3.umv(nu23, K3nu23);
                            Scalar g111 = lambda1 * (integrationOuterNormaln1 * K1nu11)/dF1;
                            Scalar g121 = lambda1 * (integrationOuterNormaln1 * K1nu21)/dF1;
                            Scalar g211 = lambda1 * (integrationOuterNormaln3 * K1nu11)/dF1;
                            Scalar g221 = lambda1 * (integrationOuterNormaln3 * K1nu21)/dF1;
                            Scalar g113 = lambda3 * (integrationOuterNormaln2 * K3nu13)/dF3;
                            Scalar g123 = lambda3 * (integrationOuterNormaln2 * K3nu23)/dF3;
                            Scalar g213 = lambda3 * (integrationOuterNormaln3 * K3nu13)/dF3;
                            Scalar g223 = lambda3 * (integrationOuterNormaln3 * K3nu23)/dF3;

                            // compute transmissibility matrix T = CA^{-1}B+F
                            Dune::FieldMatrix<Scalar,2*dim-1,2*dim-1> C(0), A(0);
                            Dune::FieldMatrix<Scalar,2*dim-1,dim> F(0), B(0);

                            // evaluate matrix C, F, A, B
                            C[0][0] = -g111;
                            C[0][2] = -g121;
                            C[1][1] = -g113;
                            C[1][2] = g123;
                            C[2][1] = -g213;
                            C[2][2] = g223;

                            F[0][0] = g111 + g121;
                            F[1][1] = g113 - g123;
                            F[2][1] = g213 - g223;

                            A[0][0] = g111;
                            A[0][2] = g121;
                            A[1][1] = g113;
                            A[1][2] = -g123;
                            A[2][0] = g211;
                            A[2][1] = -g213;
                            A[2][2] = g223 + g221;

                            B[0][0] = g111 + g121;
                            B[1][1] = g113 - g123;
                            B[2][0] = g211 + g221;
                            B[2][1] = g223 - g213;

                            // compute T
                            A.invert();
                            Dune::FieldMatrix<Scalar,2*dim-1,2*dim-1> CAinv(C.rightmultiply(A));
                            F += B.leftmultiply(CAinv);
                            Dune::FieldMatrix<Scalar,2*dim-1,dim> T(F);

                            // compute vector r
                            // evaluate r1
                            Dune::FieldVector<Scalar,2*dim-1> r1(0);
                            r1[0] = -J1 * face12vol/2.0;
                            r1[1] = -J2 * isIt34->geometry().volume()/2.0;

                            // compute r = CA^{-1}r1
                            Dune::FieldVector<Scalar,2*dim-1> r(0);
                            CAinv.umv(r1, r);

                            // assemble the global matrix this->A_ and right hand side this->f_
                            this->A_[globalIdx1][globalIdx1] += T[2][0];
                            this->A_[globalIdx1][globalIdx3] += T[2][1];
                            this->f_[globalIdx1] -= r[2];

                        }
                        // 'isIt34': Dirichlet boundary

                        else if (isIt34BcType.isDirichlet(pressEqIdx))
                        {
                            // get Dirichlet boundary value
                            problem_.dirichlet(boundValues, *isIt34);
                            Scalar g2 = boundValues[pressureIdx];

                            // compute total mobility for Dirichlet boundary 'isIt24'
                            //determine lambda at the boundary -> if no saturation is known directly at the boundary use the cell saturation
                            Scalar alambda3 = 0;
                            if (isIt34BcType.isDirichlet(satEqIdx))
                            {
                                Scalar satBound = boundValues[saturationIdx];

                                //determine phase saturations from primary saturation variable
                                Scalar satW = 0;
                                switch (saturationType_)
                                {
                                case Sw:
                                {
                                    satW = satBound;
                                    break;
                                }
                                case Sn:
                                {
                                    satW = 1 - satBound;
                                    break;
                                }
                                }

                                Scalar lambdaWBound = 0;
                                Scalar lambdaNWBound = 0;

                                lambdaWBound = MaterialLaw::krw(
                                        problem_.spatialParameters().materialLawParams(*eIt), satW)
                                        / viscosity_[wPhaseIdx];
                                lambdaNWBound = MaterialLaw::krn(
                                        problem_.spatialParameters().materialLawParams(*eIt), satW)
                                        / viscosity_[nPhaseIdx];
                                alambda3 = lambdaWBound + lambdaNWBound;
                            }
                            else
                            {
                                alambda3 = lambda3;
                            }

                            // compute normal vectors nu11,nu21; nu13, nu23;
                            FieldVector nu11(0);
                            R.umv(globalPosFace13-globalPos1 ,nu11);

                            FieldVector nu21(0);
                            R.umv(globalPos1-globalPosFace12, nu21);

                            FieldVector nu13(0);
                            R.umv(globalPos3-globalPosFace13, nu13);

                            FieldVector nu23(0);
                            R.umv(globalPos3-globalPosFace34, nu23);

                            // compute dF1, dF3 i.e., the area of quadrilateral made by normal vectors 'nu'
                            FieldVector Rnu21(0);
                            R.umv(nu21, Rnu21);
                            Scalar dF1 = fabs(nu11 * Rnu21);

                            FieldVector Rnu23(0);
                            R.umv(nu23, Rnu23);
                            Scalar dF3 = fabs(nu13 * Rnu23);

                            // compute components needed for flux calculation, denoted as 'g'
                            FieldVector K1nu11(0);
                            K1.umv(nu11, K1nu11);
                            FieldVector K1nu21(0);
                            K1.umv(nu21, K1nu21);
                            FieldVector K3nu13(0);
                            K3.umv(nu13, K3nu13);
                            FieldVector K3nu23(0);
                            K3.umv(nu23, K3nu23);
                            Scalar g111 = lambda1 * (integrationOuterNormaln1 * K1nu11)/dF1;
                            Scalar g121 = lambda1 * (integrationOuterNormaln1 * K1nu21)/dF1;
                            Scalar g211 = lambda1 * (integrationOuterNormaln3 * K1nu11)/dF1;
                            Scalar g221 = lambda1 * (integrationOuterNormaln3 * K1nu21)/dF1;
                            Scalar g213 = alambda3 * (integrationOuterNormaln3 * K3nu13)/dF3;
                            Scalar g223 = alambda3 * (integrationOuterNormaln3 * K3nu23)/dF3;

                            // compute transmissibility matrix T = CA^{-1}B+F
                            FieldMatrix C(0), A(0), F(0), B(0);

                            // evaluate matrix C, F, A, B
                            C[0][0] = -g111;
                            C[0][1] = -g121;
                            C[1][1] = g223;

                            F[0][0] = g111 + g121;
                            F[1][1] = g213 - g223;

                            A[0][0] = g111;
                            A[0][1] = g121;
                            A[1][0] = g211;
                            A[1][1] = g223 + g221;

                            B[0][0] = g111 + g121;
                            B[1][0] = g211 + g221;
                            B[1][1] = g223 - g213;

                            // compute T
                            A.invert();
                            FieldMatrix CAinv(C.rightmultiply(A));
                            F += B.leftmultiply(CAinv);
                            FieldMatrix T(F);

                            // compute vector r
                            // evaluate r1, r2
                            FieldVector r1(0), r2(0);
                            r1[1] = -g213 * g2;
                            r2[0] = -J1 * face12vol/2.0;
                            r2[1] = g213 * g2;

                            // compute r = CA^{-1}r1
                            FieldVector r(0);
                            CAinv.umv(r2, r);
                            r += r1;

                            // assemble the global matrix this->A_ and right hand side this->f_
                            this->A_[globalIdx1][globalIdx1] += T[1][0];
                            this->A_[globalIdx1][globalIdx3] += T[1][1];
                            this->f_[globalIdx1] -= r[1];

                        }
                    }
                }
                // 'isIt' is on Dirichlet boundary

                else if (isItBcType.isDirichlet(pressEqIdx))
                {
                    // get Dirichlet boundary value
                    problem_.dirichlet(boundValues, *isIt);
                    Scalar g1 = boundValues[pressureIdx];

                    // compute total mobility for Dirichlet boundary 'isIt'
                    //determine lambda at the boundary -> if no saturation is known directly at the boundary use the cell saturation
                    Scalar alambda1 = 0;
                    if (isItBcType.isDirichlet(satEqIdx))
                    {
                        Scalar satBound = boundValues[saturationIdx];

                        //determine phase saturations from primary saturation variable
                        Scalar satW = 0;
                        switch (saturationType_)
                        {
                        case Sw:
                        {
                            satW = satBound;
                            break;
                        }
                        case Sn:
                        {
                            satW = 1 - satBound;
                            break;
                        }
                        }

                        Scalar lambdaWBound = 0;
                        Scalar lambdaNWBound = 0;

                        lambdaWBound = MaterialLaw::krw(
                                problem_.spatialParameters().materialLawParams(*eIt), satW)
                                / viscosity_[wPhaseIdx];
                        lambdaNWBound = MaterialLaw::krn(
                                problem_.spatialParameters().materialLawParams(*eIt), satW)
                                / viscosity_[nPhaseIdx];
                        alambda1 = lambdaWBound + lambdaNWBound;
                    }
                    else
                    {
                        alambda1 = lambda1;
                    }

                    // 'nextisIt' is on boundary
                    if (nextisIt->boundary())
                    {
                        // get boundary condition for boundary face center of 'nextisIt'
                        BoundaryTypes nextIsItBcType;
                        problem_.boundaryTypes(nextIsItBcType, *nextisIt);

                        // 'nextisIt': Dirichlet boundary
                        if (nextIsItBcType.isDirichlet(pressEqIdx))
                        {
                            // get Dirichlet boundary value of 'nextisIt'
                            problem_.dirichlet(boundValues, *nextisIt);
                            Scalar g3 = boundValues[pressureIdx];

                            // compute total mobility for Dirichlet boundary 'nextisIt'
                            //determine lambda at the boundary -> if no saturation is known directly at the boundary use the cell saturation
                            Scalar alambda1 = 0;
                            if (nextIsItBcType.isDirichlet(satEqIdx))
                            {
                                Scalar satBound = boundValues[saturationIdx];

                                //determine phase saturations from primary saturation variable
                                Scalar satW = 0;
                                switch (saturationType_)
                                {
                                case Sw:
                                {
                                    satW = satBound;
                                    break;
                                }
                                case Sn:
                                {
                                    satW = 1 - satBound;
                                    break;
                                }
                                }

                                Scalar lambdaWBound = 0;
                                Scalar lambdaNWBound = 0;

                                lambdaWBound = MaterialLaw::krw(
                                        problem_.spatialParameters().materialLawParams(*eIt), satW)
                                        / viscosity_[wPhaseIdx];
                                lambdaNWBound = MaterialLaw::krn(
                                        problem_.spatialParameters().materialLawParams(*eIt), satW)
                                        / viscosity_[nPhaseIdx];
                                alambda1 = lambdaWBound + lambdaNWBound;
                            }
                            else
                            {
                                alambda1 = lambda1;
                            }

                            // compute normal vectors nu11,nu21;
                            FieldVector nu11(0);
                            R.umv(globalPosFace13-globalPos1 ,nu11);

                            FieldVector nu21(0);
                            R.umv(globalPos1-globalPosFace12, nu21);

                            // compute dF1 i.e., the area of quadrilateral made by normal vectors 'nu'
                            FieldVector Rnu21(0);
                            R.umv(nu21, Rnu21);
                            Scalar dF1 = fabs(nu11 * Rnu21);

                            // compute components needed for flux calculation, denoted as 'g'
                            FieldVector K1nu11(0);
                            K1.umv(nu11, K1nu11);
                            FieldVector K1nu21(0);
                            K1.umv(nu21, K1nu21);

                            Scalar g111 = alambda1 * (integrationOuterNormaln1 * K1nu11)/dF1;
                            Scalar g121 = alambda1 * (integrationOuterNormaln1 * K1nu21)/dF1;
                            Scalar g211 = alambda1 * (integrationOuterNormaln3 * K1nu11)/dF1;
                            Scalar g221 = alambda1 * (integrationOuterNormaln3 * K1nu21)/dF1;

                            // evaluate T1, T3, r1, r3
                            Scalar T1 = g111 + g121;
                            Scalar T3 = g211 + g221;
                            Scalar r1 = g111 * g1 + g121 * g3;
                            Scalar r3 = g211 * g1 + g221 * g3;

                            // assemble matrix this->A_ and right hand side this->f_
                            this->A_[globalIdx1][globalIdx1] += T1 + T3;
                            this->f_[globalIdx1] += r1 + r3;
                        }
                        // 'nextisIt': Neumann boundary

                        else if (nextIsItBcType.isNeumann(pressEqIdx))
                        {
                            // get Neumann boundary value of 'nextisIt'
                            problem_.neumann(boundValues, *nextisIt);
                            Scalar J3 = (boundValues[wPhaseIdx]/density_[wPhaseIdx]+boundValues[nPhaseIdx]/density_[nPhaseIdx]);

                            // compute normal vectors nu11,nu21;
                            FieldVector nu11(0);
                            R.umv(globalPosFace13-globalPos1 ,nu11);

                            FieldVector nu21(0);
                            R.umv(globalPos1-globalPosFace12, nu21);

                            // compute dF1 i.e., the area of quadrilateral made by normal vectors 'nu'
                            FieldVector Rnu21(0);
                            R.umv(nu21, Rnu21);
                            Scalar dF1 = fabs(nu11 * Rnu21);

                            // compute components needed for flux calculation, denoted as 'g'
                            FieldVector K1nu11(0);
                            K1.umv(nu11, K1nu11);
                            FieldVector K1nu21(0);
                            K1.umv(nu21, K1nu21);
                            Scalar g111 = alambda1 * (integrationOuterNormaln1 * K1nu11)/dF1;
                            Scalar g121 = alambda1 * (integrationOuterNormaln1 * K1nu21)/dF1;
                            Scalar g211 = alambda1 * (integrationOuterNormaln3 * K1nu11)/dF1;
                            Scalar g221 = alambda1 * (integrationOuterNormaln3 * K1nu21)/dF1;

                            // evaluate T, r
                            Scalar T = g111 - g211 * g121/g221;
                            Scalar r = -T * g1 - g121 * (-J3) * nextisIt->geometry().volume()/ (2.0 * g221);

                            // assemble matrix this->A_ and right hand side this->f_
                            this->A_[globalIdx1][globalIdx1] += T;
                            this->f_[globalIdx1] -= r;
                        }
                    }
                    // 'nextisIt' is inside

                    else
                    {
                        // neighbor cell 3
                        // access neighbor cell 3
                        ElementPointer nextisItoutside = nextisIt->outside();
                        int globalIdx3 = problem_.variables().index(*nextisItoutside);

                        CellData& cellData3 = problem_.variables().cellData(globalIdx3);

                        // neighbor cell 3 geometry type
                        //Dune::GeometryType gt3 = nextisItoutside->geometry().type();

                        // get global coordinate of neighbor cell 3 center
                        GlobalPosition
                        globalPos3 = nextisItoutside->geometry().center();

                        // get absolute permeability of neighbor cell 3
                        FieldMatrix K3(problem_.spatialParameters().intrinsicPermeability(*nextisItoutside));

                        // get total mobility of neighbor cell 3
                        Scalar lambda3 = 0;
                        lambda3 = cellData3.mobility(wPhaseIdx) + cellData3.mobility(nPhaseIdx);

                        // get the information of the face 'isIt34' between cell3 and cell4 (locally numbered)
                        IntersectionIterator isIt34 = problem_.gridView().ibegin(*nextisItoutside);
                        IntersectionIterator innernextisItEnd = problem_.gridView().iend(*nextisItoutside);
                        for (IntersectionIterator innerisIt = problem_.gridView().ibegin(*nextisItoutside);
                                innerisIt != innernextisItEnd; ++innerisIt)
                        {
                            if (innerisIt->boundary())
                            {
                                for (int i=0; i<innerisIt->geometry().corners(); ++i)
                                {
                                    GlobalPosition innerisItcorner = innerisIt->geometry().corner(i);

                                    if (innerisItcorner == corner1234)
                                    {
                                        isIt34 = innerisIt;
                                        continue;
                                    }
                                }
                            }
                        }

                        // get geometry type of face 'isIt34'
                        //Dune::GeometryType gtf34 = isIt34->geometryInInside().type();

                        // center of face in global coordinates, i.e., the midpoint of edge 'isIt34'
                        GlobalPosition
                        globalPosFace34 = isIt34->geometry().center();

                        // get face volume
                        Scalar face34vol = isIt34->geometry().volume();

                        // get outer normal vector scaled with half volume of face 'isIt34'
                        Dune::FieldVector<Scalar,dimWorld> integrationOuterNormaln2
                        = isIt34->centerUnitOuterNormal();
                        integrationOuterNormaln2
                        *= face34vol/2.0;

                        // get boundary condition for boundary face (isIt34) center
                        BoundaryTypes isIt34BcType;
                        problem_.boundaryTypes(isIt34BcType, *isIt34);

                        // 'isIt34': Dirichlet boundary
                        if (isIt34BcType.isDirichlet(pressEqIdx))
                        {
                            // get Dirichlet boundary value of 'isIt34'
                            problem_.dirichlet(boundValues, *isIt34);
                            Scalar g2 = boundValues[pressureIdx];

                            // compute total mobility for Dirichlet boundary 'isIt34'
                            //determine lambda at the boundary -> if no saturation is known directly at the boundary use the cell saturation
                            Scalar alambda3 = 0;
                            if (isIt34BcType.isDirichlet(satEqIdx))
                            {
                                Scalar satBound = boundValues[saturationIdx];

                                //determine phase saturations from primary saturation variable
                                Scalar satW = 0;
                                switch (saturationType_)
                                {
                                case Sw:
                                {
                                    satW = satBound;
                                    break;
                                }
                                case Sn:
                                {
                                    satW = 1 - satBound;
                                    break;
                                }
                                }

                                Scalar lambdaWBound = 0;
                                Scalar lambdaNWBound = 0;

                                lambdaWBound = MaterialLaw::krw(
                                        problem_.spatialParameters().materialLawParams(*eIt), satW)
                                        / viscosity_[wPhaseIdx];
                                lambdaNWBound = MaterialLaw::krn(
                                        problem_.spatialParameters().materialLawParams(*eIt), satW)
                                        / viscosity_[nPhaseIdx];
                                alambda3 = lambdaWBound + lambdaNWBound;
                            }
                            else
                            {
                                alambda3 = lambda3;
                            }

                            // compute normal vectors nu11,nu21; nu13, nu23;
                            FieldVector nu11(0);
                            R.umv(globalPosFace13-globalPos1 ,nu11);

                            FieldVector nu21(0);
                            R.umv(globalPos1-globalPosFace12, nu21);

                            FieldVector nu13(0);
                            R.umv(globalPos3-globalPosFace13, nu13);

                            FieldVector nu23(0);
                            R.umv(globalPos3-globalPosFace34, nu23);

                            // compute dF1, dF3 i.e., the area of quadrilateral made by normal vectors 'nu'
                            FieldVector Rnu21(0);
                            R.umv(nu21, Rnu21);
                            Scalar dF1 = fabs(nu11 * Rnu21);

                            FieldVector Rnu23(0);
                            R.umv(nu23, Rnu23);
                            Scalar dF3 = fabs(nu13 * Rnu23);

                            // compute components needed for flux calculation, denoted as 'g'
                            FieldVector K1nu11(0);
                            K1.umv(nu11, K1nu11);
                            FieldVector K1nu21(0);
                            K1.umv(nu21, K1nu21);
                            FieldVector K3nu13(0);
                            K3.umv(nu13, K3nu13);
                            FieldVector K3nu23(0);
                            K3.umv(nu23, K3nu23);
                            Scalar g111 = alambda1 * (integrationOuterNormaln1 * K1nu11)/dF1;
                            Scalar g121 = alambda1 * (integrationOuterNormaln1 * K1nu21)/dF1;
                            Scalar g211 = alambda1 * (integrationOuterNormaln3 * K1nu11)/dF1;
                            Scalar g221 = alambda1 * (integrationOuterNormaln3 * K1nu21)/dF1;
                            Scalar g213 = alambda3 * (integrationOuterNormaln3 * K3nu13)/dF3;
                            Scalar g223 = alambda3 * (integrationOuterNormaln3 * K3nu23)/dF3;

                            // compute the matrix T & vector r
                            FieldMatrix T(0);
                            FieldVector r(0);

                            Scalar coe = g221 + g223;

                            // evaluate matrix T
                            T[0][0] = g111 + g121 * (g223 - g211)/coe;
                            T[0][1] = -g121 * (g223 - g213)/coe;
                            T[1][0] = g223 * (g211 + g221)/coe;
                            T[1][1] = -g221 * (g223 - g213)/coe;

                            // evaluate vector r
                            r[0] = -g111 * g1 + (g1 * g121 * g211 - g2 * g213 * g121)/coe;
                            r[1] = -(g1 * g211 * g223 + g2 * g221 * g213)/coe;
                            // assemble the global matrix this->A_ and right hand side this->f_
                            this->A_[globalIdx1][globalIdx1] += T[0][0] + T[1][0];
                            this->A_[globalIdx1][globalIdx3] += T[0][1] + T[1][1];
                            this->f_[globalIdx1] -= r[0] + r[1];

                        }
                        // 'isIt34': Neumann boundary

                        else if (isIt34BcType.isNeumann(pressEqIdx))
                        {
                            // get Neumann boundary value of 'isIt34'
                            problem_.neumann(boundValues, *isIt34);
                            Scalar J2 = (boundValues[wPhaseIdx]/density_[wPhaseIdx]+boundValues[nPhaseIdx]/density_[nPhaseIdx]);

                            // compute normal vectors nu11,nu21; nu13, nu23;
                            FieldVector nu11(0);
                            R.umv(globalPosFace13-globalPos1 ,nu11);

                            FieldVector nu21(0);
                            R.umv(globalPos1-globalPosFace12, nu21);

                            FieldVector nu13(0);
                            R.umv(globalPos3-globalPosFace13, nu13);

                            FieldVector nu23(0);
                            R.umv(globalPos3-globalPosFace34, nu23);

                            // compute dF1, dF3 i.e., the area of quadrilateral made by normal vectors 'nu'
                            FieldVector Rnu21(0);
                            R.umv(nu21, Rnu21);
                            Scalar dF1 = fabs(nu11 * Rnu21);

                            FieldVector Rnu23(0);
                            R.umv(nu23, Rnu23);
                            Scalar dF3 = fabs(nu13 * Rnu23);

                            // compute components needed for flux calculation, denoted as 'g'
                            FieldVector K1nu11(0);
                            K1.umv(nu11, K1nu11);
                            FieldVector K1nu21(0);
                            K1.umv(nu21, K1nu21);
                            FieldVector K3nu13(0);
                            K3.umv(nu13, K3nu13);
                            FieldVector K3nu23(0);
                            K3.umv(nu23, K3nu23);
                            Scalar g111 = alambda1 * (integrationOuterNormaln1 * K1nu11)/dF1;
                            Scalar g121 = alambda1 * (integrationOuterNormaln1 * K1nu21)/dF1;
                            Scalar g211 = alambda1 * (integrationOuterNormaln3 * K1nu11)/dF1;
                            Scalar g221 = alambda1 * (integrationOuterNormaln3 * K1nu21)/dF1;
                            Scalar g113 = lambda3 * (integrationOuterNormaln2 * K3nu13)/dF3;
                            Scalar g123 = lambda3 * (integrationOuterNormaln2 * K3nu23)/dF3;
                            Scalar g213 = lambda3 * (integrationOuterNormaln3 * K3nu13)/dF3;
                            Scalar g223 = lambda3 * (integrationOuterNormaln3 * K3nu23)/dF3;

                            // compute the matrix T & vector r in v = A^{-1}(Bu + r1) = Tu + r
                            FieldMatrix A(0), B(0);
                            FieldVector r1(0), r(0);

                            // evaluate matrix A, B
                            A[0][0] = g113;
                            A[0][1] = -g123;
                            A[1][0] = -g213;
                            A[1][1] = g221 + g223;

                            B[0][1] = g113 - g123;
                            B[1][0] = g211 + g221;
                            B[1][1] = g223 - g213;

                            // evaluate vector r1
                            r1[0] = -J2 * isIt34->geometry().volume()/2.0;
                            r1[1] = -g211 * g1;

                            // compute T and r
                            A.invert();
                            B.leftmultiply(A);
                            FieldMatrix T(B);
                            A.umv(r1, r);

                            // assemble the global matrix this->A_ and right hand side this->f_
                            this->A_[globalIdx1][globalIdx1] += (g111 + g121 - g121 * T[1][0]) + (g211 + g221 - g221 * T[1][0]);
                            this->A_[globalIdx1][globalIdx3] += -g121 * T[1][1] - g221 * T[1][1];
                            this->f_[globalIdx1] += (g111 + g211) * g1 + (g121 + g221) * r[1];

                        }
                    }
                }
            }

        } // end all intersections

    } // end grid traversal

//    // get the number of nonzero terms in the matrix
//    Scalar num_nonzero = 0;
//
//    // determine position of matrix entries
//    for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != problem_.gridView().template end<0>(); ++eIt)
//    {
//        // cell index
//        int globalIdxI = problem_.variables().index(*eIt);
//
//        if (this->A_[globalIdxI][globalIdxI] != 0)
//        ++num_nonzero;
//
//        // run through all intersections with neighbors
//        IntersectionIterator isItBegin = problem_.gridView().ibegin(*eIt);
//        IntersectionIterator isItEnd = problem_.gridView().iend(*eIt);
//        for (IntersectionIterator isIt = isItBegin; isIt!=isItEnd; ++isIt)
//        {
//            IntersectionIterator tempisIt = isIt;
//            IntersectionIterator tempisItBegin = isItBegin;
//
//            // 'nextisIt' iterates over next codimension 1 intersection neighboring with 'isIt'
//            // sequence of next is anticlockwise of 'isIt'
//            IntersectionIterator nextisIt = ++tempisIt;
//
//            // get 'nextisIt'
//            switch (GET_PROP_VALUE(TypeTag, GridImplementation))
//            {
//                // for SGrid
//                case GridTypeIndices::sGrid:
//                {
//                    if (nextisIt == isItEnd)
//                    {
//                        nextisIt = isItBegin;
//                    }
//                    else
//                    {
//                        nextisIt = ++tempisIt;
//
//                        if (nextisIt == isItEnd)
//                        {
//                            nextisIt = ++tempisItBegin;
//                        }
//                    }
//
//                    break;
//                }
//                // for YaspGrid
//                case GridTypeIndices::yaspGrid:
//                {
//                    if (nextisIt == isItEnd)
//                    {
//                        nextisIt = isItBegin;
//                    }
//                    else
//                    {
//                        nextisIt = ++tempisIt;
//
//                        if (nextisIt == isItEnd)
//                        {
//                            nextisIt = ++tempisItBegin;
//                        }
//                    }
//
//                    break;
//                }
//                // for UGGrid
//                case GridTypeIndices::ugGrid:
//                {
//                    if (nextisIt == isItEnd)
//                    nextisIt = isItBegin;
//
//                    break;
//                }
//                default:
//                {
//                    DUNE_THROW(Dune::NotImplemented, "GridType can not be used with MPFAO implementation!");
//                }
//            }
//
//            if (isIt->neighbor())
//            {
//                // access neighbor
//                ElementPointer outside = isIt->outside();
//                int globalIdxJ = problem_.variables().index(*outside);
//
//                if (this->A_[globalIdxI][globalIdxJ] != 0)
//                ++num_nonzero;
//            }
//
//            if (isIt->neighbor() && nextisIt->neighbor())
//            {
//                // access the common neighbor of isIt's and nextisIt's outside
//                ElementPointer outside = isIt->outside();
//                ElementPointer nextisItoutside = nextisIt->outside();
//
//                IntersectionIterator innerisItEnd = problem_.gridView().iend(*outside);
//                IntersectionIterator innernextisItEnd = problem_.gridView().iend(*nextisItoutside);
//
//                for (IntersectionIterator innerisIt = problem_.gridView().ibegin(*outside);
//                        innerisIt!=innerisItEnd; ++innerisIt )
//                for (IntersectionIterator innernextisIt = problem_.gridView().ibegin(*nextisItoutside);
//                        innernextisIt!=innernextisItEnd; ++innernextisIt)
//                {
//                    if (innerisIt->neighbor() && innernextisIt->neighbor())
//                    {
//                        ElementPointer innerisItoutside = innerisIt->outside();
//                        ElementPointer innernextisItoutside = innernextisIt->outside();
//
//                        if (innerisItoutside == innernextisItoutside && innerisItoutside != isIt->inside())
//                        {
//                            int globalIdxJ = problem_.variables().index(*innerisItoutside);
//
//                            if (this->A_[globalIdxI][globalIdxJ] != 0)
//                            ++num_nonzero;
//                        }
//                    }
//                }
//            }
//        } // end of 'for' IntersectionIterator
//    } // end of 'for' ElementIterator
//
//    std::cout << "number of nonzero terms in the MPFA O-matrix on level " << problem_.gridView().grid().maxLevel() <<" nnmat: " << num_nonzero << std::endl;

    return;
}

//constitutive functions are updated once if new saturations are calculated and stored in the variables object
template<class TypeTag>
void FVMPFAOPressure2P<TypeTag>::updateMaterialLaws()
{
    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = problem_.gridView().template end<0>();
    for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eItEnd; ++eIt)
    {
        int globalIdx = problem_.variables().index(*eIt);

        CellData& cellData = problem_.variables().cellData(globalIdx);

        //determine phase saturations from primary saturation variable

        Scalar satW = cellData.saturation(wPhaseIdx);
        Scalar satNW = cellData.saturation(nPhaseIdx);

        Scalar pc = MaterialLaw::pC(problem_.spatialParameters().materialLawParams(*eIt), satW);

            cellData.setSaturation(wPhaseIdx, satW);
            cellData.setSaturation(nPhaseIdx, satNW);
            cellData.setCapillaryPressure(pc);


        // initialize mobilities
        Scalar mobilityW = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(*eIt), satW) / viscosity_[wPhaseIdx];
        Scalar mobilityNW = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(*eIt), satW) / viscosity_[nPhaseIdx];

        // initialize mobilities
        cellData.setMobility(wPhaseIdx, mobilityW);
        cellData.setMobility(nPhaseIdx, mobilityNW);

        //initialize fractional flow functions
        cellData.setFracFlowFunc(wPhaseIdx, mobilityW / (mobilityW + mobilityNW));
        cellData.setFracFlowFunc(nPhaseIdx, mobilityNW / (mobilityW + mobilityNW));
    }
    return;
}


} // end of Dune namespace
#endif
