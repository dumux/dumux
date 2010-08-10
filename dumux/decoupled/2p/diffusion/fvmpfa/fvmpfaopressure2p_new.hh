// $Id: fvmpfaopressure2p_new.hh 3784 2010-06-24 13:43:57Z bernd $
/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff                                      *
 *   Copyright (C) 2007-2010 by Yufei Cao                                    *
 *   Institute of Applied Analysis and Numerical Simulation                  *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@mathematik.uni-stuttgart.de                   *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUMUX_FVMPFAOPRESSURE2P_NEW_HH
#define DUMUX_FVMPFAOPRESSURE2P_NEW_HH

// dune environent:
#include <dune/istl/bvector.hh>

// dumux environment
#include "dumux/common/exceptions.hh"
#include <dumux/decoupled/2p/2pproperties.hh>
#include <dumux/decoupled/2p/diffusion/fvmpfa/mpfaproperties.hh>
#include "interactionvolume.hh"

/**
 * @file
 * @brief  Base class for defining an instance of a numerical diffusion model
 * @brief  MPFA O-method
 * @brief  Remark1: only for 2-D quadrilateral grid.
 * @brief  Remark2: can use UGGrid or SGrid (YaspGrid); variable 'ch' is chosen to decide which grid will be used.
 * @brief  Remark3: without capillary pressure and gravity!
 * @author Yufei Cao
 */

namespace Dumux
{
//! \ingroup diffusion
//! Base class for defining an instance of a numerical diffusion model.
/*! An interface for defining a numerical diffusion model for the
 *  solution of equations of the form
 * \f$ - \text{div}\, (\lambda K \text{grad}\, p ) = 0, \f$,
 * \f$p = g\f$ on \f$\Gamma_1\f$, and
 * \f$-\lambda K \text{grad}\, p \cdot \mathbf{n} = J\f$
 * on \f$\Gamma_2\f$. Here,
 * \f$p\f$ denotes the pressure, \f$K\f$ the absolute permeability,
 * and \f$\lambda\f$ the total mobility, possibly depending on the
 * saturation.
 Template parameters are:

 - GridView a DUNE gridView type
 - Scalar type used for return values
 */
template<class TypeTag>
class FVMPFAOPressure2P
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP(TypeTag, PTAG(ReferenceElements)) ReferenceElements;
    typedef typename ReferenceElements::Container ReferenceElementContainer;
    typedef typename ReferenceElements::ContainerFaces ReferenceElementFaceContainer;
    typedef typename ReferenceElements::ReferenceElement ReferenceElement;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;
    typedef typename SpatialParameters::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridTypeIndices)) GridTypeIndices;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum
    {
        Sw = Indices::saturationW,
        Sn = Indices::saturationNW,
        vw = Indices::velocityW,
        vn = Indices::velocityNW,
        vt = Indices::velocityTotal
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
    };
    enum
    {
        globalCorner = 2, globalEdge = 3, neumannNeumann = 0, dirichletDirichlet = 1, dirichletNeumann = 2, neumannDirichlet = 3
    };

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> FieldMatrix;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureCoefficientMatrix)) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureRHSVector)) Vector;
    typedef Dune::FieldVector<Scalar, dim> FieldVector;

    typedef Dumux::InteractionVolume<TypeTag> InteractionVolume;

    typedef std::vector<InteractionVolume> GlobalInteractionVolumeVector;

    //initializes the matrix to store the system of equations
    void initializeMatrix();

    void storeInteractionVolumeInfo();

    //function which assembles the system of equations to be solved
    void assemble();

    //solves the system of equations to get the spatial distribution of the pressure
    void solve();

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
    void updateMaterialLaws(bool first);

    void initialize(bool solveTwice = true)
    {
        updateMaterialLaws(true);

        storeInteractionVolumeInfo();

        pressure();

        updateMaterialLaws(false);

        pressure();

        return;
    }

    // serialization methods
    template<class Restarter>
    void serialize(Restarter &res)
    {
        return;
    }

    template<class Restarter>
    void deserialize(Restarter &res)
    {
        return;
    }

    //! \brief Write data files
    /*  \param name file name */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        problem().variables().addOutputVtkFields(writer);
        return;
    }

    void pressure(bool solveTwice = true)
    {
        //        Dune::Timer timer;

        //        timer.reset();
        assemble();
        //        std::cout << "assembling MPFA O-matrix on level" << problem_.gridView().grid().maxLevel() << " took " << timer.elapsed() << " seconds" << std::endl;

        //        timer.reset();
        solve();
        //        std::cout << "solving MPFA O-matrix on level" << problem_.gridView().grid().maxLevel() << " took " << timer.elapsed() << " seconds" << std::endl;

        return;
    }

    FVMPFAOPressure2P(Problem& problem) :
        problem_(problem), A_(problem.variables().gridSize(), problem.variables().gridSize(), (4 * dim + (dim - 1))
                * problem.variables().gridSize(), Matrix::random), f_(problem.variables().gridSize()),
                 interactionVolumes_(problem_.gridView().size(dim))
    {
        initializeMatrix();
    }

private:
    Problem& problem_;
    Matrix A_;
    Vector f_;

protected:
    GlobalInteractionVolumeVector interactionVolumes_;

    static const int saturationType = GET_PROP_VALUE(TypeTag, PTAG(SaturationFormulation)); //!< gives kind of saturation used (\f$ 0 = S_w\f$, \f$ 1 = S_n\f$)
    static const int velocityType_ = GET_PROP_VALUE(TypeTag, PTAG(VelocityFormulation)); //!< gives kind of velocity used (\f$ 0 = v_w\f$, \f$ 1 = v_n\f$, \f$ 2 = v_t\f$)

};

template<class TypeTag>
void FVMPFAOPressure2P<TypeTag>::initializeMatrix()
{
    // determine matrix row sizes
    ElementIterator eItBegin = problem_.gridView().template begin<0> ();
    ElementIterator eItEnd = problem_.gridView().template end<0> ();
    for (ElementIterator eIt = eItBegin; eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = problem_.variables().index(*eIt);

        // initialize row size
        int rowSize = 1;

        // run through all intersections with neighbors
        IntersectionIterator isItBegin = problem_.gridView().ibegin(*eIt);
        IntersectionIterator isItEnd = problem_.gridView().iend(*eIt);
        for (IntersectionIterator isIt = isItBegin; isIt != isItEnd; ++isIt)
        {
            IntersectionIterator tempisIt = isIt;
            IntersectionIterator tempisItBegin = isItBegin;

            // 'nextIsIt' iterates over next codimension 1 intersection neighboring with 'isIt'
            IntersectionIterator nextIsIt = ++tempisIt;

            // get 'nextIsIt'
            switch (GET_PROP_VALUE(TypeTag, PTAG(GridImplementation)))
            {
            // for SGrid
            case GridTypeIndices::sGrid:
            {
                if (nextIsIt == isItEnd)
                    nextIsIt = isItBegin;
                else
                {
                    nextIsIt = ++tempisIt;

                    if (nextIsIt == isItEnd)
                    {
                        nextIsIt = ++tempisItBegin;
                    }
                }

                break;
            }
                // for YaspGrid
            case GridTypeIndices::yaspGrid:
            {
                if (nextIsIt == isItEnd)
                {
                    nextIsIt = isItBegin;
                }
                else
                {
                    nextIsIt = ++tempisIt;

                    if (nextIsIt == isItEnd)
                    {
                        nextIsIt = ++tempisItBegin;
                    }
                }

                break;
            }
                // for UGGrid
            case GridTypeIndices::ugGrid:
            {
                if (nextIsIt == isItEnd)
                    nextIsIt = isItBegin;

                break;
            }
            default:
            {
                DUNE_THROW(Dune::NotImplemented, "GridType can not be used with MPFAO implementation!");
            }
            }

            if (isIt->neighbor())
                rowSize++;

            if (isIt->neighbor() && nextIsIt->neighbor())
                rowSize++;
        } // end of 'for' IntersectionIterator

        // set number of indices in row globalIdxI to rowSize
        A_.setrowsize(globalIdxI, rowSize);

    } // end of 'for' ElementIterator

    // indicate that size of all rows is defined
    A_.endrowsizes();

    // determine position of matrix entries
    for (ElementIterator eIt = eItBegin; eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = problem_.variables().index(*eIt);

        // add diagonal index
        A_.addindex(globalIdxI, globalIdxI);

        // run through all intersections with neighbors
        IntersectionIterator isItBegin = problem_.gridView().ibegin(*eIt);
        IntersectionIterator isItEnd = problem_.gridView().iend(*eIt);
        for (IntersectionIterator isIt = isItBegin; isIt != isItEnd; ++isIt)
        {
            IntersectionIterator tempisIt = isIt;
            IntersectionIterator tempisItBegin = isItBegin;

            // 'nextIsIt' iterates over next codimension 1 intersection neighboring with 'isIt'
            // sequence of next is anticlockwise of 'isIt'
            IntersectionIterator nextIsIt = ++tempisIt;

            // get 'nextIsIt'
            switch (GET_PROP_VALUE(TypeTag, PTAG(GridImplementation)))
            {
            // for SGrid
            case GridTypeIndices::sGrid:
            {
                if (nextIsIt == isItEnd)
                {
                    nextIsIt = isItBegin;
                }
                else
                {
                    nextIsIt = ++tempisIt;

                    if (nextIsIt == isItEnd)
                    {
                        nextIsIt = ++tempisItBegin;
                    }
                }

                break;
            }
                // for YaspGrid
            case GridTypeIndices::yaspGrid:
            {
                if (nextIsIt == isItEnd)
                {
                    nextIsIt = isItBegin;
                }
                else
                {
                    nextIsIt = ++tempisIt;

                    if (nextIsIt == isItEnd)
                    {
                        nextIsIt = ++tempisItBegin;
                    }
                }

                break;
            }
                // for UGGrid
            case GridTypeIndices::ugGrid:
            {
                if (nextIsIt == isItEnd)
                    nextIsIt = isItBegin;

                break;
            }
            default:
            {
                DUNE_THROW(Dune::NotImplemented, "GridType can not be used with MPFAO implementation!");
            }
            }

            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer outside = isIt->outside();
                int globalIdxJ = problem_.variables().index(*outside);

                // add off diagonal index
                // add index (row,col) to the matrix
                A_.addindex(globalIdxI, globalIdxJ);
            }

            if (isIt->neighbor() && nextIsIt->neighbor())
            {
                // access the common neighbor of isIt's and nextIsIt's outside
                ElementPointer outside = isIt->outside();
                ElementPointer nextisItoutside = nextIsIt->outside();

                IntersectionIterator innerisItEnd = problem_.gridView().iend(*outside);
                IntersectionIterator innernextisItEnd = problem_.gridView().iend(*nextisItoutside);

                for (IntersectionIterator innerisIt = problem_.gridView().ibegin(*outside); innerisIt != innerisItEnd; ++innerisIt)
                    for (IntersectionIterator innernextisIt = problem_.gridView().ibegin(*nextisItoutside); innernextisIt
                            != innernextisItEnd; ++innernextisIt)
                    {
                        if (innerisIt->neighbor() && innernextisIt->neighbor())
                        {
                            ElementPointer innerisItoutside = innerisIt->outside();
                            ElementPointer innernextisItoutside = innernextisIt->outside();

                            if (innerisItoutside == innernextisItoutside && innerisItoutside != isIt->inside())
                            {
                                int globalIdxJ = problem_.variables().index(*innerisItoutside);

                                A_.addindex(globalIdxI, globalIdxJ);
                            }
                        }
                    }
            }
        } // end of 'for' IntersectionIterator
    } // end of 'for' ElementIterator

    // indicate that all indices are defined, check consistency
    A_.endindices();

    return;
}
//                 Indices used in a interaction volume of the MPFA-o method
//                 ___________________________________________________
//                 |                        |                        |
//                 | nuxy: cell geometry |       nxy: face normal |
//                 |       vectors (see MPFA) |
//                 |                        |                        |
//                 |            4-----------3-----------3 |
//                 |            | --> nu43 |   nu34 <--|            |
//                 |            | |nu41 1|--> n43 ||nu32 |
//                 |            | v ^     |0 ^   v|            |
//                 |____________4__0__|n14__|__n23_|_1__2____________|
//                 |            |    1 |     0 |            |
//                 |            | ^         |1 nu23 ^ |            |
//                 |            | |nu14 0|--> n12 | |            |
//                 |            | -->nu12 |   nu21<-- |            |
//                 |            1-----------1-----------2 |
//                 |    elementnumber |inter-                  |
//                 |                        |face-                   |
//                 |                        |number |
//                 |________________________|________________________|


// only for 2-D general quadrilateral
template<class TypeTag>
void FVMPFAOPressure2P<TypeTag>::storeInteractionVolumeInfo()
{
    // introduce matrix R for vector rotation and R is initialized as zero matrix
    FieldMatrix R(0);

    // evaluate matrix R
    if (dim == 2)
    {
        R[0][1] = 1;
        R[1][0] = -1;
    }

    // run through all elements
    ElementIterator eItEnd = problem_.gridView().template end<0> ();
    for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // get common geometry information for the following computation

        // get global coordinate of cell 1 center
        const GlobalPosition& globalPos1 = eIt->geometry().center();

        // get absolute permeability of neighbor cell 1
        FieldMatrix K1(problem_.spatialParameters().intrinsicPermeability(globalPos1, *eIt));

        IntersectionIterator isIt12Begin = problem_.gridView().ibegin(*eIt);
        IntersectionIterator isIt12End = problem_.gridView().iend(*eIt);
        for (IntersectionIterator isIt12 = isIt12Begin; isIt12 != isIt12End; ++isIt12)
        {
            // intersection iterator 'nextIsIt' is used to get geometry information
            IntersectionIterator tempIsIt = isIt12;
            IntersectionIterator tempIsItBegin = isIt12Begin;

            IntersectionIterator isIt14 = ++tempIsIt;

            //get isIt14
            switch (GET_PROP_VALUE(TypeTag, PTAG(GridImplementation)))
            {
            // for SGrid
            case GridTypeIndices::sGrid:
            {
                if (isIt14 == isIt12End)
                {
                    isIt14 = isIt12Begin;
                }
                else
                {
                    isIt14 = ++tempIsIt;

                    if (isIt14 == isIt12End)
                    {
                        isIt14 = ++tempIsItBegin;
                    }
                }

                break;
            }
                // for YaspGrid
            case GridTypeIndices::yaspGrid:
            {
                if (isIt14 == isIt12End)
                {
                    isIt14 = isIt12Begin;
                }
                else
                {
                    isIt14 = ++tempIsIt;

                    if (isIt14 == isIt12End)
                    {
                        isIt14 = ++tempIsItBegin;
                    }
                }

                break;
            }
                // for UGGrid
            case GridTypeIndices::ugGrid:
            {
                if (isIt14 == isIt12End)
                    isIt14 = isIt12Begin;

                break;
            }
            default:
            {
                DUNE_THROW(Dune::NotImplemented, "GridType can not be used with MPFAO implementation!");
            }
            }

            int indexInInside12 = isIt12->indexInInside();

            // get the intersection node /bar^{x_3} between 'isIt12' and 'isIt14', denoted as 'corner1234'
            // initialization of corner1234
            GlobalPosition corner1234(0);

            int globalVertIdx1234 = 0;

            // get the global coordinate and global vertex index of corner1234
            for (int i = 0; i < isIt12->geometry().corners(); ++i)
            {
                bool finished = false;

                const GlobalPosition& isIt12corner = isIt12->geometry().corner(i);

                for (int j = 0; j < isIt14->geometry().corners(); ++j)
                {
                    const GlobalPosition& isIt14corner = isIt14->geometry().corner(j);

                    if (isIt14corner == isIt12corner)
                    {
                        corner1234 = isIt12corner;

                        const ReferenceElement& referenceElement = ReferenceElementContainer::general(eIt->geometry().type());

                        int localVertIdx = referenceElement.subEntity(indexInInside12, dim - 1, i, dim);

                        //                        std::cout<<"localIdx = "<<localVertIdx<<"\n";
                        globalVertIdx1234 = problem_.variables().index(*((*eIt).template subEntity<dim> (localVertIdx)));

                        //                        std::cout<<"vertIdx = "<<globalVertIdx1234<<"\n";
                        //                        std::cout<<"position = "<<corner1234<<"\n";

                        finished = true;
                        break;
                    }
                }

                if (finished)
                {
                    break;
                }
            }

            if (interactionVolumes_[globalVertIdx1234].isStored())
            {
                continue;
            }
            else
            {
                interactionVolumes_[globalVertIdx1234].setStored();
                //                std::cout << "vertIdx = " << globalVertIdx1234 << "\n";
            }

            //store pointer 1
            ElementPointer ePtr(*eIt);
            interactionVolumes_[globalVertIdx1234].setSubVolumeElement(ePtr, 0);
            interactionVolumes_[globalVertIdx1234].setIndexOnElement(isIt12->indexInInside(), 0, 0);
            interactionVolumes_[globalVertIdx1234].setIndexOnElement(isIt14->indexInInside(), 0, 1);

            // center of face in global coordinates, i.e., the midpoint of edge 'isIt12'
            const GlobalPosition& globalPosFace12 = isIt12->geometry().center();

            // get face volume
            Scalar faceVol12 = isIt12->geometry().volume() / 2.0;

            // get outer normal vector scaled with half volume of face 'isIt12'
            Dune::FieldVector<Scalar, dimWorld> unitOuterNormal12 = isIt12->centerUnitOuterNormal();

            // center of face in global coordinates, i.e., the midpoint of edge 'isIt14'
            const GlobalPosition& globalPosFace41 = isIt14->geometry().center();

            // get face volume
            Scalar faceVol41 = isIt14->geometry().volume() / 2.0;

            // get outer normal vector scaled with half volume of face 'isIt14': for numbering of n see Aavatsmark, Eigestad
            Dune::FieldVector<Scalar, dimWorld> unitOuterNormal14 = isIt14->centerUnitOuterNormal();

            interactionVolumes_[globalVertIdx1234].setNormalTimesPerm(unitOuterNormal12, faceVol12, K1, 0, 0);
            interactionVolumes_[globalVertIdx1234].setNormalTimesPerm(unitOuterNormal14, faceVol41, K1, 0, 1);

            // compute normal vectors nu14,nu12
            FieldVector nu14(0);
            R.mv(globalPos1 - globalPosFace12, nu14);

            FieldVector nu12(0);
            R.mv(globalPosFace41 - globalPos1, nu12);

            interactionVolumes_[globalVertIdx1234].setNu(nu12, 0, 0);
            interactionVolumes_[globalVertIdx1234].setNu(nu14, 0, 1);

            // compute dF1, the area of quadrilateral made by normal vectors 'nu'
            FieldVector Rnu12(0);
            R.umv(nu12, Rnu12);
            interactionVolumes_[globalVertIdx1234].setDF(fabs(nu14 * Rnu12), 0);

            // handle interior face
            if (isIt12->neighbor())
            {
                // access neighbor cell 2 of 'isIt12'
                ElementPointer elementPointer2 = isIt12->outside();

                //store pointer 2
                interactionVolumes_[globalVertIdx1234].setSubVolumeElement(elementPointer2, 1);
                interactionVolumes_[globalVertIdx1234].setIndexOnElement(isIt12->indexInOutside(), 1, 1);

                // get global coordinate of neighbor cell 2 center
                const GlobalPosition& globalPos2 = elementPointer2->geometry().center();

                // get absolute permeability of neighbor cell 2
                FieldMatrix K2(problem_.spatialParameters().intrinsicPermeability(globalPos2, *elementPointer2));

                interactionVolumes_[globalVertIdx1234].setNormalTimesPerm(unitOuterNormal12, faceVol12, K2, 1, 1);

                // 'isIt14' is an interior face
                if (isIt14->neighbor())
                {
                    // neighbor cell 3
                    // access neighbor cell 3
                    ElementPointer elementPointer4 = isIt14->outside();

                    //store pointer 4
                    interactionVolumes_[globalVertIdx1234].setSubVolumeElement(elementPointer4, 3);
                    interactionVolumes_[globalVertIdx1234].setIndexOnElement(isIt14->indexInOutside(), 3, 0);

                    // get basic information of cell 1,2's neighbor cell 3,4

                    // get global coordinate of neighbor cell 4 center
                    const GlobalPosition& globalPos4 = elementPointer4->geometry().center();

                    // get absolute permeability of neighbor cell 2
                    FieldMatrix K4(problem_.spatialParameters().intrinsicPermeability(globalPos4, *elementPointer4));

                    interactionVolumes_[globalVertIdx1234].setNormalTimesPerm(unitOuterNormal14, faceVol41, K4, 3, 0);

                    // cell 3
                    GlobalPosition globalPos3(0);
                    int globalIdx3 = 0;

                    GlobalPosition globalPosFace23(0);
                    GlobalPosition globalPosFace34(0);

                    IntersectionIterator isIt2End = problem_.gridView().iend(*elementPointer2);
                    IntersectionIterator isIt4End = problem_.gridView().iend(*elementPointer4);
                    for (IntersectionIterator isIt2 = problem_.gridView().ibegin(*elementPointer2); isIt2 != isIt2End; ++isIt2)
                    {
                        bool finished = false;

                        for (IntersectionIterator isIt4 = problem_.gridView().ibegin(*elementPointer4); isIt4 != isIt4End; ++isIt4)
                        {
                            if (isIt2->neighbor() && isIt4->neighbor())
                            {
                                ElementPointer elementPointer32 = isIt2->outside();
                                ElementPointer elementPointer34 = isIt4->outside();

                                // find the common neighbor cell between cell 2 and cell 3, except cell 1
                                if (elementPointer32 == elementPointer34 && elementPointer32 != eIt)
                                {
                                    //store pointer 3
                                    interactionVolumes_[globalVertIdx1234].setSubVolumeElement(elementPointer32, 2);

                                    interactionVolumes_[globalVertIdx1234].setIndexOnElement(isIt2->indexInInside(), 1, 0);
                                    interactionVolumes_[globalVertIdx1234].setIndexOnElement(isIt2->indexInOutside(), 2, 1);
                                    interactionVolumes_[globalVertIdx1234].setIndexOnElement(isIt4->indexInInside(), 3, 1);
                                    interactionVolumes_[globalVertIdx1234].setIndexOnElement(isIt4->indexInOutside(), 2, 0);

                                    // access neighbor cell 4
                                    globalIdx3 = problem_.variables().index(*elementPointer32);

                                    // get global coordinate of neighbor cell 4 center
                                    globalPos3 = elementPointer32->geometry().center();

                                    globalPosFace23 = isIt2->geometry().center();
                                    globalPosFace34 = isIt4->geometry().center();

                                    Scalar faceVol23 = isIt2->geometry().volume() / 2.0;
                                    Scalar faceVol34 = isIt4->geometry().volume() / 2.0;

                                    // get outer normal vector scaled with half volume of face : for numbering of n see Aavatsmark, Eigestad
                                    FieldVector unitOuterNormal23 = isIt2->centerUnitOuterNormal();

                                    FieldVector unitOuterNormal43 = isIt4->centerUnitOuterNormal();

                                    // get absolute permeability of neighbor cell 2
                                    FieldMatrix K3(problem_.spatialParameters().intrinsicPermeability(globalPos3, *elementPointer32));

                                    interactionVolumes_[globalVertIdx1234].setNormalTimesPerm(unitOuterNormal23, faceVol23, K2, 1, 0);
                                    interactionVolumes_[globalVertIdx1234].setNormalTimesPerm(unitOuterNormal43, faceVol34, K3, 2, 0);
                                    interactionVolumes_[globalVertIdx1234].setNormalTimesPerm(unitOuterNormal23, faceVol23, K3, 2, 1);
                                    interactionVolumes_[globalVertIdx1234].setNormalTimesPerm(unitOuterNormal43, faceVol34, K4, 3, 1);

                                    // compute normal vectors nu23, nu21; nu32, nu34; nu41, nu43;
                                    FieldVector nu23(0);
                                    R.umv(globalPosFace12 - globalPos2, nu23);

                                    FieldVector nu21(0);
                                    R.umv(globalPosFace23 - globalPos2, nu21);

                                    FieldVector nu32(0);
                                    R.umv(globalPosFace34 - globalPos3, nu32);

                                    FieldVector nu34(0);
                                    R.umv(globalPos3 - globalPosFace23, nu34);

                                    FieldVector nu41(0);
                                    R.umv(globalPos4 - globalPosFace34, nu41);

                                    FieldVector nu43(0);
                                    R.umv(globalPos4 - globalPosFace41, nu43);

                                    interactionVolumes_[globalVertIdx1234].setNu(nu23, 1, 0);
                                    interactionVolumes_[globalVertIdx1234].setNu(nu21, 1, 1);
                                    interactionVolumes_[globalVertIdx1234].setNu(nu34, 2, 0);
                                    interactionVolumes_[globalVertIdx1234].setNu(nu32, 2, 1);
                                    interactionVolumes_[globalVertIdx1234].setNu(nu41, 3, 0);
                                    interactionVolumes_[globalVertIdx1234].setNu(nu43, 3, 1);

                                    // compute dF2, dF3, dF4 i.e., the area of quadrilateral made by normal vectors 'nu'
                                    FieldVector Rnu21(0);
                                    R.umv(nu21, Rnu21);
                                    interactionVolumes_[globalVertIdx1234].setDF(fabs(nu23 * Rnu21), 1);

                                    FieldVector Rnu34(0);
                                    R.umv(nu34, Rnu34);
                                    interactionVolumes_[globalVertIdx1234].setDF(fabs(nu32 * Rnu34), 2);

                                    FieldVector Rnu43(0);
                                    R.umv(nu43, Rnu43);
                                    interactionVolumes_[globalVertIdx1234].setDF(fabs(nu41 * Rnu43), 3);

                                    finished = true;

                                    break;
                                }
                            }
                        }
                        if (finished)
                        {
                            break;
                        }
                    }
                }

                // 'isIt14' is on the boundary
                else
                {
                    BoundaryConditions::Flags bcType14 = problem_.bctypePress(globalPosFace41, *isIt14);
                    interactionVolumes_[globalVertIdx1234].setBoundary(bcType14, 3);
                    if (bcType14 == BoundaryConditions::neumann)
                    {
                        std::vector<Scalar> fluxBC = problem_.neumannPress(globalPosFace41, *isIt14);
                        fluxBC[0] *= faceVol41;
                        fluxBC[1] *= faceVol41;
                        interactionVolumes_[globalVertIdx1234].setBoundaryCondition(fluxBC, 3);
                    }
                    else if (bcType14 == BoundaryConditions::dirichlet)
                    {
                        interactionVolumes_[globalVertIdx1234].setBoundaryCondition(problem_.dirichletPress(globalPosFace41, *isIt14), 3);
                    }

                    //check for dirichlet saturation conditions
                    interactionVolumes_[globalVertIdx1234].setSatBoundDirichlet(problem_.bctypeSat(globalPosFace41, *isIt14), 3);
                    if (interactionVolumes_[globalVertIdx1234].isDirichletSatBound(3))
                    {
                        interactionVolumes_[globalVertIdx1234].setDirichletSat(problem_.dirichletSat(globalPosFace41, *isIt14), 3);
                    }

                    // get common geometry information for the following computation
                    // get the information of the face 'isIt23' between cell2 and cell4 (locally numbered)

                    // center of face in global coordinates, i.e., the midpoint of edge 'isIt23'
                    GlobalPosition globalPosFace23(0);

                    // get face volume
                    Scalar faceVol23 = 0;

                    // get outer normal vector scaled with half volume of face 'isIt23'
                    Dune::FieldVector<Scalar, dimWorld> unitOuterNormal23(0);

                    bool finished = false;

                    IntersectionIterator isIt2End = problem_.gridView().iend(*elementPointer2);
                    for (IntersectionIterator isIt2 = problem_.gridView().ibegin(*elementPointer2); isIt2 != isIt2End; ++isIt2)
                    {
                        if (isIt2->boundary())
                        {
                            for (int i = 0; i < isIt2->geometry().corners(); ++i)
                            {
                                const GlobalPosition& corner2 = isIt2->geometry().corner(i);

                                if (corner2 == corner1234)
                                {
                                    interactionVolumes_[globalVertIdx1234].setIndexOnElement(isIt2->indexInInside(), 1, 0);

                                    globalPosFace23 = isIt2->geometry().center();

                                    faceVol23 = isIt2->geometry().volume() / 2.0;

                                    BoundaryConditions::Flags bcType23 = problem_.bctypePress(globalPosFace23, *isIt2);
                                    interactionVolumes_[globalVertIdx1234].setBoundary(bcType23, 1);
                                    if (bcType23 == BoundaryConditions::neumann)
                                    {
                                        std::vector<Scalar> fluxBC = problem_.neumannPress(globalPosFace23, *isIt2);
                                        fluxBC[0] *= faceVol23;
                                        fluxBC[1] *= faceVol23;
                                        interactionVolumes_[globalVertIdx1234].setBoundaryCondition(fluxBC, 1);
                                    }
                                    else if (bcType23 == BoundaryConditions::dirichlet)
                                    {
                                        interactionVolumes_[globalVertIdx1234].setBoundaryCondition(problem_.dirichletPress(
                                                globalPosFace23, *isIt2), 1);
                                    }
                                    interactionVolumes_[globalVertIdx1234].setBoundary(InteractionVolume::outside, 2);

                                    //check for dirichlet saturation conditions
                                    interactionVolumes_[globalVertIdx1234].setSatBoundDirichlet(
                                            problem_.bctypeSat(globalPosFace23, *isIt2), 1);
                                    if (interactionVolumes_[globalVertIdx1234].isDirichletSatBound(1))
                                    {
                                        interactionVolumes_[globalVertIdx1234].setDirichletSat(problem_.dirichletSat(globalPosFace23,
                                                *isIt2), 1);
                                    }

                                    unitOuterNormal23 = isIt2->centerUnitOuterNormal();

                                    interactionVolumes_[globalVertIdx1234].setNormalTimesPerm(unitOuterNormal23, faceVol23, K2, 1, 0);

                                    // compute normal vectors nu23, nu21;
                                    FieldVector nu23(0);
                                    R.umv(globalPosFace12 - globalPos2, nu23);

                                    FieldVector nu21(0);
                                    R.umv(globalPosFace23 - globalPos2, nu21);

                                    interactionVolumes_[globalVertIdx1234].setNu(nu23, 1, 0);
                                    interactionVolumes_[globalVertIdx1234].setNu(nu21, 1, 1);

                                    // compute dF2 i.e., the area of quadrilateral made by normal vectors 'nu'
                                    FieldVector Rnu21(0);
                                    R.umv(nu21, Rnu21);
                                    interactionVolumes_[globalVertIdx1234].setDF(fabs(nu23 * Rnu21), 1);

                                    finished = true;

                                    break;
                                }
                            }
                        }
                        if (finished)
                        {
                            break;
                        }
                    }
                    if (!finished)
                    {
                        DUNE_THROW(Dune::NotImplemented, "boundary shape not available as interaction volume shape");
                    }
                }
            }

            // handle boundary face 'isIt12'
            else
            {
                BoundaryConditions::Flags bcType12 = problem_.bctypePress(globalPosFace12, *isIt12);
                interactionVolumes_[globalVertIdx1234].setBoundary(bcType12, 0);
                if (bcType12 == BoundaryConditions::neumann)
                {
                    std::vector<Scalar> fluxBC = problem_.neumannPress(globalPosFace12, *isIt12);
                    fluxBC[0] *= faceVol12;
                    fluxBC[1] *= faceVol12;
                    interactionVolumes_[globalVertIdx1234].setBoundaryCondition(fluxBC, 0);
                }
                else if (bcType12 == BoundaryConditions::dirichlet)
                {
                    interactionVolumes_[globalVertIdx1234].setBoundaryCondition(problem_.dirichletPress(globalPosFace12, *isIt12), 0);
                }

                //check for dirichlet saturation conditions
                interactionVolumes_[globalVertIdx1234].setSatBoundDirichlet(problem_.bctypeSat(globalPosFace12, *isIt12), 0);
                if (interactionVolumes_[globalVertIdx1234].isDirichletSatBound(0))
                {
                    interactionVolumes_[globalVertIdx1234].setDirichletSat(problem_.dirichletSat(globalPosFace12, *isIt12), 0);
                }

                // 'isIt14' is on boundary
                if (isIt14->boundary())
                {
                    BoundaryConditions::Flags bcType41 = problem_.bctypePress(globalPosFace41, *isIt14);
                    interactionVolumes_[globalVertIdx1234].setBoundary(bcType41, 3);
                    if (bcType41 == BoundaryConditions::neumann)
                    {
                        std::vector<Scalar> fluxBC = problem_.neumannPress(globalPosFace41, *isIt14);
                        fluxBC[0] *= faceVol41;
                        fluxBC[1] *= faceVol41;
                        interactionVolumes_[globalVertIdx1234].setBoundaryCondition(fluxBC, 3);
                    }
                    else if (bcType41 == BoundaryConditions::dirichlet)
                    {
                        interactionVolumes_[globalVertIdx1234].setBoundaryCondition(problem_.dirichletPress(globalPosFace41, *isIt14), 3);
                    }
                    interactionVolumes_[globalVertIdx1234].setBoundary(InteractionVolume::outside, 1);
                    interactionVolumes_[globalVertIdx1234].setBoundary(InteractionVolume::outside, 2);

                    //check for dirichlet saturation conditions
                    interactionVolumes_[globalVertIdx1234].setSatBoundDirichlet(problem_.bctypeSat(globalPosFace41, *isIt14), 3);
                    if (interactionVolumes_[globalVertIdx1234].isDirichletSatBound(3))
                    {
                        interactionVolumes_[globalVertIdx1234].setDirichletSat(problem_.dirichletSat(globalPosFace41, *isIt14), 3);
                    }
                }

                // 'isIt14' is inside
                else
                {
                    // neighbor cell 3
                    // access neighbor cell 3
                    ElementPointer elementPointer4 = isIt14->outside();
                    interactionVolumes_[globalVertIdx1234].setIndexOnElement(isIt14->indexInOutside(), 3, 0);

                    //store pointer 4
                    interactionVolumes_[globalVertIdx1234].setSubVolumeElement(elementPointer4, 3);

                    // get global coordinate of neighbor cell 3 center
                    const GlobalPosition& globalPos4 = elementPointer4->geometry().center();

                    bool finished = false;

                    // get the information of the face 'isIt34' between cell3 and cell4 (locally numbered)
                    IntersectionIterator isIt4End = problem_.gridView().iend(*elementPointer4);
                    for (IntersectionIterator isIt4 = problem_.gridView().ibegin(*elementPointer4); isIt4 != isIt4End; ++isIt4)
                    {
                        if (isIt4->boundary())
                        {
                            for (int i = 0; i < isIt4->geometry().corners(); ++i)
                            {
                                const GlobalPosition& corner4 = isIt4->geometry().corner(i);

                                if (corner4 == corner1234)
                                {
                                    interactionVolumes_[globalVertIdx1234].setIndexOnElement(isIt4->indexInInside(), 3, 1);

                                    const GlobalPosition& globalPosFace34 = isIt4->geometry().center();

                                    Scalar faceVol34 = isIt4->geometry().volume() / 2.0;

                                    BoundaryConditions::Flags bcType34 = problem_.bctypePress(globalPosFace34, *isIt4);
                                    interactionVolumes_[globalVertIdx1234].setBoundary(bcType34, 2);
                                    if (bcType34 == BoundaryConditions::neumann)
                                    {
                                        std::vector<Scalar> fluxBC = problem_.neumannPress(globalPosFace34, *isIt4);
                                        fluxBC[0] *= faceVol34;
                                        fluxBC[1] *= faceVol34;
                                        interactionVolumes_[globalVertIdx1234].setBoundaryCondition(fluxBC, 2);
                                    }
                                    else if (bcType34 == BoundaryConditions::dirichlet)
                                    {
                                        interactionVolumes_[globalVertIdx1234].setBoundaryCondition(problem_.dirichletPress(
                                                globalPosFace34, *isIt4), 2);
                                    }
                                    interactionVolumes_[globalVertIdx1234].setBoundary(InteractionVolume::outside, 1);

                                    //check for dirichlet saturation conditions
                                    interactionVolumes_[globalVertIdx1234].setSatBoundDirichlet(
                                            problem_.bctypeSat(globalPosFace34, *isIt4), 2);
                                    if (interactionVolumes_[globalVertIdx1234].isDirichletSatBound(2))
                                    {
                                        interactionVolumes_[globalVertIdx1234].setDirichletSat(problem_.dirichletSat(globalPosFace34,
                                                *isIt4), 2);
                                    }

                                    FieldVector unitOuterNormal43 = isIt4->centerUnitOuterNormal();

                                    // get absolute permeability of neighbor cell 2
                                    FieldMatrix K4(problem_.spatialParameters().intrinsicPermeability(globalPos4, *elementPointer4));

                                    interactionVolumes_[globalVertIdx1234].setNormalTimesPerm(unitOuterNormal14, faceVol41, K4, 3, 0);
                                    interactionVolumes_[globalVertIdx1234].setNormalTimesPerm(unitOuterNormal43, faceVol34, K4, 3, 1);

                                    // compute normal vectors nu41, nu43;
                                    FieldVector nu41(0);
                                    R.umv(globalPos4 - globalPosFace34, nu41);

                                    FieldVector nu43(0);
                                    R.umv(globalPos4 - globalPosFace41, nu43);

                                    interactionVolumes_[globalVertIdx1234].setNu(nu41, 3, 0);
                                    interactionVolumes_[globalVertIdx1234].setNu(nu43, 3, 1);

                                    // compute dF1, dF3 i.e., the area of quadrilateral made by normal vectors 'nu'
                                    FieldVector Rnu43(0);
                                    R.umv(nu43, Rnu43);
                                    interactionVolumes_[globalVertIdx1234].setDF(fabs(nu41 * Rnu43), 3);

                                    finished = true;

                                    break;
                                }
                            }
                        }
                        if (finished)
                        {
                            break;
                        }
                    }
                    if (!finished)
                    {
                        DUNE_THROW(Dune::NotImplemented, "boundary shape not available as interaction volume shape");
                    }
                }
            }

        } // end all intersections
    } // end grid traversal

    return;
}

// only for 2-D general quadrilateral
template<class TypeTag>
void FVMPFAOPressure2P<TypeTag>::assemble()
{
    // initialization: set global matrix A_ to zero
    A_ = 0;
    f_ = 0;

    // run through all elements
    VertexIterator vItEnd = problem_.gridView().template end<dim> ();
    for (VertexIterator vIt = problem_.gridView().template begin<dim> (); vIt != vItEnd; ++vIt)
    {
        int globalVertIdx = problem_.variables().index(*vIt);

        std::vector<int> bcTypeFace(2 * dim);
        bcTypeFace[0] = interactionVolumes_[globalVertIdx].getBoundaryType(0);
        bcTypeFace[1] = interactionVolumes_[globalVertIdx].getBoundaryType(1);
        bcTypeFace[2] = interactionVolumes_[globalVertIdx].getBoundaryType(2);
        bcTypeFace[3] = interactionVolumes_[globalVertIdx].getBoundaryType(3);

        if (bcTypeFace[0] == InteractionVolume::inside && bcTypeFace[1] == InteractionVolume::inside && bcTypeFace[2]
                == InteractionVolume::inside && bcTypeFace[3] == InteractionVolume::inside)
        {

            ElementPointer& elementPointer1 = interactionVolumes_[globalVertIdx].getSubVolumeElement(0);
            ElementPointer& elementPointer2 = interactionVolumes_[globalVertIdx].getSubVolumeElement(1);
            ElementPointer& elementPointer3 = interactionVolumes_[globalVertIdx].getSubVolumeElement(2);
            ElementPointer& elementPointer4 = interactionVolumes_[globalVertIdx].getSubVolumeElement(3);

            // get global coordinate of cell centers
            const GlobalPosition& globalPos1 = elementPointer1->geometry().center();
            const GlobalPosition& globalPos2 = elementPointer2->geometry().center();
            const GlobalPosition& globalPos3 = elementPointer3->geometry().center();
            const GlobalPosition& globalPos4 = elementPointer4->geometry().center();

            // cell volumes
            Scalar volume1 = elementPointer1->geometry().volume();
            Scalar volume2 = elementPointer2->geometry().volume();
            Scalar volume3 = elementPointer3->geometry().volume();
            Scalar volume4 = elementPointer4->geometry().volume();

            // cell index
            int globalIdx1 = problem_.variables().index(*elementPointer1);
            int globalIdx2 = problem_.variables().index(*elementPointer2);
            int globalIdx3 = problem_.variables().index(*elementPointer3);
            int globalIdx4 = problem_.variables().index(*elementPointer4);

            // evaluate right hand side
            std::vector<Scalar> source(problem_.source(globalPos1, *elementPointer1));
            f_[globalIdx1] += volume1 / (4.0) * (source[wPhaseIdx] + source[nPhaseIdx]);
            source = problem_.source(globalPos2, *elementPointer2);
            f_[globalIdx2] += volume2 / (4.0) * (source[wPhaseIdx] + source[nPhaseIdx]);
            source = problem_.source(globalPos3, *elementPointer3);
            f_[globalIdx3] += volume3 / (4.0) * (source[wPhaseIdx] + source[nPhaseIdx]);
            source = problem_.source(globalPos4, *elementPointer4);
            f_[globalIdx4] += volume4 / (4.0) * (source[wPhaseIdx] + source[nPhaseIdx]);

            //compute total mobility of cell 1
            Scalar lambda1(problem_.variables().mobilityWetting(globalIdx1));
            lambda1 += problem_.variables().mobilityNonwetting(globalIdx1);

            Scalar lambda2(problem_.variables().mobilityWetting(globalIdx2));
            lambda2 += problem_.variables().mobilityNonwetting(globalIdx2);

            Scalar lambda3(problem_.variables().mobilityWetting(globalIdx3));
            lambda3 += problem_.variables().mobilityNonwetting(globalIdx3);

            Scalar lambda4(problem_.variables().mobilityWetting(globalIdx4));
            lambda4 += problem_.variables().mobilityNonwetting(globalIdx4);

            Scalar gn12nu14 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 0, 1);
            Scalar gn12nu12 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 0, 0);
            Scalar gn14nu14 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 1, 1);
            Scalar gn14nu12 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 1, 0);
            Scalar gn12nu23 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda2, 1, 1, 0);
            Scalar gn12nu21 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda2, 1, 1, 1);
            Scalar gn23nu23 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda2, 1, 0, 0);
            Scalar gn23nu21 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda2, 1, 0, 1);
            Scalar gn43nu32 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda3, 2, 0, 1);
            Scalar gn43nu34 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda3, 2, 0, 0);
            Scalar gn23nu32 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda3, 2, 1, 1);
            Scalar gn23nu34 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda3, 2, 1, 0);
            Scalar gn43nu41 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda4, 3, 1, 0);
            Scalar gn43nu43 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda4, 3, 1, 1);
            Scalar gn14nu41 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda4, 3, 0, 0);
            Scalar gn14nu43 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda4, 3, 0, 1);

            // compute transmissibility matrix T = CA^{-1}B+F
            Dune::FieldMatrix<Scalar, 2 * dim, 2 * dim> C(0), F(0), A(0), B(0);

            // evaluate matrix C, F, A, B
            C[0][0] = -gn12nu12;
            C[0][3] = -gn12nu14;
            C[1][0] = gn23nu21;
            C[1][1] = -gn23nu23;
            C[2][1] = gn43nu32;
            C[2][2] = gn43nu34;
            C[3][2] = -gn14nu43;
            C[3][3] = gn14nu41;

            F[0][0] = gn12nu12 + gn12nu14;
            F[1][1] = -gn23nu21 + gn23nu23;
            F[2][2] = -gn43nu34 - gn43nu32;
            F[3][3] = gn14nu43 - gn14nu41;

            A[0][0] = gn12nu12 + gn12nu21;
            A[0][1] = -gn12nu23;
            A[0][3] = gn12nu14;
            A[1][0] = -gn23nu21;
            A[1][1] = gn23nu23 + gn23nu32;
            A[1][2] = gn23nu34;
            A[2][1] = -gn43nu32;
            A[2][2] = -gn43nu34 - gn43nu43;
            A[2][3] = gn43nu41;
            A[3][0] = -gn14nu12;
            A[3][2] = gn14nu43;
            A[3][3] = -gn14nu41 - gn14nu14;

            B[0][0] = gn12nu12 + gn12nu14;
            B[0][1] = gn12nu21 - gn12nu23;
            B[1][1] = -gn23nu21 + gn23nu23;
            B[1][2] = gn23nu34 + gn23nu32;
            B[2][2] = -gn43nu34 - gn43nu32;
            B[2][3] = -gn43nu43 + gn43nu41;
            B[3][0] = -gn14nu12 - gn14nu14;
            B[3][3] = gn14nu43 - gn14nu41;

            // compute T
            A.invert();
            F += C.rightmultiply(B.leftmultiply(A));
            Dune::FieldMatrix<Scalar, 2 * dim, 2 * dim> T(F);

            // assemble the global matrix A_ and right hand side f
            A_[globalIdx1][globalIdx1] += T[0][0] + T[3][0];
            A_[globalIdx1][globalIdx2] += T[0][1] + T[3][1];
            A_[globalIdx1][globalIdx3] += T[0][2] + T[3][2];
            A_[globalIdx1][globalIdx4] += T[0][3] + T[3][3];

            A_[globalIdx2][globalIdx1] += -T[0][0] + T[1][0];
            A_[globalIdx2][globalIdx2] += -T[0][1] + T[1][1];
            A_[globalIdx2][globalIdx3] += -T[0][2] + T[1][2];
            A_[globalIdx2][globalIdx4] += -T[0][3] + T[1][3];

            A_[globalIdx3][globalIdx1] -= T[1][0] + T[2][0];
            A_[globalIdx3][globalIdx2] -= T[1][1] + T[2][1];
            A_[globalIdx3][globalIdx3] -= T[1][2] + T[2][2];
            A_[globalIdx3][globalIdx4] -= T[1][3] + T[2][3];

            A_[globalIdx4][globalIdx1] += T[2][0] - T[3][0];
            A_[globalIdx4][globalIdx2] += T[2][1] - T[3][1];
            A_[globalIdx4][globalIdx3] += T[2][2] - T[3][2];
            A_[globalIdx4][globalIdx4] += T[2][3] - T[3][3];
        }

        // at least one face on boundary!
        else
        {
            std::vector<int> interactionVolFaces(0);
            for (int faceIdx = 0; faceIdx < 2 * dim; faceIdx++)
            {
                if (bcTypeFace[faceIdx] != InteractionVolume::outside)
                {
                    interactionVolFaces.push_back(faceIdx);
                }
            }

            int numInteractionVolFaces = interactionVolFaces.size();

            switch (numInteractionVolFaces)
            {
            case globalCorner:
            {
                ElementPointer& elementPointer1 = interactionVolumes_[globalVertIdx].getSubVolumeElement(0);

                // get global coordinate of cell centers
                const GlobalPosition& globalPos1 = elementPointer1->geometry().center();

                // cell volumes
                Scalar volume1 = elementPointer1->geometry().volume();

                // cell index
                int globalIdx1 = problem_.variables().index(*elementPointer1);

                // evaluate right hand side
                std::vector<Scalar> source(problem_.source(globalPos1, *elementPointer1));
                f_[globalIdx1] += volume1 / (4.0) * (source[wPhaseIdx] + source[nPhaseIdx]);

                //get the densities
                Scalar densityW = this->problem().variables().densityWetting(globalIdx1);
                Scalar densityNW = this->problem().variables().densityNonwetting(globalIdx1);

                //compute total mobility of cell 1
                Scalar lambda1(problem_.variables().mobilityWetting(globalIdx1));
                lambda1 += problem_.variables().mobilityNonwetting(globalIdx1);

                Scalar gn12nu14 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 0, 1);
                Scalar gn12nu12 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 0, 0);
                Scalar gn14nu14 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 1, 1);
                Scalar gn14nu12 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 1, 0);

                //neumann - neumann
                if (bcTypeFace[interactionVolFaces[0]] == BoundaryConditions::neumann && bcTypeFace[interactionVolFaces[1]]
                        == BoundaryConditions::neumann)
                {
                    // compute transmissibility matrix T = CA^{-1}B+F
                    Dune::FieldMatrix<Scalar, dim, dim> C(0), F(0), A(0), B(0);
                    Dune::FieldVector<Scalar, dim> fB(0);

                    // evaluate matrix C, F, A, B
                    C[0][0] = -gn12nu12;
                    C[0][1] = -gn12nu14;

                    F[0][0] = gn12nu12 + gn12nu14;

                    A[0][0] = gn12nu12;
                    A[0][1] = gn12nu14;
                    A[1][0] = -gn14nu12;
                    A[1][1] = -gn14nu14;

                    B[0][0] = gn12nu12 + gn12nu14;
                    B[1][0] = -gn14nu12 - gn14nu14;

                    // get neumann boundary value
                    std::vector<Scalar> J(interactionVolumes_[globalVertIdx].template getBoundaryCondition<BoundaryConditions::neumann> (
                            interactionVolFaces[0]));
                    Scalar J1 = (J[wPhaseIdx] / densityW + J[nPhaseIdx] / densityNW);

                    J.swap(interactionVolumes_[globalVertIdx].template getBoundaryCondition<BoundaryConditions::neumann> (
                            interactionVolFaces[1]));
                    Scalar J4 = (J[wPhaseIdx] / densityW + J[nPhaseIdx] / densityNW);

                    fB[0] = J1;
                    fB[1] = -J4;

                    // compute T
                    A.invert();
                    C.rightmultiply(A);

                    F += B.leftmultiply(C);
                    Dune::FieldMatrix<Scalar, dim, dim> T(F);
                    Dune::FieldVector<Scalar, dim> r(0);
                    C.mv(fB, r);

                    // assemble the global matrix A_ and right hand side f
                    A_[globalIdx1][globalIdx1] += T[0][0] + T[1][0];

                    f_[globalIdx1] -= r[0] + r[1];
                }
                //dirichlet - dirichlet
                else if (bcTypeFace[interactionVolFaces[0]] == BoundaryConditions::dirichlet && bcTypeFace[interactionVolFaces[1]]
                        == BoundaryConditions::dirichlet)
                {
                    Scalar g1 = interactionVolumes_[globalVertIdx].template getBoundaryCondition<BoundaryConditions::dirichlet> (
                            interactionVolFaces[0]);
                    Scalar g4 = interactionVolumes_[globalVertIdx].template getBoundaryCondition<BoundaryConditions::dirichlet> (
                            interactionVolFaces[1]);

                    // assemble the global matrix A_ and right hand side f
                    A_[globalIdx1][globalIdx1] += gn12nu12 + gn12nu14 + gn14nu12 + gn14nu14;

                    f_[globalIdx1] -= gn12nu12 * g1 + gn12nu14 * g4 + gn14nu12 * g1 + gn14nu14 * g4;
                }
                //neumann - dirichlet
                else if (bcTypeFace[interactionVolFaces[0]] == BoundaryConditions::neumann && bcTypeFace[interactionVolFaces[1]]
                        == BoundaryConditions::dirichlet)
                {
                    // get neumann boundary value
                    std::vector<Scalar> J(interactionVolumes_[globalVertIdx].template getBoundaryCondition<BoundaryConditions::neumann> (
                            interactionVolFaces[0]));
                    Scalar J1 = (J[wPhaseIdx] / densityW + J[nPhaseIdx] / densityNW);
                    // get dirichlet boundary value
                    Scalar g4 = interactionVolumes_[globalVertIdx].template getBoundaryCondition<BoundaryConditions::dirichlet> (
                            interactionVolFaces[1]);

                    A_[globalIdx1][globalIdx1] += gn14nu14 - gn14nu12 * gn12nu14 / gn12nu12;

                    f_[globalIdx1] -= (gn14nu12 * gn12nu14 / gn12nu12 - gn14nu14) * g4 + gn14nu12 / gn12nu12 * J1;
                }
                //dirichlet - neumann
                else if (bcTypeFace[interactionVolFaces[0]] == BoundaryConditions::dirichlet && bcTypeFace[interactionVolFaces[1]]
                        == BoundaryConditions::neumann)
                {
                    // get dirichlet boundary value
                    Scalar g1 = interactionVolumes_[globalVertIdx].template getBoundaryCondition<BoundaryConditions::dirichlet> (
                            interactionVolFaces[0]);

                    // get neumann boundary value
                    std::vector<Scalar> J(interactionVolumes_[globalVertIdx].template getBoundaryCondition<BoundaryConditions::neumann> (
                            interactionVolFaces[1]));
                    Scalar J4 = (J[wPhaseIdx] / densityW + J[nPhaseIdx] / densityNW);

                    A_[globalIdx1][globalIdx1] += gn12nu12 - gn12nu14 * gn14nu12 / gn14nu14;

                    f_[globalIdx1] -= (gn12nu14 * gn14nu12 / gn14nu14 - gn12nu12) * g1 + gn12nu14 / gn14nu14 * J4;
                }
                else
                {
                    std::cout << interactionVolFaces[0] << ", " << interactionVolFaces[1] << ", " << interactionVolFaces[2] << ", "
                            << interactionVolFaces[3] << "\n";
                    DUNE_THROW(Dune::NotImplemented, "Boundary combination not supported in MPFA implementation");
                }

                break;
            }
            case globalEdge:
            {
                //numbering assures that either faces 2 and 4 or faces 1 and 3 are on the boundary
                //neumann - neumann

                if (bcTypeFace[interactionVolFaces[0]] == BoundaryConditions::neumann && bcTypeFace[interactionVolFaces[1]]
                        == BoundaryConditions::neumann)
                {
                    ElementPointer& elementPointer1 = interactionVolumes_[globalVertIdx].getSubVolumeElement(0);
                    ElementPointer& elementPointer4 = interactionVolumes_[globalVertIdx].getSubVolumeElement(3);

                    // get global coordinate of cell centers
                    const GlobalPosition& globalPos1 = elementPointer1->geometry().center();
                    const GlobalPosition& globalPos4 = elementPointer4->geometry().center();

                    // cell volumes
                    Scalar volume1 = elementPointer1->geometry().volume();
                    Scalar volume4 = elementPointer4->geometry().volume();

                    // cell index
                    int globalIdx1 = problem_.variables().index(*elementPointer1);
                    int globalIdx4 = problem_.variables().index(*elementPointer4);

                    // evaluate right hand side
                    std::vector<Scalar> source(problem_.source(globalPos1, *elementPointer1));
                    f_[globalIdx1] += volume1 / (4.0) * (source[wPhaseIdx] + source[nPhaseIdx]);
                    source = problem_.source(globalPos4, *elementPointer4);
                    f_[globalIdx4] += volume4 / (4.0) * (source[wPhaseIdx] + source[nPhaseIdx]);

                    //get the densities
                    Scalar densityW = this->problem().variables().densityWetting(globalIdx1);
                    Scalar densityNW = this->problem().variables().densityNonwetting(globalIdx1);

                    //compute total mobility of cell 1
                    Scalar lambda1(problem_.variables().mobilityWetting(globalIdx1));
                    lambda1 += problem_.variables().mobilityNonwetting(globalIdx1);

                    Scalar lambda4(problem_.variables().mobilityWetting(globalIdx4));
                    lambda4 += problem_.variables().mobilityNonwetting(globalIdx4);

                    Scalar gn12nu14 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 0, 1);
                    Scalar gn12nu12 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 0, 0);
                    Scalar gn14nu14 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 1, 1);
                    Scalar gn14nu12 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 1, 0);
                    Scalar gn43nu41 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda4, 3, 1, 0);
                    Scalar gn43nu43 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda4, 3, 1, 1);
                    Scalar gn14nu41 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda4, 3, 0, 0);
                    Scalar gn14nu43 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda4, 3, 0, 1);

                    // compute transmissibility matrix T = CA^{-1}B+F
                    Dune::FieldMatrix<Scalar, 2 * dim - 1, 2 * dim - 1> C(0), A(0);
                    Dune::FieldMatrix<Scalar, 2 * dim - 1, 2 * dim - 2> F(0), B(0);
                    Dune::FieldVector<Scalar, 2 * dim - 1> fB(0);

                    // evaluate matrix C, F, A, B
                    C[0][0] = -gn12nu12;
                    C[0][2] = -gn12nu14;
                    C[1][1] = -gn43nu43;
                    C[1][2] = gn43nu41;
                    C[2][1] = -gn14nu43;
                    C[2][2] = gn14nu41;

                    F[0][0] = gn12nu12 + gn12nu14;
                    F[1][1] = gn43nu43 - gn43nu41;
                    F[2][1] = gn14nu43 - gn14nu41;

                    A[0][0] = gn12nu12;
                    A[0][2] = gn12nu14;
                    A[1][1] = gn43nu43;
                    A[1][2] = -gn43nu41;
                    A[2][0] = -gn14nu12;
                    A[2][1] = gn14nu43;
                    A[2][2] = -gn14nu41 - gn14nu14;

                    B[0][0] = gn12nu12 + gn12nu14;
                    B[1][1] = gn43nu43 - gn43nu41;
                    B[2][0] = -gn14nu12 - gn14nu14;
                    B[2][1] = gn14nu43 - gn14nu41;

                    // get neumann boundary value
                    std::vector<Scalar> J(interactionVolumes_[globalVertIdx].template getBoundaryCondition<BoundaryConditions::neumann> (
                            interactionVolFaces[0]));

                    Scalar J1 = (J[wPhaseIdx] / densityW + J[nPhaseIdx] / densityNW);

                    J.swap(interactionVolumes_[globalVertIdx].template getBoundaryCondition<BoundaryConditions::neumann> (
                            interactionVolFaces[1]));
                    Scalar J4 = (J[wPhaseIdx] / densityW + J[nPhaseIdx] / densityNW);

                    fB[0] = -J1;
                    fB[1] = -J4;

                    Dune::FieldVector<Scalar, 2 * dim - 1> r(0);
                    // compute T
                    A.invert();
                    A.mv(fB, r);
                    C.mv(r, r);

                    Dune::FieldMatrix<Scalar, 2 * dim - 1, 2 * dim - 2> AinvB(A.rightmultiplyany(B));
                    Dune::FieldMatrix<Scalar, 2 * dim - 1, 2 * dim - 2> T(AinvB.leftmultiplyany(C));
                    T += F;

                    C.mv(fB, r);

                    // assemble the global matrix A_ and right hand side f
                    //flux 4
                    A_[globalIdx1][globalIdx1] += T[0][0] + T[2][0];
                    A_[globalIdx1][globalIdx4] += T[0][1] + T[2][1];

                    //fluxes 1 and 2
                    A_[globalIdx4][globalIdx1] += T[1][0] - T[2][0];
                    A_[globalIdx4][globalIdx4] += T[1][1] - T[2][1];

                    f_[globalIdx1] -= (r[0] + r[2]);
                    f_[globalIdx4] -= (r[1] - r[2]);

                }
                //neumann - neumann
                else if (bcTypeFace[interactionVolFaces[1]] == BoundaryConditions::neumann && bcTypeFace[interactionVolFaces[2]]
                        == BoundaryConditions::neumann)
                {
                    ElementPointer& elementPointer1 = interactionVolumes_[globalVertIdx].getSubVolumeElement(0);
                    ElementPointer& elementPointer2 = interactionVolumes_[globalVertIdx].getSubVolumeElement(1);

                    // get global coordinate of cell centers
                    const GlobalPosition& globalPos1 = elementPointer1->geometry().center();
                    const GlobalPosition& globalPos2 = elementPointer2->geometry().center();

                    // cell volumes
                    Scalar volume1 = elementPointer1->geometry().volume();
                    Scalar volume2 = elementPointer2->geometry().volume();

                    // cell index
                    int globalIdx1 = problem_.variables().index(*elementPointer1);
                    int globalIdx2 = problem_.variables().index(*elementPointer2);

                    // evaluate right hand side
                    std::vector<Scalar> source(problem_.source(globalPos1, *elementPointer1));
                    f_[globalIdx1] += volume1 / (4.0) * (source[wPhaseIdx] + source[nPhaseIdx]);
                    source = problem_.source(globalPos2, *elementPointer2);
                    f_[globalIdx2] += volume2 / (4.0) * (source[wPhaseIdx] + source[nPhaseIdx]);

                    //get the densities
                    Scalar densityW = this->problem().variables().densityWetting(globalIdx1);
                    Scalar densityNW = this->problem().variables().densityNonwetting(globalIdx1);

                    //compute total mobility of cell 1
                    Scalar lambda1(problem_.variables().mobilityWetting(globalIdx1));
                    lambda1 += problem_.variables().mobilityNonwetting(globalIdx1);

                    Scalar lambda2(problem_.variables().mobilityWetting(globalIdx2));
                    lambda2 += problem_.variables().mobilityNonwetting(globalIdx2);

                    Scalar gn12nu14 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 0, 1);
                    Scalar gn12nu12 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 0, 0);
                    Scalar gn14nu14 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 1, 1);
                    Scalar gn14nu12 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 1, 0);
                    Scalar gn12nu23 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda2, 1, 1, 0);
                    Scalar gn12nu21 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda2, 1, 1, 1);
                    Scalar gn23nu23 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda2, 1, 0, 0);
                    Scalar gn23nu21 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda2, 1, 0, 1);

                    // compute transmissibility matrix T = CA^{-1}B+F
                    Dune::FieldMatrix<Scalar, 2 * dim - 1, 2 * dim - 1> C(0), A(0);
                    Dune::FieldMatrix<Scalar, 2 * dim - 1, 2 * dim - 2> F(0), B(0);
                    Dune::FieldVector<Scalar, 2 * dim - 1> fB(0);

                    // evaluate matrix C, F, A, B
                    C[0][0] = -gn12nu12;
                    C[0][2] = -gn12nu14;
                    C[1][0] = gn23nu21;
                    C[1][1] = -gn23nu23;
                    C[2][0] = -gn14nu12;
                    C[2][2] = -gn14nu14;

                    F[0][0] = gn12nu12 + gn12nu14;
                    F[1][1] = -gn23nu21 + gn23nu23;
                    F[2][0] = gn14nu12 + gn14nu14;

                    A[0][0] = gn12nu12 + gn12nu21;
                    A[0][1] = -gn12nu23;
                    A[0][2] = gn12nu14;
                    A[1][0] = -gn23nu21;
                    A[1][1] = gn23nu23;
                    A[2][0] = gn14nu12;
                    A[2][2] = gn14nu14;

                    B[0][0] = gn12nu12 + gn12nu14;
                    B[0][1] = gn12nu21 - gn12nu23;
                    B[1][1] = -gn23nu21 + gn23nu23;
                    B[2][0] = gn14nu12 + gn14nu14;

                    // get neumann boundary value
                    std::vector<Scalar> J(interactionVolumes_[globalVertIdx].template getBoundaryCondition<BoundaryConditions::neumann> (
                            interactionVolFaces[1]));
                    Scalar J1 = (J[wPhaseIdx] / densityW + J[nPhaseIdx] / densityNW);

                    J.swap(interactionVolumes_[globalVertIdx].template getBoundaryCondition<BoundaryConditions::neumann> (
                            interactionVolFaces[2]));
                    Scalar J2 = (J[wPhaseIdx] / densityW + J[nPhaseIdx] / densityNW);

                    fB[1] = -J2;
                    fB[2] = -J1;

                    Dune::FieldVector<Scalar, 2 * dim - 1> r(0);
                    // compute T
                    A.invert();
                    A.mv(fB, r);
                    C.mv(r, r);

                    Dune::FieldMatrix<Scalar, 2 * dim - 1, 2 * dim - 2> AinvB(A.rightmultiplyany(B));
                    Dune::FieldMatrix<Scalar, 2 * dim - 1, 2 * dim - 2> T(AinvB.leftmultiplyany(C));

                    T += F;

                    C.mv(fB, r);

                    //flux 1
                    // assemble the global matrix A_ and right hand side f
                    A_[globalIdx1][globalIdx1] += T[0][0] + T[2][0];
                    A_[globalIdx1][globalIdx2] += T[0][1] + T[2][1];

                    //fluxes 1 and 2
                    A_[globalIdx2][globalIdx1] += -T[0][0] + T[1][0];
                    A_[globalIdx2][globalIdx2] += -T[0][1] + T[1][1];

                    f_[globalIdx1] -= (r[0] + r[2]);
                    f_[globalIdx2] -= (-r[0] + r[1]);

                }
                //dirichlet- dirichlet
                else if (bcTypeFace[interactionVolFaces[0]] == BoundaryConditions::dirichlet && bcTypeFace[interactionVolFaces[1]]
                        == BoundaryConditions::dirichlet)
                {
                    ElementPointer& elementPointer1 = interactionVolumes_[globalVertIdx].getSubVolumeElement(0);
                    ElementPointer& elementPointer4 = interactionVolumes_[globalVertIdx].getSubVolumeElement(3);

                    // get global coordinate of cell centers
                    const GlobalPosition& globalPos1 = elementPointer1->geometry().center();
                    const GlobalPosition& globalPos4 = elementPointer4->geometry().center();

                    // cell volumes
                    Scalar volume1 = elementPointer1->geometry().volume();
                    Scalar volume4 = elementPointer4->geometry().volume();

                    // cell index
                    int globalIdx1 = problem_.variables().index(*elementPointer1);
                    int globalIdx4 = problem_.variables().index(*elementPointer4);

                    // evaluate right hand side
                    std::vector<Scalar> source(problem_.source(globalPos1, *elementPointer1));
                    f_[globalIdx1] += volume1 / (4.0) * (source[wPhaseIdx] + source[nPhaseIdx]);
                    source = problem_.source(globalPos4, *elementPointer4);
                    f_[globalIdx4] += volume4 / (4.0) * (source[wPhaseIdx] + source[nPhaseIdx]);

                    //compute total mobility of cell 1
                    Scalar lambda1(problem_.variables().mobilityWetting(globalIdx1));
                    lambda1 += problem_.variables().mobilityNonwetting(globalIdx1);

                    Scalar lambda4(problem_.variables().mobilityWetting(globalIdx4));
                    lambda4 += problem_.variables().mobilityNonwetting(globalIdx4);

                    Scalar gn12nu14 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 0, 1);
                    Scalar gn12nu12 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 0, 0);
                    Scalar gn14nu14 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 1, 1);
                    Scalar gn14nu12 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 1, 0);
                    Scalar gn43nu41 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda4, 3, 1, 0);
                    Scalar gn43nu43 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda4, 3, 1, 1);
                    Scalar gn14nu41 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda4, 3, 0, 0);
                    Scalar gn14nu43 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda4, 3, 0, 1);

                    Scalar g1 = interactionVolumes_[globalVertIdx].template getBoundaryCondition<BoundaryConditions::dirichlet> (
                            interactionVolFaces[0]);
                    Scalar g4 = interactionVolumes_[globalVertIdx].template getBoundaryCondition<BoundaryConditions::dirichlet> (
                            interactionVolFaces[1]);

                    Dune::FieldMatrix<Scalar, 2 * dim - 1, 2 * dim - 2> T(0);

                    Scalar b = 1 / (gn14nu41 + gn14nu14);

                    T[0][0] = gn12nu12 + gn12nu14 - gn12nu14 * b * (gn14nu12 + gn14nu14);
                    T[0][1] = -gn12nu14 * b * (gn14nu41 - gn14nu43);
                    T[1][0] = gn43nu41 * b * (gn14nu12 + gn14nu14);
                    T[1][1] = gn43nu43 - gn43nu41 + gn43nu41 * b * (gn14nu41 - gn14nu43);
                    T[2][0] = gn14nu14 * b * (gn14nu12 + gn14nu14);
                    T[2][1] = gn14nu43 - gn14nu41 + gn14nu41 * b * (gn14nu41 - gn14nu43);

                    Dune::FieldVector<Scalar, 2 * dim - 1> r(0);

                    r[0] = -gn12nu12 * g1 - gn12nu14 * b * (gn14nu43 * g4 - gn14nu12 * g1);
                    r[1] = -gn43nu43 * g4 + gn43nu41 * b * (gn14nu43 * g4 - gn14nu12 * g1);
                    r[2] = -gn14nu43 * g4 + gn14nu41 * b * (gn14nu43 * g4 - gn14nu12 * g1);

                    //                    //fluxes 1 and 4
                    A_[globalIdx1][globalIdx1] += T[0][0] + T[2][0];
                    A_[globalIdx1][globalIdx4] += T[0][1] + T[2][1];

                    //fluxes 1 and 2
                    A_[globalIdx4][globalIdx1] += T[1][0] - T[2][0];
                    A_[globalIdx4][globalIdx4] += T[1][1] - T[2][1];

                    f_[globalIdx1] -= (r[0] + r[2]);
                    f_[globalIdx4] -= (r[1] - r[2]);

                }
                //dirichlet - dirichlet
                else if (bcTypeFace[interactionVolFaces[1]] == BoundaryConditions::dirichlet && bcTypeFace[interactionVolFaces[2]]
                        == BoundaryConditions::dirichlet)
                {
                    ElementPointer& elementPointer1 = interactionVolumes_[globalVertIdx].getSubVolumeElement(0);
                    ElementPointer& elementPointer2 = interactionVolumes_[globalVertIdx].getSubVolumeElement(1);

                    // get global coordinate of cell centers
                    const GlobalPosition& globalPos1 = elementPointer1->geometry().center();
                    const GlobalPosition& globalPos2 = elementPointer2->geometry().center();

                    // cell volumes
                    Scalar volume1 = elementPointer1->geometry().volume();
                    Scalar volume2 = elementPointer2->geometry().volume();

                    // cell index
                    int globalIdx1 = problem_.variables().index(*elementPointer1);
                    int globalIdx2 = problem_.variables().index(*elementPointer2);

                    // evaluate right hand side
                    std::vector<Scalar> source(problem_.source(globalPos1, *elementPointer1));
                    f_[globalIdx1] += volume1 / (4.0) * (source[wPhaseIdx] + source[nPhaseIdx]);
                    source = problem_.source(globalPos2, *elementPointer2);
                    f_[globalIdx2] += volume2 / (4.0) * (source[wPhaseIdx] + source[nPhaseIdx]);

                    //compute total mobility of cell 1
                    Scalar lambda1(problem_.variables().mobilityWetting(globalIdx1));
                    lambda1 += problem_.variables().mobilityNonwetting(globalIdx1);

                    Scalar lambda2(problem_.variables().mobilityWetting(globalIdx2));
                    lambda2 += problem_.variables().mobilityNonwetting(globalIdx2);

                    Scalar gn12nu14 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 0, 1);
                    Scalar gn12nu12 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 0, 0);
                    Scalar gn14nu14 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 1, 1);
                    Scalar gn14nu12 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 1, 0);
                    Scalar gn12nu23 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda2, 1, 1, 0);
                    Scalar gn12nu21 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda2, 1, 1, 1);
                    Scalar gn23nu23 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda2, 1, 0, 0);
                    Scalar gn23nu21 = interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda2, 1, 0, 1);

                    Scalar g1 = interactionVolumes_[globalVertIdx].template getBoundaryCondition<BoundaryConditions::dirichlet> (
                            interactionVolFaces[1]);
                    Scalar g2 = interactionVolumes_[globalVertIdx].template getBoundaryCondition<BoundaryConditions::dirichlet> (
                            interactionVolFaces[2]);

                    Dune::FieldMatrix<Scalar, 2 * dim - 1, 2 * dim - 2> T(0);

                    Scalar b = 1 / (gn12nu12 + gn12nu21);

                    T[0][0] = gn12nu12 + gn12nu14 - gn12nu12 * b * (gn12nu14 + gn12nu12);
                    T[0][1] = -gn12nu12 * b * (gn12nu21 - gn12nu23);
                    T[1][0] = gn23nu21 * b * (gn12nu14 + gn12nu12);
                    T[1][1] = -gn23nu21 + gn23nu23 + gn23nu21 * b * (gn12nu21 - gn12nu23);
                    T[2][0] = gn14nu12 + gn14nu14 - gn14nu12 * b * (gn12nu14 + gn12nu12);
                    T[2][1] = -gn14nu12 * b * (gn12nu21 - gn12nu23);

                    Dune::FieldVector<Scalar, 2 * dim - 1> r(0);

                    r[0] = -gn12nu14 * g1 - gn12nu12 * b * (gn12nu23 * g2 - gn12nu14 * g1);
                    r[1] = -gn23nu23 * g2 + gn23nu21 * b * (gn12nu23 * g2 - gn12nu14 * g1);
                    r[2] = -gn14nu14 * g1 - gn14nu12 * b * (gn12nu23 * g2 - gn12nu14 * g1);

                    //fluxes 1 and 4
                    A_[globalIdx1][globalIdx1] += T[0][0] + T[2][0];
                    A_[globalIdx1][globalIdx2] += T[0][1] + T[2][1];

                    //fluxes 1 and 2
                    A_[globalIdx2][globalIdx1] += -T[0][0] + T[1][0];
                    A_[globalIdx2][globalIdx2] += -T[0][1] + T[1][1];

                    f_[globalIdx1] -= (r[0] + r[2]);
                    f_[globalIdx2] -= (-r[0] + r[1]);
                }
                else
                {
                    DUNE_THROW(Dune::NotImplemented, "Boundary combination not supported in MPFA implementation");
                }

                break;
            }
            }

        } // end boundaries

    } // end vertex iterator

    // get the number of nonzero terms in the matrix
    //    Scalar num_nonzero = 0;
    //
    //    // determine position of matrix entries
    //    for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt != problem_.gridView().template end<0> (); ++eIt)
    //    {
    //        // cell index
    //        int globalIdxI = problem_.variables().index(*eIt);
    //
    //        if (A_[globalIdxI][globalIdxI] != 0)
    //            ++num_nonzero;
    //
    //        // run through all intersections with neighbors
    //        IntersectionIterator isItBegin = problem_.gridView().ibegin(*eIt);
    //        IntersectionIterator isItEnd = problem_.gridView().iend(*eIt);
    //        for (IntersectionIterator isIt = isItBegin; isIt != isItEnd; ++isIt)
    //        {
    //            IntersectionIterator tempIsIt = isIt;
    //            IntersectionIterator tempIsItBegin = isItBegin;
    //
    //            // 'nextIsIt' iterates over next codimension 1 intersection neighboring with 'isIt'
    //            // sequence of next is anticlockwise of 'isIt'
    //            IntersectionIterator nextIsIt = ++tempIsIt;
    //
    //            // get 'nextIsIt'
    //            switch (GET_PROP_VALUE(TypeTag, PTAG(GridImplementation)))
    //            {
    //            // for SGrid
    //            case GridTypeIndices::sGrid:
    //            {
    //                if (nextIsIt == isItEnd)
    //                {
    //                    nextIsIt = isItBegin;
    //                }
    //                else
    //                {
    //                    nextIsIt = ++tempIsIt;
    //
    //                    if (nextIsIt == isItEnd)
    //                    {
    //                        nextIsIt = ++tempIsItBegin;
    //                    }
    //                }
    //
    //                break;
    //            }
    //                // for YaspGrid
    //            case GridTypeIndices::yaspGrid:
    //            {
    //                if (nextIsIt == isItEnd)
    //                {
    //                    nextIsIt = isItBegin;
    //                }
    //                else
    //                {
    //                    nextIsIt = ++tempIsIt;
    //
    //                    if (nextIsIt == isItEnd)
    //                    {
    //                        nextIsIt = ++tempIsItBegin;
    //                    }
    //                }
    //
    //                break;
    //            }
    //                // for UGGrid
    //            case GridTypeIndices::ugGrid:
    //            {
    //                if (nextIsIt == isItEnd)
    //                    nextIsIt = isItBegin;
    //
    //                break;
    //            }
    //            default:
    //            {
    //                DUNE_THROW(Dune::NotImplemented, "GridType can not be used with MPFAO implementation!");
    //            }
    //            }
    //
    //            if (isIt->neighbor())
    //            {
    //                // access neighbor
    //                ElementPointer outside = isIt->outside();
    //                int globalIdxJ = problem_.variables().index(*outside);
    //
    //                if (A_[globalIdxI][globalIdxJ] != 0)
    //                    ++num_nonzero;
    //            }
    //
    //            if (isIt->neighbor() && nextIsIt->neighbor())
    //            {
    //                // access the common neighbor of isIt's and nextIsIt's outside
    //                ElementPointer outside = isIt->outside();
    //                ElementPointer nextisItoutside = nextIsIt->outside();
    //
    //                IntersectionIterator innerisItEnd = problem_.gridView().iend(*outside);
    //                IntersectionIterator innernextisItEnd = problem_.gridView().iend(*nextisItoutside);
    //
    //                for (IntersectionIterator innerisIt = problem_.gridView().ibegin(*outside); innerisIt != innerisItEnd; ++innerisIt)
    //                    for (IntersectionIterator innernextisIt = problem_.gridView().ibegin(*nextisItoutside); innernextisIt
    //                            != innernextisItEnd; ++innernextisIt)
    //                    {
    //                        if (innerisIt->neighbor() && innernextisIt->neighbor())
    //                        {
    //                            ElementPointer innerisItoutside = innerisIt->outside();
    //                            ElementPointer innernextisItoutside = innernextisIt->outside();
    //
    //                            if (innerisItoutside == innernextisItoutside && innerisItoutside != isIt->inside())
    //                            {
    //                                int globalIdxJ = problem_.variables().index(*innerisItoutside);
    //
    //                                if (A_[globalIdxI][globalIdxJ] != 0)
    //                                    ++num_nonzero;
    //                            }
    //                        }
    //                    }
    //            }
    //        } // end of 'for' IntersectionIterator
    //    } // end of 'for' ElementIterator
    //
    //    std::cout << "number of nonzero terms in the MPFA O-matrix on level " << problem_.gridView().grid().maxLevel() << " nnmat: "
    //            << num_nonzero << std::endl;

    return;
}

template<class TypeTag>
void FVMPFAOPressure2P<TypeTag>::solve()
{
    typedef typename GET_PROP(TypeTag, PTAG(SolverParameters)) SolverParameters;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressurePreconditioner)) Preconditioner;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureSolver)) Solver;

    Dune::MatrixAdapter<Matrix, Vector, Vector> op(A_); // make linear operator from A_
    Dune::InverseOperatorResult result;

    double reduction = SolverParameters::reductionSolver;
    int maxItSolver = SolverParameters::maxIterationNumberSolver;
    int iterPreconditioner = SolverParameters::iterationNumberPreconditioner;
    int verboseLevelSolver = SolverParameters::verboseLevelSolver;
    double relaxation = SolverParameters::relaxationPreconditioner;

    if (verboseLevelSolver)
    std::cout << "FVMPFAOPressure2P: solve for pressure" << std::endl;

    Preconditioner preconditioner(A_, iterPreconditioner, relaxation);
    Solver solver(op, preconditioner, reduction, maxItSolver, verboseLevelSolver);
    solver.apply(problem_.variables().pressure(), f_, result);

    //                        printmatrix(std::cout, A_, "global stiffness matrix", "row", 11, 3);
    //                        printvector(std::cout, f_, "right hand side", "row", 200, 1, 3);
    //        printvector(std::cout, (problem_.variables().pressure()), "pressure", "row", 200, 1, 3);

    return;
}

//constitutive functions are updated once if new saturations are calculated and stored in the variables object
template<class TypeTag>
void FVMPFAOPressure2P<TypeTag>::updateMaterialLaws(bool first = false)
{
    FluidState fluidState;

    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = problem_.gridView().template end<0> ();
    for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // get geometry type
        Dune::GeometryType gt = eIt->geometry().type();

        // get cell center in reference element
        const LocalPosition &localPos = ReferenceElementContainer::general(gt).position(0, 0);

        // get global coordinate of cell center
        const GlobalPosition& globalPos = eIt->geometry().global(localPos);

        int globalIdx = problem_.variables().index(*eIt);

        Scalar temperature = problem_.temperature(globalPos, *eIt);
        Scalar referencePressure = problem_.referencePressure(globalPos, *eIt);

        //determine phase saturations from primary saturation variable
        Scalar satW = 0;
        switch (saturationType)
        {
        case Sw:
        {
            satW = problem_.variables().saturation()[globalIdx];
            break;
        }
        case Sn:
        {
            satW = 1 - problem_.variables().saturation()[globalIdx];
            break;
        }
        }

        problem_.variables().capillaryPressure(globalIdx) = MaterialLaw::pC(
                problem_.spatialParameters().materialLawParams(globalPos, *eIt), satW);

        Scalar densityW = 0;
        Scalar densityNW = 0;
        Scalar viscosityW = 0;
        Scalar viscosityNW = 0;

        fluidState.update(satW, referencePressure, referencePressure, temperature);

        densityW = FluidSystem::phaseDensity(wPhaseIdx, temperature, referencePressure, fluidState);
        densityNW = FluidSystem::phaseDensity(nPhaseIdx, temperature, referencePressure, fluidState);

        viscosityW = FluidSystem::phaseViscosity(wPhaseIdx, temperature, referencePressure, fluidState);
        viscosityNW = FluidSystem::phaseViscosity(nPhaseIdx, temperature, referencePressure, fluidState);

        Scalar relPermW = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(globalPos, *eIt), satW);
        Scalar relPermNW = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(globalPos, *eIt), satW);

        Scalar mobilityW = relPermW / viscosityW;
        Scalar mobilityNW = relPermNW / viscosityNW;

        // initialize mobilities
        problem_.variables().mobilityWetting(globalIdx) = mobilityW;
        problem_.variables().mobilityNonwetting(globalIdx) = mobilityNW;

        // initialize densities
        problem_.variables().densityWetting(globalIdx) = densityW;
        problem_.variables().densityNonwetting(globalIdx) = densityNW;

        // initialize viscosities
        problem_.variables().viscosityWetting(globalIdx) = viscosityW;
        problem_.variables().viscosityNonwetting(globalIdx) = viscosityNW;

        problem_.spatialParameters().update(satW, *eIt);
    }
    return;
}

}
// end of Dune namespace
#endif
