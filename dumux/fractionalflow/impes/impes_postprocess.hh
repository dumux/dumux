// $Id: impes.hh 972 2009-01-12 10:15:57Z lauser $

#ifndef DUNE_IMPESPOSTPROCESS_HH
#define DUNE_IMPESPOSTPROCESS_HH

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include "dumux/fractionalflow/impes/impes.hh"

/**
 * @file
 * @brief
 * @author Markus Wolff
 */

namespace Dune
{
/**
 * \ingroup fracflow
 * @brief IMplicit Pressure Explicit Saturation (IMPES) scheme for the solution of
 * coupled diffusion/transport problems
 */

template<class Grid, class Diffusion, class Transport, class VC> class IMPESPostProcess: public IMPES<
        Grid, Diffusion, Transport, VC>
{
    enum
    {
        dim = Grid::dimension, dimWorld = Grid::dimensionworld
    };

    typedef Dune::FractionalFlow<Grid, Diffusion, Transport, VC> FractionalFlow;
    typedef Dune::IMPES<Grid, Diffusion, Transport, VC> IMPES;

typedef    typename Grid::LevelGridView GridView;
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;
    typedef typename Grid::template Codim<1>::EntityPointer FacePointer;

    typedef typename Grid::template Codim<0>::HierarchicIterator HierarchicIterator;

    typedef typename FractionalFlow::Scalar Scalar;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;
    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;
    typedef Dune::BlockVector< FieldVector<FieldVector<Scalar, dim>, 2*dim> > VelType;

    Scalar calcQT()
    {
        const GridView& gridView = grid_.levelView(transLevel_);

        Scalar fluxCoarse = 0;

        ElementIterator eItCoarseEnd = gridView.template end<0>();
        for (ElementIterator eItCoarse = gridView.template begin<0>(); eItCoarse != eItCoarseEnd; ++eItCoarse)
        {
            int globalIdxICoarse = indexSetCoarse_.index(*eItCoarse);

            IntersectionIterator isItEndCoarse = gridView.template iend(*eItCoarse);
            for (IntersectionIterator isItCoarse = gridView.template ibegin(*eItCoarse); isItCoarse!=isItEndCoarse; ++isItCoarse)
            {
                // boundary face
                if (isItCoarse->boundary())
                {
                    int faceNumberCoarse = isItCoarse->numberInSelf();

                    Scalar faceAreaCoarse = (*isItCoarse).intersectionGlobal().volume();
                    // get normal vector scaled with volume

                    GeometryType faceGTCoarse = isItCoarse->intersectionSelfLocal().type();

                    const FieldVector<Scalar,dim-1>& faceLocalCoarse
                    = ReferenceElements<Scalar,dim-1>::general(faceGTCoarse).position(0,0);

                    const GlobalPosition& globalPosFaceCoarse= isItCoarse->intersectionGlobal().global(faceLocalCoarse); // globalPosFace coordinate of face center

                    const LocalPosition&
                    localPosFaceCoarse = ReferenceElements<Scalar,dim>::general(faceGTCoarse).position(faceNumberCoarse,1);

                    //get boundary condition for boundary face center
                    BoundaryConditions::Flags bctypePressCoarse = this->transProblem.bctypePress(globalPosFaceCoarse, *eItCoarse, localPosFaceCoarse);
                    BoundaryConditions::Flags bctypeSatCoarse = this->transProblem.bctypeSat(globalPosFaceCoarse, *eItCoarse, localPosFaceCoarse);

                    if (bctypePressCoarse == BoundaryConditions::dirichlet && bctypeSatCoarse == BoundaryConditions::neumann)
                    {
//                        std::cout<<"facearea = "<<faceAreaCoarse<<"velocity = "<<velocity_[globalIdxICoarse][faceNumberCoarse][0]<<std::endl;
                        fluxCoarse += velocity_[globalIdxICoarse][faceNumberCoarse][0]*faceAreaCoarse;
                    }
                }
            }// end intersection iteration
        } // end grid traversal

        return fluxCoarse;
    }

    Scalar calcQO()
    {
        const GridView& gridView = grid_.levelView(transLevel_);

        FieldVector<Scalar,dim> fluxCoarse(0);

        ElementIterator eItCoarseEnd = gridView.template end<0>();
        for (ElementIterator eItCoarse = gridView.template begin<0>(); eItCoarse != eItCoarseEnd; ++eItCoarse)
        {
            int globalIdxICoarse = indexSetCoarse_.index(*eItCoarse);

            if (transLevel_ == diffLevel_)
            {
                //get some cell properties
                GeometryType gt = eItCoarse->geometry().type();
                const LocalPosition&
                localPos = ReferenceElements<Scalar,dim>::general(gt).position(0,0);
                const GlobalPosition& globalPos = eItCoarse->geometry().global(localPos); //globalPosFace coordinates of cell center
                int globalIdxI = indexSetFine_.index(*eItCoarse); // index of fine-scale cell

                // run through all intersections with neighbors and boundary
                IntersectionIterator isItEnd = gridView.template iend(*eItCoarse);
                for (IntersectionIterator isIt = gridView.template ibegin(*eItCoarse); isIt!=isItEnd; ++isIt)
                {
                    // boundary face
                    if (isIt->boundary())
                    {
                        GeometryType faceGT = isIt->intersectionSelfLocal().type();

                        const FieldVector<Scalar,dim-1>& faceLocal
                        = ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

                        const GlobalPosition& globalPosFace= isIt->intersectionGlobal().global(faceLocal); // globalPosFace coordinate of face center

                        int faceNumberFine = isIt->numberInSelf();

                        // get pressure and permeability and total mobility in fine-scale element
                        Scalar pressI = pressure_[globalIdxI];

                        FieldMatrix permeabilityI = this->transProblem.soil().K(globalPos,*eItCoarse,localPos);

                        Scalar satI = saturation_[globalIdxICoarse];
                        Scalar lambdaI = this->transProblem.materialLaw().mobN((1-satI),globalPos,*eItCoarse,localPos);

                        Scalar faceAreaFine=isIt->intersectionGlobal().volume(); // volume of face

                        FieldVector<Scalar,dim> unitOuterNormal
                        = isIt->unitOuterNormal(faceLocal); // normal vector of unit length

                        const LocalPosition&
                        localPosFace = ReferenceElements<Scalar,dim>::general(faceGT).position(faceNumberFine,1);

                        FieldVector<Scalar,dim> velocityFine(0);

                        //get boundary condition for boundary face center
                        BoundaryConditions::Flags bctypePress = this->transProblem.bctypePress(globalPosFace, *eItCoarse, localPosFace);
                        BoundaryConditions::Flags bctypeSat = this->transProblem.bctypeSat(globalPosFace, *eItCoarse, localPosFace);

                        if (bctypePress == BoundaryConditions::dirichlet && bctypeSat == BoundaryConditions::neumann)
                        {
                            // distance vector between barycenters
                            FieldVector<Scalar,dimWorld> distVec = globalPos - globalPosFace;
                            Scalar dist = distVec.two_norm();

                            //normalise distVec for multiplication with the permeability
                            distVec /= dist;

                            // compute directed permeability vector permeabilityI.n
                            FieldVector<Scalar,dim> normalPermeabilityI(0);
                            permeabilityI.umv(distVec, normalPermeabilityI);

                            // compute averaged total mobility
                            Scalar meanLambda = lambdaI;

                            Scalar pressBound = this->transProblem.dirichletPress(globalPosFace, *eItCoarse, localPosFace);

                            velocityFine=normalPermeabilityI;
                            velocityFine *= (meanLambda * (pressBound-pressI) / dist);

                            fluxCoarse += (velocityFine*=faceAreaFine);
                        }
                    }
                }
            }
            else
            {
                HierarchicIterator eItEnd = eItCoarse-> hend(diffLevel_);
                for (HierarchicIterator eIt = eItCoarse->hbegin(diffLevel_); eIt != eItEnd; ++eIt)
                {
                    //only iterat through difflevel!!!
                    if (eIt->level() != diffLevel_) continue;

                    //get some cell properties
                    GeometryType gt = eIt->geometry().type();
                    const LocalPosition&
                    localPos = ReferenceElements<Scalar,dim>::general(gt).position(0,0);
                    const GlobalPosition& globalPos = eIt->geometry().global(localPos); //globalPosFace coordinates of cell center
                    int globalIdxI = indexSetFine_.index(*eIt); // index of fine-scale cell

                    // run through all intersections with neighbors and boundary
                    IntersectionIterator isItEnd = gridView.template iend(*eIt);
                    for (IntersectionIterator isIt = gridView.template ibegin(*eIt); isIt!=isItEnd; ++isIt)
                    {
                        // boundary face
                        if (isIt->boundary())
                        {
                            GeometryType faceGT = isIt->intersectionSelfLocal().type();

                            const FieldVector<Scalar,dim-1>& faceLocal
                            = ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

                            const GlobalPosition& globalPosFace= isIt->intersectionGlobal().global(faceLocal); // globalPosFace coordinate of face center

                            int faceNumberFine = isIt->numberInSelf();

                            // get pressure and permeability and total mobility in fine-scale element
                            Scalar pressI = pressure_[globalIdxI];

                            FieldMatrix permeabilityI = this->transProblem.soil().K(globalPos,*eIt,localPos);

                            Scalar satI = saturation_[globalIdxICoarse];
                            Scalar lambdaI = this->transProblem.materialLaw().mobN((1-satI),globalPos,*eIt,localPos);

                            Scalar faceAreaFine=isIt->intersectionGlobal().volume(); // volume of face

                            FieldVector<Scalar,dim> unitOuterNormal
                            = isIt->unitOuterNormal(faceLocal); // normal vector of unit length

                            const LocalPosition&
                            localPosFace = ReferenceElements<Scalar,dim>::general(faceGT).position(faceNumberFine,1);

                            FieldVector<Scalar,dim> velocityFine(0);

                            //get boundary condition for boundary face center
                            BoundaryConditions::Flags bctypePress = this->transProblem.bctypePress(globalPosFace, *eIt, localPosFace);
                            BoundaryConditions::Flags bctypeSat = this->transProblem.bctypeSat(globalPosFace, *eIt, localPosFace);

                            if (bctypePress == BoundaryConditions::dirichlet && bctypeSat == BoundaryConditions::neumann)
                            {
                                // distance vector between barycenters
                                FieldVector<Scalar,dimWorld> distVec = globalPos - globalPosFace;
                                Scalar dist = distVec.two_norm();

                                //normalise distVec for multiplication with the permeability
                                distVec /= dist;

                                // compute directed permeability vector permeabilityI.n
                                FieldVector<Scalar,dim> normalPermeabilityI(0);
                                permeabilityI.umv(distVec, normalPermeabilityI);

                                // compute averaged total mobility
                                Scalar meanLambda = lambdaI;

                                Scalar pressBound = this->transProblem.dirichletPress(globalPosFace, *eIt, localPosFace);

                                velocityFine=normalPermeabilityI;
                                velocityFine *= (meanLambda * (pressBound-pressI) / dist);

                                fluxCoarse += (velocityFine*=faceAreaFine);
                            }
                        }
                    }

                    // end intersection traversal
                }// end hierarchic iteration
            }

        } // end grid traversal
        return fluxCoarse[0];
    }

    void calcPoreVolume()
    {
        typedef typename Element::Geometry Geometry;

        const GridView& gridView = grid_.levelView(diffLevel_);

        totalPoreVolume_ = 0;

        ElementIterator eItEnd = gridView.template end<0>();
        ElementIterator eItBegin = gridView.template begin<0>();
        for (ElementIterator eIt = eItBegin; eIt != eItEnd; ++eIt)
        {
            const Geometry& geometry = eIt->geometry();

            GeometryType gt = geometry.type();

            // cell center in reference element
            const LocalPosition& localPos = ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

            const GlobalPosition& globalPos = geometry.global(localPos);

            // cell volume, assume linear map here
            Scalar volume = eIt->geometry().integrationElement(localPos)
            *Dune::ReferenceElements<Scalar,dim>::general(gt).volume();

            totalPoreVolume_ += (volume*this->transProblem.soil().porosity(globalPos,*eIt,localPos));
        }
        return;
    }

public:
    typedef typename IMPES::RepresentationType RepresentationType;

    void postProcessUpdate(Scalar t, Scalar dt)
    {
        k_++;
        Scalar qT = calcQT();

        for(int i=0;i<size_;i++)
        {
            pVI_[i] += (qT*dt)/totalPoreVolume_;
        }
        if (k_ % modulo_ == 0)
        {
            Scalar qO = calcQO();
            oilCut_=qO/qT;
//            std::cout<<"qO = "<<qO<<std::endl;
        }
//        std::cout<<"qT = "<<qT<<std::endl;

        return;
    }

    void vtkout(const char* name, int k) const
    {
        VTKWriter<GridView> vtkwriter(grid_.levelView(transLevel_));
        char fname[128];
        sprintf(fname, "%s-postprocess-%05d", name, k);
        vtkwriter.addCellData(pVI_, "PVI");
        vtkwriter.addCellData(oilCut_, "oil cut");
        vtkwriter.write(fname, VTKOptions::ascii);

        this->variables.vtkout(name, k);
        return;
    }

    //! Construct an IMPES object.
    IMPESPostProcess(Diffusion& diffusion, Transport& transport, int modulo = 1, int flag = 0, int nIt = 2,
            Scalar maxDef = 1e-5, Scalar om = 1) :
    IMPES(diffusion, transport, flag, nIt, maxDef, om), modulo_(modulo),
    diffLevel_(this->variables.diffLevel), transLevel_(this->variables.transLevel),size_(this->variables.transSize),
    velocity_(this->variables.velocity),pressure_(this->variables.pressure), saturation_(this->variables.saturation), grid_(this->variables.grid),
    pVI_(size_), oilCut_(size_),indexSetCoarse_(grid_.levelIndexSet(transLevel_)),indexSetFine_(grid_.levelIndexSet(diffLevel_)), k_(0)
    {
        calcPoreVolume();
        pVI_=0;
        oilCut_=1;
    }

private:
    int modulo_;
    int diffLevel_;
    int transLevel_;
    int size_;
    const VelType& velocity_;
    const RepresentationType& pressure_;
    const RepresentationType& saturation_;
    const Grid& grid_;
    RepresentationType pVI_;
    RepresentationType oilCut_;
    Scalar totalPoreVolume_;
    const IndexSet& indexSetCoarse_;
    const IndexSet& indexSetFine_;
    int k_;
};
}
#endif
