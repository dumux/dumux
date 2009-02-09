// $Id$

#ifndef DUNE_RANDOMPERMEABILITY_HH
#define DUNE_RANDOMPERMEABILITY_HH

#include<iostream>
#include<vector>
#include<set>
#include<map>
#include<stdio.h>
#include<stdlib.h>

#include <dune/disc/functions/p0function.hh>

namespace Dune
{
/*! \brief providing the absolute permeability.
 *
 *  The class Permeability is derived from the template argument BV which usually
 *  respresents a block vector. The values for the permeability should be already set by
 *  the constructor.
 */
template<class Grid>
class RandomPermeability
{
    template<int dim>
    struct ElementLayout
    {
        bool contains(Dune::GeometryType gt)
        {
            return gt.dim() == dim;
        }
    };

    enum
    {
        dim = Grid::dimension, dimWorld = Grid::dimensionworld
    };
typedef    typename Grid::ctype Scalar;
    typedef typename Grid::LeafGridView GridView;
    typedef P0Function<GridView,Scalar,1> PermType;
    typedef BlockVector<FieldVector<Scalar,1> > RepresentationType;
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::IndexSet IndexSet;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<Grid,IndexSet,ElementLayout> ElementMapper;

    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*! \brief Constructor.
     *
     *  \param size number of degrees of freedom
     *  \param grid pointer to grid
     *  \param mapper pointer to mapper
     *
     *  The constructor already sets the entries to
     *  permeability values which are described by the functor \a permFunc,
     *  as for example given by PermeabilityBall or RandomPermeability. If the bool
     *  \a permFunc.random is set to true, \a permFunc is expected to set all entries in one call.
     *  Otherwise, a traversal over the cells is done, and \a permFunc should return the
     *  permeability at the cell center.
     */
    RandomPermeability(const Grid& grid, const char* name = "permeab.dat", const bool create = true)
    : grid_(grid), fileName_(name), createNew_(create), perm_(grid_.leafView()), permLoc_(0),
    elementMapper_(grid_, grid_.leafIndexSet())
    {
        typedef typename GridView::template Codim<0>::Iterator ElementIterator;

        const GridView& gridView(grid_.leafView());
        ElementIterator eItEnd = gridView.template end<0>();

        char* pwd(getenv("PWD"));
        char startCommand[220];
        strcpy(startCommand, "find ");
        strcat(startCommand, pwd);
        strcat(startCommand, "/");
        char systemCommand[220];

        char simSetDir[221];
        int foundSimSet = 0;
        int k = 0;
        while (!foundSimSet && k++ < 5)
        {
            strcpy(systemCommand, startCommand);
            strcat(systemCommand, " -type f -name simset -exec dirname {} \\; >simsetloc.txt");
            system(systemCommand);
            std::ifstream simSetLoc("simSetLoc.txt");
            simSetLoc.seekg (0, std::ios::end);
            int length = simSetLoc.tellg();
            simSetLoc.seekg (0, std::ios::beg);
            if (length> 0)
            {
                foundSimSet = 1;
                simSetLoc.getline(simSetDir, 220);
            }
            simSetLoc.close();
            strcat(startCommand, "../");
        }

        if (createNew_)
        {
            // SIMSET creates random permeabilities for given coordinates, so the coordinates of the center of gravity of each element
            // are written to a file 'SIMKOR'
            // open output stream for simset output file name
            char namefileName[100];
            strcpy(namefileName, simSetDir);
            strcat(namefileName, "/SIMNAM");
            std::ofstream namefile(namefileName);
            // Choose simset output filename
            namefile << fileName_ << std::endl;
            namefile.close();
            // open output stream for simset input file
            char outfileName[100];
            strcpy(outfileName, simSetDir);
            strcat(outfileName, "/SIMKOR");
            std::ofstream outfile(outfileName);
            for (ElementIterator eIt = gridView.template begin<0>(); eIt != eItEnd; ++eIt)
            {
                Dune::GeometryType gt = eIt->geometry().type();

                const LocalPosition&
                localPos = Dune::ReferenceElements<Scalar,dim>::general(gt).position(0,0);

                // get global coordinate of cell center
                const GlobalPosition& globalPos = eIt->geometry().global(localPos);

                outfile << globalPos[0] << "\t" << globalPos[1] << std::endl;
            }
            outfile.close();
            strcpy(systemCommand, "cd ");
            strcat(systemCommand, simSetDir);
            strcat(systemCommand, "; ./simset; cd $OLDPWD");
            system(systemCommand);
        }

        // open input stream for simset output file
        char concd[100];
        strcpy (concd, simSetDir);
        strcat(concd, "/");
        std::ifstream infile(strcat(concd, fileName_));
        std::cout << "Read permeability data from " << concd << std::endl;
        for (ElementIterator eIt = gridView.template begin<0>(); eIt != eItEnd; ++eIt)
        {
            int globalIdxI = elementMapper_.map(*eIt);
            Scalar dummy1, dummy2, permi;
            char zeile [221];
            infile.getline(zeile, 220);
            std::istringstream ist(zeile);
            ist >> dummy1 >> dummy2 >> permi;
            (*perm_)[globalIdxI] = pow(10.0, permi);
        }
        infile.close();
    }

    //! return const reference to permeability vector
    const RepresentationType& operator* () const
    {
        return (*perm_);
    }

    //! return reference to permeability vector
    RepresentationType& operator* ()
    {
        return (*perm_);
    }

    Dune::FieldMatrix<Scalar,dim,dim>& K (const Element& e)
    {
        int elemId = elementMapper_.map(e);
        Scalar permE = (*perm_)[elemId];

        for (int i = 0; i < dim; i++)
        permLoc_[i][i] = permE;

        return permLoc_;
    }

    void vtkout (const char* name, const Grid& grid_) const
    {
        Dune::VTKWriter<typename Grid::LeafGridView>
        vtkwriter(grid_.leafView());
        vtkwriter.addCellData(*perm_, "absolute permeability");
        int size = (*perm_).size();
        RepresentationType logPerm(size);
        for (int i = 0; i < size; i++)
        logPerm[i] = log10((*perm_)[i]);
        vtkwriter.addCellData(logPerm, "logarithm of permeability");
        vtkwriter.write(name, Dune::VTKOptions::ascii);
    }

private:
    const Grid& grid_;
    PermType perm_;
    Dune::FieldMatrix<Scalar,dim,dim> permLoc_;
    const bool createNew_;
    const char* fileName_;
    ElementMapper elementMapper_;
};

/*! \brief providing the absolute permeability for cells on given level and their children.
 *
 *  Unlike the class RandomPermeability
 *  which provides the permeability for the leaf grid_, in LevelRandomPermeability the
 *  permeability field is provided on a given grid_ level \f$ l \f$. The permeability
 *  of the level-\f$ l \f$ elements is also inherited to their children.
 */
template<class Grid>
class LevelRandomPermeability
{
    template<int dim>
    struct ElementLayout
    {
        bool contains (Dune::GeometryType gt)
        {
            return gt.dim() == dim;
        }
    };

    enum
    {   dim = Grid::dimension, dimWorld = Grid::dimensionworld};
    typedef typename Grid::LevelGridView GridView;
    typedef typename Grid::ctype Scalar;
    typedef LevelP0Function<Grid,Scalar,1> PermType;
    typedef BlockVector<FieldVector<Scalar,1> > RepresentationType;
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::Traits::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IndexSet IndexSet;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<Grid,IndexSet,ElementLayout> ElementMapper;

    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*! \brief Constructor.
     *
     *  \param g a grid_ objectof type Grid
     *  \param lev the level on which the permeability is to be provided
     *  \param name the name of the file in the simset-directory in which the permeabilities are to be stored.
     *  \param create set true if new field shall be created, set false if permeabilities shall be read from specified file.
     *
     *  The constructor already sets the entries to
     *  permeability values which are described by the functor \a permFunc,
     *  as for example given by PermeabilityBall or RandomPermeability. If the bool
     *  \a permFunc.random is set to true, \a permFunc is expected to set all entries in one call.
     *  Otherwise, a traversal over the cells is done, and \a permFunc should return the
     *  permeability at the cell center.
     */
    LevelRandomPermeability(const Grid& grid, const int lev, const char* name = "permeab.dat", const bool create = true)
    : grid_(grid), level_(lev), fileName_(name), createNew_(create), perm_(grid_,lev), permLoc_(0),
    elementMapper_(grid_, grid_.levelIndexSet(level_))
    {
        if (level_> grid_.maxLevel() ) DUNE_THROW(Dune::Exception,"Level specified for permeability data is higher than maximum grid_ level!");
        typedef typename GridView::template Codim<0>::Iterator ElementIterator;

        const GridView& gridView(grid_.levelView(level()));
        ElementIterator eItEnd = gridView.template end<0>();

        char* pwd(getenv("PWD"));
        char startCommand[220];
        strcpy(startCommand, "find ");
        strcat(startCommand, pwd);
        strcat(startCommand, "/");
        char systemCommand[220];

        char simSetDir[221];
        int foundSimSet = 0;
        int k = 0;
        while (!foundSimSet && k++ < 5)
        {
            strcpy(systemCommand, startCommand);
            strcat(systemCommand, " -type f -name simset -exec dirname {} \\; >simSetLoc.txt");
            system(systemCommand);
            std::ifstream simSetLoc("simSetLoc.txt");
            simSetLoc.seekg (0, std::ios::end);
            int length = simSetLoc.tellg();
            simSetLoc.seekg (0, std::ios::beg);
            if (length> 0)
            {
                foundSimSet = 1;
                simSetLoc.getline(simSetDir, 220);
            }
            simSetLoc.close();
            strcat(startCommand, "../");
        }

        if (createNew_)
        {
            // SIMSET creates random permeabilities for given coordinates, so the coordinates of the center of gravity of each element
            // are written to a file 'SIMKOR'
            // open output stream for simset output file name
            char namefileName[100];
            strcpy(namefileName, simSetDir);
            strcat(namefileName, "/SIMNAM");
            std::ofstream namefile(namefileName);
            // Choose simset output filename
            namefile << fileName_ << std::endl;
            namefile.close();
            // open output stream for simset input file
            char outfileName[100];
            strcpy(outfileName, simSetDir);
            strcat(outfileName, "/SIMKOR");
            std::ofstream outfile(outfileName);
            for (ElementIterator eIt = gridView.template begin<0>(); eIt != eItEnd; ++eIt)
            {
                Dune::GeometryType gt = eIt->geometry().type();

                const LocalPosition&
                localPos = Dune::ReferenceElements<Scalar,dim>::general(gt).position(0,0);

                // get global coordinate of cell center
                const GlobalPosition& globalPos = eIt->geometry().global(localPos);

                outfile << globalPos[0] << "\t" << globalPos[1] << std::endl;
            }
            outfile.close();
            strcpy(systemCommand, "cd ");
            strcat(systemCommand, simSetDir);
            strcat(systemCommand, "; ./simset; cd $OLDPWD");
            system(systemCommand);
        }

        // open input stream for simset output file
        char concd[100];
        strcpy (concd, simSetDir);
        strcat(concd, "/");
        std::ifstream infile(strcat(concd, fileName_));
        std::cout << "Read permeability data from " << concd << std::endl;
        for (ElementIterator eIt = gridView.template begin<0>(); eIt != eItEnd; ++eIt)
        {
            int globalIdxI = elementMapper_.map(*eIt);
            Scalar dummy1, dummy2, permi;
            char zeile [221];
            infile.getline(zeile, 220);
            std::istringstream ist(zeile);
            ist >> dummy1 >> dummy2 >> permi;
            (*perm_)[globalIdxI] = pow(10.0, permi);
        }
        infile.close();
    }

    //! return const reference to permeability vector
    const RepresentationType& operator* () const
    {
        return (*perm_);
    }

    //! return reference to permeability vector
    RepresentationType& operator* ()
    {
        return (*perm_);
    }

    //! \brief return reference to permeability tensor of specified cell.
    /** \param e cell of level\f$ l \f$ or higher
     *
     */
    Dune::FieldMatrix<Scalar,dim,dim>& K (const Element& e)
    {
        int le = e.level();
        int elemId;
        if (le < level_) DUNE_THROW(Dune::Exception, "Level of element lower than level of permeability discretisation, permeability not uniquely defined");
        else if (le> level_)
        {
            ElementPointer f = e.father();
            le = f->level();
            while (le> level_)
            {
                f = f->father();
                le = f->level();
            }
            elemId = elementMapper_.map(*f);
        }
        else elemId = elementMapper_.map(e);
        Scalar permE = (*perm_)[elemId];

        for (int i = 0; i < dim; i++)
        permLoc_[i][i] = permE;

        return permLoc_;
    }

    void vtkout (const char* name, const Grid& grid_) const
    {
        Dune::VTKWriter<typename Grid::LeafGridView>
        vtkwriter(grid_.leafView());
        int size = (*perm_).size();
        vtkwriter.addCellData(*perm_, "absolute permeability");
        RepresentationType logPerm(size);
        for (int i = 0; i < size; i++)
        logPerm[i] = log10((*perm_)[i]);
        vtkwriter.addCellData(logPerm, "logarithm of permeability");
        vtkwriter.write(name, Dune::VTKOptions::ascii);
    }

    int level()
    {
        return level_;
    }

private:
    const Grid& grid_;
    const int level_;
    const char* fileName_;
    const bool createNew_;
    PermType perm_;
    Dune::FieldMatrix<Scalar,dim,dim> permLoc_;
    ElementMapper elementMapper_;
};
}

#endif

