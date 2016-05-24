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
 *
 * \brief Creating random fields using gstat
 */
#ifndef DUMUX_GSTAT_RANDOM_FIELD_HH
#define DUMUX_GSTAT_RANDOM_FIELD_HH

#if !HAVE_GSTAT
#warning Gstat has not been found by CMake, please set the variable GSTAT_EXECUTABLE manually.
// #define GSTAT_EXECUTABLE PATH
#endif

#include <iostream>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/common.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/istl/bvector.hh>

namespace Dumux
{

/*!
 * \brief Creating random fields using gstat
 *
 * gstat is an open source software tool which can (among other things) generate
 * geostatistical random fields (see <a href="www.gstat.org">www.gstat.org</a>).
 *
 * To use this class, unpack and install the zipped gstat tarball from the website
 * or use the script installexternal.sh provided with dumux.
 */
template<class GridView, class Scalar>
class GstatRandomField
{
    enum { dim = GridView::dimension,
           dimWorld = GridView::dimensionworld};

    typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > DataVector;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Traits::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IndexSet IndexSet;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, Dune::MCMGElementLayout> ElementMapper;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*!\brief Constructor.
     *
     * Creates a new field with random variables, if desired.
     * Otherwise creates a data field from already available data.
     *
     * For the a rand field generation three files are necessary.
     * A \a gstatControlFile in which all commands and in/output files for gstat are specified.
     * A \a gstatInputFile contains all coordinates (cell centers) of the grid, so that
     * gstat can perform its random realization. The filename must be same as in the gstatControlFile.
     * A \a gstatOutputFile in which gstat writes the random values to this file.
     * The filename must be same as in the gstatControlFile.
     *
     * \param gridView the used gridView
     * \param gstatControlFile name of control file for gstat
     * \param gstatInputFile name of input file for gstat
     * \param gstatOutputFile name of the gstat output file
     * \param createNew set true to create a new field
     */
    GstatRandomField(const GridView& gridView,
                     const std::string gstatControlFile,
                     const std::string gstatInputFile = "gstatInput.txt",
                     const std::string gstatOutputFile = "permeab.dat",
                     const bool createNew = true)
    : gridView_(gridView),
      data_(gridView.size(0)),
      powField_(false),
      elementMapper_(gridView)
    {
        ElementIterator eItEnd = gridView_.template end<0>();

#if HAVE_GSTAT
        if (createNew)
        {
            std::ofstream gstatInput(gstatInputFile);
            for (ElementIterator eIt = gridView_.template begin<0>(); eIt != eItEnd; ++eIt)
            {
                // get global coordinates of cell centers
                const GlobalPosition& globalPos = eIt->geometry().center();

                for (int dimIdx = 0; dimIdx < dim; dimIdx++)
                {
                    gstatInput << globalPos[dimIdx] << " " << std::flush;
                }
                gstatInput << std::endl;
            }
            gstatInput.close();

            std::string syscom;
            syscom = GSTAT_EXECUTABLE;
            syscom += " ";
            syscom += gstatControlFile;

            if (!gstatInput.good())
            {
                DUNE_THROW(Dune::IOError,
                           "Reading the gstat control file: " << gstatControlFile << " failed." << std::endl);
            }

            if (system(syscom.c_str()))
            {
                DUNE_THROW(Dune::IOError, "Executing gstat failed.");
            }
        }
#endif

        std::ifstream gstatOutput(gstatOutputFile);
        if (!gstatOutput.good())
        {
            DUNE_THROW(Dune::IOError,
                        "Reading from file: " << gstatOutputFile << " failed." << std::endl);
        }

        std::string line;
        std::getline(gstatOutput, line);

        Scalar dummy1, dummy2, dummy3, dataValue;
        for (ElementIterator eIt = gridView_.template begin<0>(); eIt != eItEnd; ++eIt)
        {
            std::getline(gstatOutput, line);
            std::istringstream curLine(line);
            if (dim == 1)
                curLine >> dummy1 >> dataValue;
            else if (dim == 2)
                curLine >> dummy1 >> dummy2 >> dataValue;
            else // dim = 3
                curLine >> dummy1 >> dummy2 >> dummy3 >> dataValue;
            data_[elementMapper_.index(*eIt)] = dataValue;
        }
        gstatOutput.close();
    }

    //! \brief Use the input as an exponent
    void createPowField()
    {
        for (unsigned int i = 0; i < data_.size(); ++i)
        {
            data_[i] = std::pow(10.0, data_[i]);
        }
        powField_ = true;
    }

    //! \brief Return an entry of the data vector
    Scalar data(const Element& e)
    {
        return data_[elementMapper_.index(e)];
    }

    //! \brief Write the data to a vtk file
    void writeVtk(const std::string vtkName,
                  const std::string dataName = "data") const
    {
        Dune::VTKWriter<GridView> vtkwriter(gridView_);
        vtkwriter.addCellData(data_, dataName);

        DataVector logPerm(data_.size());
        if (powField_)
        {
            for (unsigned int i = 0; i < data_.size(); ++i)
            {
                logPerm[i] = std::log10(data_[i]);
            }
            vtkwriter.addCellData(logPerm, "logarithm of " + dataName);
        }
        vtkwriter.write(vtkName, Dune::VTK::OutputType::ascii);
    }

private:
    const GridView& gridView_;
    DataVector data_;
    bool powField_;
    ElementMapper elementMapper_;
};

}

#endif
