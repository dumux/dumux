// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*!
 * \file
 * \ingroup InputOutput
 * \brief Base class to write pvd-files which contains a list of all collected vtk-files.
 *         This is a modified version of DUNE's pvd writer which takes a VTKWriter as template
 *         argument making it more general.
 */
#ifndef DUMUX_VTKSEQUENCEWRITER_HH
#define DUMUX_VTKSEQUENCEWRITER_HH

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <memory>

#include <dune/grid/io/file/vtk/common.hh>
#include <dune/common/path.hh>


namespace Dumux {

  /*!
   * \ingroup InputOutput
   * \brief Base class to write pvd-files which contains a list of all collected vtk-files.
   *         This is a modified version of DUNE's pvd writer which takes a VTKWriter as template
   *         argument making it more general.
   *
   * Derive from this class to write pvd-file suitable for easy visualization with
   * <a href="http://www.vtk.org/">The Visualization Toolkit (VTK)</a>.
   *
   * \tparam VTKWriter The VTKWriter class
   *
   */
  template<class VTKWriter>
  class VTKSequenceWriter
  {
    std::shared_ptr<VTKWriter > vtkWriter_;
    std::vector<double> timesteps_;
    std::string name_,path_,extendpath_;
    int rank_;
    int size_;
  public:
    /*! \brief Set up the VTKSequenceWriter class
     *
     * \param 				vtkWriter Writer object used to write the individual time step data files
     * \param name 			Base name of the output files.  This should not
     *                   	contain any directory part and not filename
     *                   	extensions.  It will be used both for each processes
     *                   	piece as well as the parallel collection file.
     * \param path 			Directory where to put the parallel collection
     *                   	(.pvtu/.pvtp) file.  If it is relative, it is taken
     *                   	relative to the current directory
     * \param extendpath 	Directory where to put the piece file (.vtu/.vtp) of
     *                   	this process.  If it is relative, it is taken
     *                   	relative to the directory denoted by path
     * \param rank 			Process number in a multi-process setting
     * \param size 			Total number of processes
     */
    explicit VTKSequenceWriter( std::shared_ptr<VTKWriter > vtkWriter,
                                    const std::string& name,
                                    const std::string& path,
                                    const std::string& extendpath,
                                    int rank,
                                    int size)
      : vtkWriter_(vtkWriter),
        name_(name), path_(path),
        extendpath_(extendpath),
        rank_(rank),
        size_(size)
    {}

    ~VTKSequenceWriter() {}


    /*!
     * \brief Writes VTK data for the given time,
     * \param time The time(step) for the data to be written.
     * \param type VTK output type.
     */
    void write (double time, Dune::VTK::OutputType type = Dune::VTK::ascii)
    {
      /* remember current time step */
      unsigned int count = timesteps_.size();
      timesteps_.push_back(time);

      /* write VTK file */
      if(size_==1)
        vtkWriter_->write(Dune::concatPaths(path_,seqName(count)),type);
      else
        vtkWriter_->pwrite(seqName(count), path_,extendpath_,type);

      /* write pvd file ... only on rank 0 */
      if (rank_==0) {
        std::ofstream pvdFile;
        pvdFile.exceptions(std::ios_base::badbit | std::ios_base::failbit |
                           std::ios_base::eofbit);
        std::string pvdname = name_ + ".pvd";
        pvdFile.open(pvdname.c_str());
        pvdFile << "<?xml version=\"1.0\"?> \n"
                << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"" << Dune::VTK::getEndiannessString() << "\"> \n"
                << "<Collection> \n";
        for (unsigned int i=0; i<=count; i++)
        {
          // filename
          std::string piecepath;
          std::string fullname;
          if(size_==1) {
            piecepath = path_;
            fullname = vtkWriter_->getSerialPieceName(seqName(i), piecepath);
          }
          else {
            piecepath = Dune::concatPaths(path_, extendpath_);
            fullname = vtkWriter_->getParallelHeaderName(seqName(i), piecepath, size_);
          }
          pvdFile << "<DataSet timestep=\"" << timesteps_[i]
                  << "\" group=\"\" part=\"0\" name=\"\" file=\""
                  << fullname << "\"/> \n";
        }
        pvdFile << "</Collection> \n"
                << "</VTKFile> \n" << std::flush;
        pvdFile.close();
      }
    }
  private:

    // create sequence name
    std::string seqName(unsigned int count) const
    {
      std::stringstream n;
      n.fill('0');
      n << name_ << "-" << std::setw(5) << count;
      return n.str();
    }
  };

} // end namespace Dumux

#endif
