// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup InputOutput
 * \brief Interface for passing data sets to a file and plotting them, if gnuplot
 *        is installed.
 *
 * The data sets for a specific window have to be passed by the addDataSet function
 * and then plotted by using the plot function.
 */
#ifndef DUMUX_GNUPLOT_INTERFACE_HH
#define DUMUX_GNUPLOT_INTERFACE_HH

#if !HAVE_GNUPLOT
// Gnuplot has not been found by CMake, no output possible.
#define GNUPLOT_EXECUTABLE "/usr/bin/gnuplot"
#endif

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/stdstreams.hh>

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief Interface for passing data sets to a file and plotting them, if gnuplot
 *        is installed.
 */
template<class Scalar>
class GnuplotInterface
{
public:
    using StringVector = std::vector<std::string>;
    enum class CurveType
    { function, file, data };
    using CurveTypeVector = std::vector<CurveType>;

    //! \brief The constructor
    explicit GnuplotInterface(bool persist = true) :
        pipeInteractive_(0), pipeImage_(0),
        openPlotWindow_(true), persist_(persist), createImage_(true),
        terminalType_("x11"), outputDirectory_("./"),
        datafileSeparator_(' '), linetype_("solid"),
        xRangeIsSet_(false), yRangeIsSet_(false),
        xLabel_(""), yLabel_(""),
        gnuplotPath_(GNUPLOT_EXECUTABLE)
    {
        open(persist_);
        resetPlot();
    }

    //! \brief The destructor
    ~GnuplotInterface()
    {
        close();
    }

    /*!
     * \brief Plots the files for a specific window number, writes a gnuplot and png file.
     *
     * \param filename The name of the output file
     */
    void plot(const std::string &filename = "")
    {
        // set correct terminal and general options
        std::string plot = "set datafile separator \'" + std::string(1, datafileSeparator_) + "\'\n";

        // set the labels and axes ranges
        plot += "set xlabel \"" + xLabel_ + "\"\n";
        plot += "set ylabel \"" + yLabel_ + "\"\n";
        if (xRangeIsSet_)
            plot += "set xrange [" + toStringWithPrecision(xRangeMin_) + ":" + toStringWithPrecision(xRangeMax_) + "]" + "\n";
        if (yRangeIsSet_)
            plot += "set yrange [" + toStringWithPrecision(yRangeMin_) + ":" + toStringWithPrecision(yRangeMax_) + "]" + "\n";

        // set user defined options
        plot += plotOptions_ + "\n";

        // plot curves
        plot += "plot";
        std::string plotCommandForFile(plot);
        for (unsigned int i = 0; i < curve_.size(); ++i)
        {
            if (curveType_[i] == CurveType::function)
            {
                plot += + " " + curve_[i] + " " + curveOptions_[i];
                plotCommandForFile += + " " + curve_[i] + " " + curveOptions_[i];
            }
            else
            {
                plot += + " '" + outputDirectory_ + curve_[i] + "' " + curveOptions_[i];
                plotCommandForFile += + " '" + curve_[i] + "' " + curveOptions_[i];
            }

            if (i < curve_.size()-1)
            {
                plot += ",\\";
                plotCommandForFile += ",\\";
            }
            plot += "\n";
            plotCommandForFile += "\n";
        }

        // initialize the interactive plot
        std::string interactivePlot = "reset\n";

        // set the terminal if the defaults were overwritten
        if (terminalType_.compare("x11") != 0 || linetype_.compare("solid") != 0)
            interactivePlot += "set term " + terminalType_ + " " + linetype_ + " " + " \n";

        // add the plot command and plot
        interactivePlot += plot;
        if (openPlotWindow_)
            executeGnuplot(interactivePlot, pipeInteractive_);

        // create a gnuplot file if a filename is specified
        if (filename.compare("") != 0)
        {
            std::string filePlot = "reset\n";
            filePlot += "set term pngcairo size 800,600 " + linetype_ + " \n";
            filePlot += "set output \"" + filename + ".png\"\n";
            filePlot += plot;
            std::string gnuplotFileName = outputDirectory_ + filename + ".gp";
            std::ofstream file;
            file.open(gnuplotFileName);
            file << filePlot;
            file.close();

            // create the image if desired
            if (createImage_)
                executeGnuplot(filePlot, pipeImage_);
        }
    }

    /*!
     * \brief Resets all gnuplot parameters
     */
    void resetAll(const bool persist = true)
    {
        close();
        open(persist);
        resetPlot();
    }

    /*!
     * \brief Deletes all plots from a plotting window and resets user-defined options
     */
    void resetPlot()
    {
        curve_.clear();
        curveOptions_.clear();
        plotOptions_ = "";
    }

    /*!
     * \brief Opens gnuplot
     */
    void open(const bool persist = true)
    {
        if (persist)
            pipeInteractive_ = popen((gnuplotPath_ + " -persist").c_str(), "w"); // "w" - writing
        else
            pipeInteractive_ = popen((gnuplotPath_).c_str(), "w");

        // the image pipe should not persist
        pipeImage_ = popen((gnuplotPath_).c_str(), "w");
    }

    /*!
     * \brief Closes gnuplot
     */
    void close()
    {
        if (pclose(pipeInteractive_) == -1 || pclose(pipeImage_) == -1)
            assert("Could not close pipe to Gnuplot!");
    }

    /*!
     * \brief Adds a function to list of plots
     *
     * \param function Function to be plotted
     * \param options Specific gnuplot options passed to this plot
     */
    void addFunctionToPlot(const std::string& function,
                           const std::string& options = "with lines")
    {
        curve_.push_back(function);
        curveOptions_.push_back(options);
        curveType_.push_back(CurveType::function);
    }

    /*!
     * \brief Adds a file to list of plots
     *
     * \param fileName Name and path of the file to be plotted
     * \param options Specific gnuplot options passed to this plot
     */
    void addFileToPlot(const std::string& fileName,
                       const std::string& options = "with lines")
    {
        curve_.push_back(fileName);
        curveOptions_.push_back(options);
        curveType_.push_back(CurveType::file);
    }

    /*!
     * \brief Adds a data set and writes a data file
     *
     * The title of the plot can be changed by setting the title in the options
     *
     * \param x Vector containing the x-axis data points
     * \param y Vector containing the y-axis data points
     * \param fileName The name of the written data file
     * \param options Specific gnuplot options passed to this plot
     */
    void addDataSetToPlot(const std::vector<Scalar>& x,
                          const std::vector<Scalar>& y,
                          const std::string& fileName,
                          const std::string& options = "with lines")
    {
        if (x.empty() || y.empty())
            DUNE_THROW(Dune::InvalidStateException, "Data vectors have to contain data!");

        if (x.size() > y.size())
            DUNE_THROW(Dune::InvalidStateException, "Non-matching data field sizes!");

        if (x.size() != y.size())
            std::cout << "GnuplotInterface warning: Added data fields of different size! "
                      << "Only plotting the first " << x.size() << " elements.\n";

        // write data to file
        std::ofstream file;
        file.open(outputDirectory_ + fileName);
        for (unsigned int i = 0; i < x.size(); i++)
        {
            checkNumber(x[i], "x[i] i=" + std::to_string(i) + " in " + fileName);
            checkNumber(y[i], "y[i] i=" + std::to_string(i) + " in " + fileName);
            file << x[i] << datafileSeparator_ << y[i] << std::endl;
        }
        file.close();

        // adding file to list of plotted lines
        curve_.push_back(fileName);
        curveOptions_.push_back(options);
        curveType_.push_back(CurveType::data);
    }

    /*!
     * \brief Sets the label for the x-axis
     *
     * \param label The label of the x-axis
     */
    void setXlabel(const std::string& label)
    {
        xLabel_ = label;
    }

    /*!
     * \brief Sets the label for the y-axis
     *
     * \param label The label of the y-axis
     */
    void setYlabel(const std::string& label)
    {
        yLabel_ = label;
    }

    /*!
     * \brief Sets the range for the x-axis
     *
     * \param min The lowest plotted value for the x-axis
     * \param max The highest plotted value for the x-axis
     */
    void setXRange(Scalar min, Scalar max)
    {
        xRangeMin_ = min;
        xRangeMax_ = max;
        xRangeIsSet_ = true;
    }

    /*!
     * \brief Sets the range for the y-axis
     *
     * \param min The lowest plotted value for the y-axis
     * \param max The highest plotted value for the y-axis
     */
    void setYRange(Scalar min, Scalar max)
    {
        yRangeMin_ = min;
        yRangeMax_ = max;
        yRangeIsSet_ = true;
    }

    /*!
     * \brief Sets additional user-defined options
     *
     * \param option Additional line of option in gnuplot language
     */
    void setOption(const std::string& option)
    {
        plotOptions_ += option + "\n";
    }

    /*!
     * \brief Define whether the gnuplot window should be opened
     *
     * \param openPlotWindow Open gnuplot or not
     */
    void setOpenPlotWindow(bool openPlotWindow)
    {
        openPlotWindow_ = openPlotWindow;
    }

    /*!
     * \brief Define whether gnuplot should create .png files
     *
     * \param createImage Create an image or not
     */
    void setCreateImage(bool createImage)
    {
        createImage_ = createImage;
    }

    /*!
     * \brief Sets the datafile separator
     *
     * \param separator The separator sign between two data columns
     */
    void setDatafileSeparator(char separator)
    {
        datafileSeparator_ = separator;
    }

    /*!
     * \brief Sets the terminal used for interactive outpu
     *
     * \param terminal The user-specified terminal
     */
    void setTerminalType(std::string terminal)
    {
        terminalType_ = terminal;
    }

    /*!
     * \brief Sets the output directory for data and gnuplot files
     *
     * \param outputDirectory The user-specified terminal
     */
    void setOutputDirectory(const std::string& outputDirectory)
    {
        outputDirectory_ = outputDirectory + "/";
    }

    /*!
     * \brief Use dashed (true) or solid (false) lines
     *
     * \param dashed Use dashed lines
     */
    void useDashedLines(bool dashed)
    {
        linetype_ = dashed ? "dashed" : "solid";
    }

private:
    // Give plot command to gnuplot
    void executeGnuplot(const std::string& plotCommand, std::FILE * pipe) const
    {
#ifdef HAVE_GNUPLOT
        fputs((plotCommand + "\n").c_str(), pipe);
        fflush(pipe);
#else
        std::cerr << "Warning: Gnuplot has not been found by CMake, no image generation or interactive display possible." << std::endl;
        std::cerr << "Note: The data and the gnuplot instruction file will still be created." << std::endl;
#endif
    }

    // Check validity of number
    void checkNumber(Scalar number, const std::string& text = "") const
    {
        using std::isnan;
        using std::isinf;
        if (isnan(number))
            Dune::dwarn << "warning: " << text << " is not a number, adjust your data range" << std::endl;
        if (isinf(number))
            Dune::dwarn << "warning: " << text << " is infinity, adjust your data range" << std::endl;
    }

    // Convert string with higher precision
    template <typename T>
    std::string toStringWithPrecision(const T value, const int n = 8)
    {
        std::ostringstream out;
        out << std::setprecision(n) << value;
        return out.str();
    }

    std::FILE * pipeInteractive_;
    std::FILE * pipeImage_;
    bool openPlotWindow_;
    bool persist_;
    bool createImage_;
    std::string terminalType_;
    std::string outputDirectory_;
    char datafileSeparator_;
    std::string linetype_;
    StringVector curve_;
    StringVector curveOptions_;
    CurveTypeVector curveType_;
    Scalar xRangeMin_;
    Scalar xRangeMax_;
    Scalar xRangeIsSet_;
    Scalar yRangeMin_;
    Scalar yRangeMax_;
    Scalar yRangeIsSet_;
    std::string xLabel_;
    std::string yLabel_;
    std::string plotOptions_;
    std::string gnuplotPath_;
};
} // end namespace Dumux
#endif
