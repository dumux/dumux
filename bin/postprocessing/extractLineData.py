import argparse
import csv
import sys
import os

# parse arguments
parser = argparse.ArgumentParser(
  prog='pvpython ' + sys.argv[0],
  description='Extract data from the paraview plotOverLine filter.'
)
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inputDirectory', help="Directory in which the .vtu are located")
parser.add_argument('-f', '--files', nargs='+', required=True, help="vtu files to be processed")
parser.add_argument('-o', '--outputDirectory', help="Directory to which the data files are written")
parser.add_argument('-p1', '--point1', type=float, nargs=3, required=True, help='Coordinates of the first point (in 3D)')
parser.add_argument('-p2', '--point2', type=float, nargs=3, required=True, help='Coordinates of the second point (in 3D)')
parser.add_argument('-r', '--resolution', type=int, default=1000, help='Resolution of the line (number of data points written to data file)')
parser.add_argument('-v', '--verbosity', type=int, default=2, help='Verbosity of the output. 1 = print progress. 2 = print data columns')
args = vars(parser.parse_args())

from paraview.simple import *

# import locations
inDirectory = args['inputDirectory']
outDirectory = args['outputDirectory']
if not inDirectory == '':
    inDirectory += '/'
if not outDirectory == '':
    outDirectory += '/'
    if not os.path.exists(outDirectory):
        os.makedirs(outDirectory)

# loop over all vtk files
counter = 1
for curFile in args['files']:
    # print progress to command line
    fileWithoutPath = os.path.basename(curFile)
    basename = os.path.splitext(fileWithoutPath)[0]
    if args['verbosity'] == 1:
        print("Processing file ({}/{}): {}".format(counter, len(args['files']), fileWithoutPath))
    counter += 1

    # load vtk file
    vtkFile = XMLUnstructuredGridReader(FileName=inDirectory+curFile)
    SetActiveSource(vtkFile)

    # apply and configure PlotOverLine filter
    plotOverLine = PlotOverLine(Source="High Resolution Line Source")
    plotOverLine.Source.Resolution = args['resolution']
    plotOverLine.Source.Point1 = args['point1']
    plotOverLine.Source.Point2 = args['point2']

    # write output to csv writer
    csvFile = outDirectory + basename + '.csv'
    writer = CreateWriter(csvFile, plotOverLine)
    writer.UpdatePipeline()

    # print the parameters and the column numbers
    if args['verbosity'] == 2:
        with open(csvFile) as f:
            print csvFile
            reader = csv.reader(f)
            paramList = list(reader)
            paramCounter=1
            for param in paramList[0]:
                print "%-2i   %s" % (paramCounter, param)
                paramCounter += 1
