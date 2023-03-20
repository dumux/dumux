import argparse
import csv
import sys
import os

# parse arguments
parser = argparse.ArgumentParser(
  prog='\033[1m\033[94m' + 'pvpython' + '\033[0m' + ' ' + sys.argv[0],
  description='Extract data from the paraview plotOverLine filter.'
)
parser.add_argument('-f', '--files', nargs='+', required=True, help="vtu files to be processed")
parser.add_argument('-o', '--outputDirectory', default='', help="Directory to which the data files are written")
parser.add_argument('-of', '--outFile', default='', help="Basename of the written csv file")
parser.add_argument('-p1', '--point1', type=float, nargs=3, required=True, help='Coordinates of the first point (in 3D)')
parser.add_argument('-p2', '--point2', type=float, nargs=3, required=True, help='Coordinates of the second point (in 3D)')
parser.add_argument('-r', '--resolution', type=int, default=1000, help='Resolution of the line (number of data points written to data file)')
parser.add_argument('-v', '--verbosity', type=int, default=2, help='Verbosity of the output. 1 = print progress. 2 = print data columns')
args = vars(parser.parse_args())

try:
    from paraview.simple import *
except:
    raise ImportError("`paraview.simple` not found. Make sure using pvpython instead of python.")

# import locations
outDirectory = args['outputDirectory']
if outDirectory.strip():
    os.makedirs(outDirectory, exist_ok=True)

# loop over all vtk files
counter = 0
for curFile in args['files']:

    # if no output directory was specified, use the directory of the given file
    curOutDirectory = outDirectory
    if curOutDirectory == '':
        curOutDirectory = os.path.dirname( os.path.abspath(curFile) )

    # if no basename was specified, reuse the file name for the .csv file
    if args['outFile'] == '':
        csvFileName = os.path.join(curOutDirectory, os.path.splitext(os.path.basename(curFile))[0] + ".csv")
    elif len(args['files']) > 1:
        csvFileName = os.path.join(curOutDirectory, args['outFile'] + "_" + str(counter) + ".csv")
    else:
        csvFileName = os.path.join(curOutDirectory, args['outFile'] + ".csv")
    counter += 1

    # print progress to command line
    if args['verbosity'] == 1:
        print("Processing file ({}/{}): {}".format(counter, len(args['files']), os.path.basename(curFile)))

    # load vtk file
    if os.path.splitext(curFile)[1] == ".vtp":
        vtkFile = XMLPolyDataReader(FileName=curFile)
    else:
        vtkFile = XMLUnstructuredGridReader(FileName=curFile)
    SetActiveSource(vtkFile)

    # apply and configure PlotOverLine filter
    plotOverLine = PlotOverLine(Source="High Resolution Line Source")
    plotOverLine.Source.Resolution = args['resolution']
    plotOverLine.Source.Point1 = args['point1']
    plotOverLine.Source.Point2 = args['point2']

    # write output to csv writer
    writer = CreateWriter(csvFileName, plotOverLine)
    writer.UpdatePipeline()

    # print the parameters and the column numbers
    if args['verbosity'] == 2:
        with open(csvFileName) as csvFile:
            print(csvFileName)
            paramList = list(csv.reader(csvFile))[0]
            for i, param in enumerate(paramList):
                print("{:>5}  {}".format(i+1, param))
