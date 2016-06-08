import argparse
import csv
import fileinput
import os
import sys

# parse arguments
parser = argparse.ArgumentParser(
  prog='pvpython ' + sys.argv[0],
  description='Extract data from the paraview probeLocation and plotOverTime filters.'
)
parser.add_argument('-i', '--inputDirectory', help="Directory in which the pvd file is located")
parser.add_argument('-f', '--files', nargs='+', required=True, help="pvd files to be processed")
parser.add_argument('-o', '--outputDirectory', help="Directory to which the .csv files are written")
parser.add_argument('-p', '--point', type=float, nargs=3, required=True, help='Coordinates of the probed point (in 3D)')
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

# loop over all pvd files
counter = 1
for curFile in args['files']:
    # print progress to command line
    fileWithoutPath = os.path.basename(curFile)
    basename = os.path.splitext(fileWithoutPath)[0]
    if args['verbosity'] == 1:
        print("Processing file ({}/{}): {}".format(counter, len(args['files']), inDirectory+curFile))
    counter += 1

    # load pvd file
    pvdFile = PVDReader(FileName=inDirectory+curFile)

    # Extract Point and Probe at a location
    selectionSource = IDSelectionSource( ContainingCells=0, InsideOut=False, FieldType='POINT', IDs=0 )
    ExtractSelection = ExtractSelection(Selection=selectionSource, Input=pvdFile)
    ExtractSelection.UpdatePipeline()
    selectionData = servermanager.Fetch(ExtractSelection)

    probeLocation = ProbeLocation()
    probeLocation.Input = pvdFile
    pointSource = probeLocation.ProbeType
    pointSource.Center.SetData(args['point'])
    probeLocation.UpdatePipeline()

    # Parse the extracted source and plot over time
    selectionSourceprobeLocationation = IDSelectionSource( ContainingCells=0, InsideOut=False, FieldType='POINT', IDs=[0, 0] )
    plotSelectionOverTime = PlotSelectionOverTime(OnlyReportSelectionStatistics=False)
    plotSelectionOverTime.Selection = selectionSourceprobeLocationation
    plotSelectionOverTime.Input = probeLocation

    # write output to csv writer
    csvFile = outDirectory + basename + '.csv'
    writer = CreateWriter(csvFile, plotSelectionOverTime)
    writer.UpdatePipeline()

    # print the parameters and the column numbers
    if args['verbosity'] == 2:
        with open(outDirectory + basename + '0.csv') as file:
            print outDirectory + basename + '.csv'
            reader = csv.reader(file)
            paramList = list(reader)
            paramCounter=1
            for param in paramList[0]:
                print "%-2i   %s" % (paramCounter, param)
                paramCounter += 1

    # create a proper csv file with semicolons as separators
    f = open(outDirectory + basename + '0.csv', 'r')
    filedata = f.read()
    f.close()
    os.remove(outDirectory + basename + '0.csv')
    newdata = filedata.replace(',', ';')
    f = open(csvFile,'w')
    f.write(newdata)
    f.close()
