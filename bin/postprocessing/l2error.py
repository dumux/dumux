import argparse
import csv
import sys

#Auxiliary function that provides a handy parser
parser = argparse.ArgumentParser(prog='python ' + sys.argv[0], description='Calculate the l2 error of csv data files.')
parser.add_argument('-f1', '--reference', type=str, required=True, help='Reference csv-file')
parser.add_argument('-f2', '--newSolution', type=str, required=True, help='NewSolution csv-file')
parser.add_argument('-xMin', '--xMin', type=float, required=False, default=-1e99, help='Restrict data to x>xMin')
parser.add_argument('-xMax', '--xMax', type=float, required=False, default=1e99, help='Restrict data to x>xMax')
parser.add_argument('-x1', '--xData1', type=int, required=True, help='Column index of x data in reference')
parser.add_argument('-x2', '--xData2', type=int, required=True, help='Column index of x data in reference')
parser.add_argument('-y1', '--yData1', type=int, required=True, help='Column index of y data in newSolution')
parser.add_argument('-y2', '--yData2', type=int, required=True, help='Column index of y data in newSolution')
parser.add_argument('-p', '--percent', action='store_true', help='Print errors in percent')
parser.add_argument('-f', '--force', action='store_true', help='Ignore \'not-matching\' errors')
parser.add_argument('-v', '--verbose', action='store_true', help='Verbosity of the script')
args = vars(parser.parse_args())

with open(args['reference'], 'rb') as referenceFile:
  reader = csv.reader(referenceFile)
  reference = list(reader)
with open(args['newSolution'], 'rb') as newSolutionFile:
  reader = csv.reader(newSolutionFile)
  newSolution = list(reader)

if (reference[0][args['xData1']] != newSolution[0][args['xData2']] and not args['force']):
    print "X-Identifier not equal: ref=", reference[0][args['xData1']], ",new=", newSolution[0][args['xData2']], ". Aborting! (Use -f to continue anyway)"
    exit (1)

if (reference[0][args['yData1']] != newSolution[0][args['yData2']] and not args['force']):
    print "Y-Identifier not equal. ref=", reference[0][args['yData1']], ",new=", newSolution[0][args['yData2']], ". Aborting! (Use -f to continue anyway)"
    exit (2)

if (len(reference) != len(newSolution)):
    print "Length of reference and newSolution not equal: ref=", len(reference), ",new=", len(newSolution), ". Aborting!"
    exit (3)

distanceOld = 0.0
sumError = 0.0
sumReference = 0.0
sumDistance = 0.0
numPoints = 0
for i in range(1,len(reference)):
    coord_ref = float(reference[i][args['xData1']])
    coord_newSolution = float(newSolution[i][args['xData2']])
    if (coord_ref != coord_newSolution):
        print "Coordinates not equal: ref=", coord_ref, ",new=", coord_newSolution, ". Aborting!"
        exit (4)
    if (coord_ref < float(args['xMin']) or coord_ref > float(args['xMax'])):
        continue

    if (i == 1):
        distance = 0.5*(float(reference[2][args['xData1']]) - float(reference[1][args['xData1']]))
    elif (i == len(reference)-1):
        distance = 0.5*(float(reference[len(reference)-1][args['xData1']]) - float(reference[len(reference)-2][args['xData1']]))
    else:
        distance = 0.5*(float(reference[i+1][args['xData1']]) - float(reference[i-1][args['xData1']]))
    sumError += ((float(reference[i][args['yData1']])-float(newSolution[i][args['yData2']]))**2)*distance
    sumReference += ((float(reference[i][args['yData1']]))**2)*distance
    sumDistance += distance
    numPoints += 1

if (numPoints < 999 and not args['force']):
    print "Warning: numPoints=", numPoints, " is low, could result in bad the error approximation. (Use -f to suppress this warning)"

l2normAbs = (sumError/sumDistance)**0.5 # numPoints is needed, resulting from the equidistant integration
l2normRel = (sumError/sumReference)**0.5 # numPoints cancels out for equidistant integration

if (args['percent']):
    print "L2_Error_in_%: ", "{0:.5f}%".format(l2normAbs*100), "Rel_L2_Error_in_%: ", "{0:.5f}%".format(l2normRel*100)
else:
    print "L2_Error: ", "{0:.5e}".format(l2normAbs), " Rel_L2_Error: ", "{0:.5e}".format(l2normRel)
