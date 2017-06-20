#### TODO:
# - different colors for legend
# - read-in pvds with time outputs
# - read-in multiple vtus, e.g. for multidomain
# - rendering method 2d and 3d

# parse arguments
import argparse
import os
import sys

bool = ['True','False']
parameterType = ['CELLS','POINTS']
legendOrientation = ['Horizontal','Vertical']
parser = argparse.ArgumentParser(
  prog='\033[1m\033[94m' + 'pvbatch' + '\033[0m' + ' ' + sys.argv[0],
  description='Export a screenshot of a standard 2D plot. To change the color palette, change the default in the paraview GUI.'
)
# on/off-type features
offscreen = parser.add_mutually_exclusive_group(required=False)
offscreen.add_argument('--offscreen', dest='offscreen', action='store_true', help='Enable offscreen rendering (large pixel size)')
offscreen.add_argument('--no-offscreen', dest='offscreen', action='store_false', help='Disable offscreen rendering (low pixel sizes)')
parser.set_defaults(offscreen=False)
showaxesgrid = parser.add_mutually_exclusive_group(required=False)
showaxesgrid.add_argument('--showAxesGrid', dest='showAxesGrid', action='store_true', help='Show the axes grid for the domain')
showaxesgrid.add_argument('--no-showAxesGrid', dest='showAxesGrid', action='store_false', help='Do not show the axes grid for the domain')
parser.set_defaults(showAxesGrid=False)
showlegend = parser.add_mutually_exclusive_group(required=False)
showlegend.add_argument('--showLegend', dest='showLegend', action='store_true', help='Show the parameter legend/range')
showlegend.add_argument('--no-showLegend', dest='showLegend', action='store_false', help='Do not show the parameter legend/range')
parser.set_defaults(showLegend=True)
showorientaxes = parser.add_mutually_exclusive_group(required=False)
showorientaxes.add_argument('--showOrientationAxes', dest='showOrientationAxes', action='store_true', help='Show the orientation axis')
showorientaxes.add_argument('--no-showOrientationAxes', dest='showOrientationAxes', action='store_false', help='Do not the orientation axis')
parser.set_defaults(showOrientationAxes=False)
# more complicated features
parser.add_argument('-f', '--files', nargs='+', required=True, help="vtu files to be processed")
parser.add_argument('-o', '--outputDirectory', default='', help="Directory to which the pcitures are written")
parser.add_argument('-of', '--outFile', default='', help="Basename of the written png file")
parser.add_argument('-p', '--parameter', default='', help='The name of the parameter to be plotted')
parser.add_argument('-pc', '--parameterComponent', type=int, default=-1, help='Plot only a specific component of a vector (default: magnitude=-1)')
parser.add_argument('-pr', '--parameterRange', type=float, nargs=2, default=[0, 0], help='Adjustment of the color rane (default: min/max)')
parser.add_argument('-pt', '--parameterType', choices=parameterType, default='CELLS', help='The type of the data field (CELLS or POINTS)')
parser.add_argument('-lo', '--legendOrientation', choices=legendOrientation, default='Horizontal', help='The name of the parameter to be plotted')
parser.add_argument('-lp', '--legendPosition', type=float, nargs=2, default=[0.25, 0.0], help='The position of the legend')
parser.add_argument('-lz', '--legendZoom', type=float, nargs=2, default=[0.5, 0.5], help='The zoom of the legend')
parser.add_argument('-lt', '--legendTitle', default='', help="Title of the legend")
parser.add_argument('-lct', '--legendComponentTitle', default='none', help="Title of the legend component")
parser.add_argument('--size', type=int, nargs=2, default=[600, 400], help="The pixel size of the png file (default: 600x400)")
parser.add_argument('--scale', type=float, nargs=3, default=[1.0, 1.0, 1.0], help="Scaling factors in each direction")
parser.add_argument('--whiteBackground', dest='whiteBackground', action='store_true', help="Sets a white background")
parser.add_argument('-v', '--verbosity', type=int, default=2, help='Verbosity of the output')
args = vars(parser.parse_args())

try:
    from paraview.simple import *
except ImportError:
    print("`paraview.simple` not found. Make sure using pvbatch.")

# import locations
commonOutDirectory = False
outDirectory = args['outputDirectory']
if not outDirectory == '':
    outDirectory += '/'
    commonOutDirectory = True
    if not os.path.exists(outDirectory):
        os.makedirs(outDirectory)

# loop over all vtu files
counter = 1
for curFile in args['files']:
    # print progress to command line
    fileWithoutPath = os.path.basename(curFile)
    if not commonOutDirectory:
        abspath = os.path.abspath(curFile)
        outDirectory = os.path.dirname(abspath) + '/'
    basename = os.path.splitext(fileWithoutPath)[0]
    if args['verbosity'] == 1:
        print("Processing file ({}/{}): {}".format(counter, len(args['files']), fileWithoutPath))
    counter += 1

    # read vtu file and print available parameters
    vtuFile = XMLUnstructuredGridReader(FileName=curFile)
    if args['parameter'] == '':
        print "\nNo parameter was specified, use '-p PARAMETER' to specify it. Available parameters are:"
        if args['parameterType'] == 'CELLS':
            print vtuFile.CellArrayStatus
        else:
            print vtuFile.PointArrayStatus
        exit(1)

    # get active view
    renderView1 = GetActiveView()
    if not renderView1:
        # When using the ParaView UI, the View will be present, not otherwise.
        renderView1 = CreateRenderView()

    # print additional help message for large picture sizes
    if (args['size'][0] > 1024 or args['size'][1] > 1024) and args['offscreen'] == False:
        print "\nIt seems like you want to export a picture greater then your actual screen size. Use:"
        print "pvbatch --use-offscreen-rendering SCRIPT OPTIONS --offscreen"
        exit(2)
    renderView1.ViewSize = args['size']

    if args['showOrientationAxes'] == False:
        renderView1.OrientationAxesVisibility = 0

    if args['showAxesGrid'] == True:
        renderView1.AxesGrid.Visibility = 1

    # show data in view
    vtuFileDisplay = Show(vtuFile, renderView1)
    vtuFileDisplay.Scale = args['scale']

    # reset view to fit data
    renderView1.ResetCamera()

    # set scalar coloring
    ColorBy(vtuFileDisplay, (args['parameterType'], args['parameter']))

    # show color bar/color legend
    if args['showLegend'] == True:
        vtuFileDisplay.SetScalarBarVisibility(renderView1, True)

    # get color transfer function/color map for the parameter
    parameterLUT = GetColorTransferFunction(args['parameter'])

    # plot only a specific vector component
    if args['parameterComponent'] != -1:
        parameterLUT.VectorMode = 'Component'
        parameterLUT.VectorComponent = args['parameterComponent']
        #if args['parameterRange'][0] == 0 and args['parameterRange'][1] == 0:
            #vtuFileDisplay.RescaleTransferFunctionToDataRange(False)
        if args['showLegend'] == True:
            velocityLUTColorBar = GetScalarBar(parameterLUT, renderView1)
            velocityLUTColorBar.Title = args['parameter']
            velocityLUTColorBar.ComponentTitle = str(args['parameterComponent'])

    # adjust the range of the legend
    if args['parameterRange'][0] != 0 or args['parameterRange'][1] != 0:
        mean = (args['parameterRange'][0] + args['parameterRange'][1]) / 2
        parameterLUT.RGBPoints = [args['parameterRange'][0], 0.231373, 0.298039, 0.752941,
                                  mean, 0.865003, 0.865003, 0.865003,
                                  args['parameterRange'][1], 0.705882, 0.0156863, 0.14902]
    if args['parameterRange'][0] == 0 and args['parameterRange'][1] == 0:
        vtuFileDisplay.RescaleTransferFunctionToDataRange(False)

    # the placement and size of the legend
    legend = GetScalarBar(parameterLUT, renderView1)
    legend.Position = args['legendPosition']
    legend.Position2 = args['legendZoom']
    legend.Orientation = args['legendOrientation']

    # rename the legend if desired
    if args['legendTitle'] != '':
        legend.Title = args['legendTitle']
    if args['legendComponentTitle'] != 'none':
        legend.ComponentTitle = args['legendComponentTitle']

    # set a white background color and black color for fonts and the grid
    if args['whiteBackground'] == True:
        renderView1.Background = [255, 255, 255]
        legend.TitleColor = [0.0, 0.0, 0.0]
        legend.LabelColor = [0.0, 0.0, 0.0]
        if args['showAxesGrid'] == True:
            renderView1.AxesGrid.GridColor = [0.0, 0.0, 0.0]
            renderView1.AxesGrid.XTitleColor = [0.0, 0.0, 0.0]
            renderView1.AxesGrid.YTitleColor = [0.0, 0.0, 0.0]
            renderView1.AxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
            renderView1.AxesGrid.XLabelColor = [0.0, 0.0, 0.0]
            renderView1.AxesGrid.YLabelColor = [0.0, 0.0, 0.0]
            renderView1.AxesGrid.ZLabelColor = [0.0, 0.0, 0.0]

    # current camera placement for renderView1
    renderView1.InteractionMode = '2D'
    #renderView1.CameraPosition = [5.0, 0.12345, 0.0]
    #renderView1.CameraFocalPoint = [5.0, 0.12345, 0.0]

    # uncomment the following to render all views
    RenderAllViews()
    UseOffscreenRenderingForScreenshots = 1

    # save screenshot
    if not args['outFile'] == '':
        basename = args['outFile']
    SaveScreenshot(outDirectory + basename + '.png')
