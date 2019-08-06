# import visit_utils, we will use it to help encode our movie
import sys
sys.path.append("/soft/visualization/visit/2.13.2/linux-x86_64/lib/site-packages")

from visit_utils import *
from visit import *

LaunchNowin()

OpenDatabase("movie.visit")

AddPlot("Pseudocolor", "phi", 1, 1)
AddPlot("Subset", "patches", 1, 1)
SetActivePlots(1)
SetActivePlots(1)
SubsetAtts = SubsetAttributes()
SubsetAtts.colorType = SubsetAtts.ColorByMultipleColors  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
SubsetAtts.colorTableName = "Default"
SubsetAtts.invertColorTable = 0
SubsetAtts.legendFlag = 1
SubsetAtts.lineStyle = SubsetAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
SubsetAtts.lineWidth = 0
SubsetAtts.singleColor = (0, 0, 0, 255)
SubsetAtts.SetMultiColor(0, (255, 0, 0, 255))
SubsetAtts.SetMultiColor(1, (0, 255, 0, 255))
SubsetAtts.SetMultiColor(2, (0, 0, 255, 255))
SubsetAtts.SetMultiColor(3, (0, 255, 255, 255))
SubsetAtts.SetMultiColor(4, (255, 0, 255, 255))
SubsetAtts.SetMultiColor(5, (255, 255, 0, 255))
SubsetAtts.SetMultiColor(6, (255, 135, 0, 255))
SubsetAtts.SetMultiColor(7, (255, 0, 135, 255))
SubsetAtts.SetMultiColor(8, (168, 168, 168, 255))
SubsetAtts.SetMultiColor(9, (255, 68, 68, 255))
SubsetAtts.SetMultiColor(10, (99, 255, 99, 255))
SubsetAtts.SetMultiColor(11, (99, 99, 255, 255))
SubsetAtts.SetMultiColor(12, (40, 165, 165, 255))
SubsetAtts.SetMultiColor(13, (255, 99, 255, 255))
SubsetAtts.SetMultiColor(14, (255, 255, 99, 255))
SubsetAtts.SetMultiColor(15, (255, 170, 99, 255))
SubsetAtts.SetMultiColor(16, (170, 79, 255, 255))
SubsetAtts.subsetNames = ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17")
SubsetAtts.opacity = 1
SubsetAtts.wireframe = 1
SubsetAtts.drawInternal = 0
SubsetAtts.smoothingLevel = 0
SubsetAtts.pointSize = 0.05
SubsetAtts.pointType = SubsetAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
SubsetAtts.pointSizeVarEnabled = 0
SubsetAtts.pointSizeVar = "default"
SubsetAtts.pointSizePixels = 2
SetPlotOptions(SubsetAtts)
DrawPlots()

# get the number of timesteps
nts = TimeSliderGetNStates()

# set basic save options
swatts = SaveWindowAttributes()
#
# The 'family' option controls if visit automatically adds a frame number to 
# the rendered files. For this example we will explicitly manage the output name.
#
swatts.family = 0
#
# select PNG as the output file format
#
swatts.format = swatts.PNG 
#
# set the width of the output image
#
swatts.width = 1024 
#
# set the height of the output image
#
swatts.height = 1024

#the encoder expects file names with an integer sequence
# 0,1,2,3 .... N-1

file_idx = 0

for ts in range(0,nts,1): # look at every frame
    # Change to the next timestep
    TimeSliderSetState(ts)
    #before we render the result, explicitly set the filename for this render
    swatts.fileName = "amr_advection_%04d.png" % file_idx
    SetSaveWindowAttributes(swatts)
    # render the image to a PNG file
    SaveWindow()
    file_idx +=1

    ################
    # use visit_utils.encoding to encode these images into a "wmv" movie
    #
    # The encoder looks for a printf style pattern in the input path to identify the frames of the movie. 
    # The frame numbers need to start at 0. 
    # 
    # The encoder selects a set of decent encoding settings based on the extension of the
    # the output movie file (second argument). In this case we will create a "wmv" file. 
    # 
    # Other supported options include ".mpg", ".mov". 
    #   "wmv" is usually the best choice and plays on all most all platforms (Linux ,OSX, Windows). 
    #   "mpg" is lower quality, but should play on any platform.
    #
    # 'fdup' controls the number of times each frame is duplicated. 
    #  Duplicating the frames allows you to slow the pace of the movie to something reasonable. 
    #
    ################

    input_pattern = "amr_advection_%04d.png"
    output_movie = "amr_advection.mp4"
    encoding.encode(input_pattern,output_movie,fdup=4)
