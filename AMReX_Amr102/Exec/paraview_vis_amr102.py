# trace generated using paraview version 5.8.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *

import glob

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get all the plotfiles
AllPlotFiles = sorted(glob.glob("plt" + "[0-9]"*5))

# create a new 'AMReX/BoxLib Grid Reader'
plt00 = AMReXBoxLibGridReader(FileNames=AllPlotFiles)
plt00.CellArrayStatus = []

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Properties modified on plt00
plt00.CellArrayStatus = ['phi', 'proc', 'vfrac', 'xvel', 'yvel', 'zvel']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1166, 1176]

# get layout
layout1 = GetLayout()

# show data in view
plt00Display = Show(plt00, renderView1, 'AMRRepresentation')

# trace defaults for the display properties.
plt00Display.Representation = 'Outline'
plt00Display.ColorArrayName = [None, '']
plt00Display.OSPRayScaleFunction = 'PiecewiseFunction'
plt00Display.SelectOrientationVectors = 'None'
plt00Display.ScaleFactor = 0.1
plt00Display.SelectScaleArray = 'None'
plt00Display.GlyphType = 'Arrow'
plt00Display.GlyphTableIndexArray = 'None'
plt00Display.GaussianRadius = 0.005
plt00Display.SetScaleArray = [None, '']
plt00Display.ScaleTransferFunction = 'PiecewiseFunction'
plt00Display.OpacityArray = [None, '']
plt00Display.OpacityTransferFunction = 'PiecewiseFunction'
plt00Display.DataAxesGrid = 'GridAxesRepresentation'
plt00Display.PolarAxes = 'PolarAxesRepresentation'
plt00Display.ScalarOpacityUnitDistance = 0.04436647145156464

# reset view to fit data
renderView1.ResetCamera()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Slice'
slice1 = Slice(Input=plt00)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [0.5, 0.5, 0.0625]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [0.5, 0.5, 0.0625]

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=slice1.SliceType)

# Properties modified on slice1.SliceType
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

# show data in view
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = [None, '']
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 0.1
slice1Display.SelectScaleArray = 'None'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'None'
slice1Display.GaussianRadius = 0.005
slice1Display.SetScaleArray = [None, '']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = [None, '']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.PolarAxes = 'PolarAxesRepresentation'

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(slice1Display, ('FIELD', 'vtkBlockColors'))

# get color transfer function/color map for 'vtkBlockColors'
vtkBlockColorsLUT = GetColorTransferFunction('vtkBlockColors')

# get opacity transfer function/opacity map for 'vtkBlockColors'
vtkBlockColorsPWF = GetOpacityTransferFunction('vtkBlockColors')

# set scalar coloring
ColorBy(slice1Display, ('CELLS', 'phi'))

# rescale color and/or opacity maps used to include current data range
slice1Display.RescaleTransferFunctionToDataRange(True, False)

# get color transfer function/color map for 'phi'
phiLUT = GetColorTransferFunction('phi')

# get opacity transfer function/opacity map for 'phi'
phiPWF = GetOpacityTransferFunction('phi')

# create a new 'AMReX/BoxLib Particles Reader'
plt00_1 = AMReXBoxLibParticlesReader(FileNames=AllPlotFiles)
plt00_1.PointArrayStatus = ['id', 'cpu', 'real_comp0', 'real_comp1', 'real_comp2', 'real_comp3']

# show data in view
plt00_1Display = Show(plt00_1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
plt00_1Display.Representation = 'Surface'
plt00_1Display.ColorArrayName = [None, '']
plt00_1Display.OSPRayScaleArray = 'cpu'
plt00_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
plt00_1Display.SelectOrientationVectors = 'None'
plt00_1Display.ScaleFactor = 0.029962979569768036
plt00_1Display.SelectScaleArray = 'None'
plt00_1Display.GlyphType = 'Arrow'
plt00_1Display.GlyphTableIndexArray = 'None'
plt00_1Display.GaussianRadius = 0.0014981489784884018
plt00_1Display.SetScaleArray = ['POINTS', 'cpu']
plt00_1Display.ScaleTransferFunction = 'PiecewiseFunction'
plt00_1Display.OpacityArray = ['POINTS', 'cpu']
plt00_1Display.OpacityTransferFunction = 'PiecewiseFunction'
plt00_1Display.DataAxesGrid = 'GridAxesRepresentation'
plt00_1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
plt00_1Display.ScaleTransferFunction.Points = [2.0, 0.0, 0.5, 0.0, 3.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
plt00_1Display.OpacityTransferFunction.Points = [2.0, 0.0, 0.5, 0.0, 3.0, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(plt00_1Display, ('FIELD', 'vtkBlockColors'))

# create a new 'Glyph'
glyph1 = Glyph(Input=plt00_1,
    GlyphType='Arrow')
glyph1.OrientationArray = ['POINTS', 'No orientation array']
glyph1.ScaleArray = ['POINTS', 'No scale array']
glyph1.ScaleFactor = 0.029962979569768036
glyph1.GlyphTransform = 'Transform2'

# Properties modified on glyph1
glyph1.GlyphType = 'Sphere'
glyph1.ScaleFactor = 0.01
glyph1.GlyphMode = 'Every Nth Point'
glyph1.Stride = 100

# show data in view
glyph1Display = Show(glyph1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.ColorArrayName = [None, '']
glyph1Display.OSPRayScaleArray = 'Normals'
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.SelectOrientationVectors = 'None'
glyph1Display.ScaleFactor = 0.030791985988616946
glyph1Display.SelectScaleArray = 'None'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GlyphTableIndexArray = 'None'
glyph1Display.GaussianRadius = 0.0015395992994308473
glyph1Display.SetScaleArray = ['POINTS', 'Normals']
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityArray = ['POINTS', 'Normals']
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
glyph1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
glyph1Display.ScaleTransferFunction.Points = [-0.9749279618263245, 0.0, 0.5, 0.0, 0.9749279618263245, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
glyph1Display.OpacityTransferFunction.Points = [-0.9749279618263245, 0.0, 0.5, 0.0, 0.9749279618263245, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(glyph1Display, ('FIELD', 'vtkBlockColors'))

# create a new 'XML Partitioned Polydata Reader'
ebpvtp = XMLPartitionedPolydataReader(FileName=['eb.pvtp'])

# show data in view
ebpvtpDisplay = Show(ebpvtp, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
ebpvtpDisplay.Representation = 'Surface'
ebpvtpDisplay.ColorArrayName = [None, '']
ebpvtpDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
ebpvtpDisplay.SelectOrientationVectors = 'None'
ebpvtpDisplay.ScaleFactor = 0.019999998807907107
ebpvtpDisplay.SelectScaleArray = 'None'
ebpvtpDisplay.GlyphType = 'Arrow'
ebpvtpDisplay.GlyphTableIndexArray = 'None'
ebpvtpDisplay.GaussianRadius = 0.0009999999403953552
ebpvtpDisplay.SetScaleArray = [None, '']
ebpvtpDisplay.ScaleTransferFunction = 'PiecewiseFunction'
ebpvtpDisplay.OpacityArray = [None, '']
ebpvtpDisplay.OpacityTransferFunction = 'PiecewiseFunction'
ebpvtpDisplay.DataAxesGrid = 'GridAxesRepresentation'
ebpvtpDisplay.PolarAxes = 'PolarAxesRepresentation'

# update the view to ensure updated data information
renderView1.Update()

# current camera placement for renderView1
renderView1.CameraPosition = [0.806487938976244, 0.6098478265765159, 2.3240231832552687]
renderView1.CameraFocalPoint = [0.4999999999999996, 0.49999999999999956, 0.06249999999999997]
renderView1.CameraViewUp = [-0.004859468445493335, 0.9988423431014976, -0.04785769733216954]
renderView1.CameraParallelScale = 0.7159515667518355

# save animation
SaveAnimation('amr102.avi', renderView1, ImageResolution=[1024, 1024],
    FrameRate=35,
    FrameWindow=[0, len(AllPlotFiles)-1])

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [0.806487938976244, 0.6098478265765159, 2.3240231832552687]
renderView1.CameraFocalPoint = [0.4999999999999996, 0.49999999999999956, 0.06249999999999997]
renderView1.CameraViewUp = [-0.004859468445493335, 0.9988423431014976, -0.04785769733216954]
renderView1.CameraParallelScale = 0.7159515667518355

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
