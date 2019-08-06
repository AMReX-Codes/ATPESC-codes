# trace generated using paraview version 5.6.1
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
import glob

AllPlotFiles = sorted(glob.glob("plt[0-9][0-9][0-9][0-9][0-9]"))

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Partitioned Polydata Reader'
ebpvtp = XMLPartitionedPolydataReader(FileName=['eb.pvtp'])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
renderView1.ViewSize = [1400, 801]

# show data in view
ebpvtpDisplay = Show(ebpvtp, renderView1)

# trace defaults for the display properties.
ebpvtpDisplay.Representation = 'Surface'
ebpvtpDisplay.ColorArrayName = [None, '']
ebpvtpDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
ebpvtpDisplay.SelectOrientationVectors = 'None'
ebpvtpDisplay.ScaleFactor = 0.14000000208616256
ebpvtpDisplay.SelectScaleArray = 'None'
ebpvtpDisplay.GlyphType = 'Arrow'
ebpvtpDisplay.GlyphTableIndexArray = 'None'
ebpvtpDisplay.GaussianRadius = 0.007000000104308128
ebpvtpDisplay.SetScaleArray = [None, '']
ebpvtpDisplay.ScaleTransferFunction = 'PiecewiseFunction'
ebpvtpDisplay.OpacityArray = [None, '']
ebpvtpDisplay.OpacityTransferFunction = 'PiecewiseFunction'
ebpvtpDisplay.DataAxesGrid = 'GridAxesRepresentation'
ebpvtpDisplay.SelectionCellLabelFontFile = ''
ebpvtpDisplay.SelectionPointLabelFontFile = ''
ebpvtpDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
ebpvtpDisplay.DataAxesGrid.XTitleFontFile = ''
ebpvtpDisplay.DataAxesGrid.YTitleFontFile = ''
ebpvtpDisplay.DataAxesGrid.ZTitleFontFile = ''
ebpvtpDisplay.DataAxesGrid.XLabelFontFile = ''
ebpvtpDisplay.DataAxesGrid.YLabelFontFile = ''
ebpvtpDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
ebpvtpDisplay.PolarAxes.PolarAxisTitleFontFile = ''
ebpvtpDisplay.PolarAxes.PolarAxisLabelFontFile = ''
ebpvtpDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
ebpvtpDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# reset view to fit data
renderView1.ResetCamera()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# create a new 'VisItBoxlib3DReader'
header = VisItBoxlib3DReader(FileName=['plt00000/Header'])
header.MeshStatus = ['Mesh']
header.MaterialStatus = []
header.CellArrayStatus = []

# Properties modified on header
header.CellArrayStatus = ['proc']

# show data in view
headerDisplay = Show(header, renderView1)

# trace defaults for the display properties.
headerDisplay.Representation = 'Outline'
headerDisplay.ColorArrayName = [None, '']
headerDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
headerDisplay.SelectOrientationVectors = 'None'
headerDisplay.ScaleFactor = 0.1999999761581421
headerDisplay.SelectScaleArray = 'None'
headerDisplay.GlyphType = 'Arrow'
headerDisplay.GlyphTableIndexArray = 'None'
headerDisplay.GaussianRadius = 0.009999998807907105
headerDisplay.SetScaleArray = [None, '']
headerDisplay.ScaleTransferFunction = 'PiecewiseFunction'
headerDisplay.OpacityArray = [None, '']
headerDisplay.OpacityTransferFunction = 'PiecewiseFunction'
headerDisplay.DataAxesGrid = 'GridAxesRepresentation'
headerDisplay.SelectionCellLabelFontFile = ''
headerDisplay.SelectionPointLabelFontFile = ''
headerDisplay.PolarAxes = 'PolarAxesRepresentation'
headerDisplay.ScalarOpacityUnitDistance = 0.05194153242079811

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
headerDisplay.DataAxesGrid.XTitleFontFile = ''
headerDisplay.DataAxesGrid.YTitleFontFile = ''
headerDisplay.DataAxesGrid.ZTitleFontFile = ''
headerDisplay.DataAxesGrid.XLabelFontFile = ''
headerDisplay.DataAxesGrid.YLabelFontFile = ''
headerDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
headerDisplay.PolarAxes.PolarAxisTitleFontFile = ''
headerDisplay.PolarAxes.PolarAxisLabelFontFile = ''
headerDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
headerDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Slice'
slice1 = Slice(Input=header)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [0.6249998807907104, 0.9999998807907104, 0.25]

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=slice1.SliceType)

# Properties modified on slice1.SliceType
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

# Properties modified on slice1.SliceType
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

# show data in view
slice1Display = Show(slice1, renderView1)

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = [None, '']
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 0.1999999761581421
slice1Display.SelectScaleArray = 'None'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'None'
slice1Display.GaussianRadius = 0.009999998807907105
slice1Display.SetScaleArray = [None, '']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = [None, '']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.SelectionCellLabelFontFile = ''
slice1Display.SelectionPointLabelFontFile = ''
slice1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
slice1Display.DataAxesGrid.XTitleFontFile = ''
slice1Display.DataAxesGrid.YTitleFontFile = ''
slice1Display.DataAxesGrid.ZTitleFontFile = ''
slice1Display.DataAxesGrid.XLabelFontFile = ''
slice1Display.DataAxesGrid.YLabelFontFile = ''
slice1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
slice1Display.PolarAxes.PolarAxisTitleFontFile = ''
slice1Display.PolarAxes.PolarAxisLabelFontFile = ''
slice1Display.PolarAxes.LastRadialAxisTextFontFile = ''
slice1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(slice1Display, ('FIELD', 'vtkBlockColors'))

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'vtkBlockColors'
vtkBlockColorsLUT = GetColorTransferFunction('vtkBlockColors')

# get opacity transfer function/opacity map for 'vtkBlockColors'
vtkBlockColorsPWF = GetOpacityTransferFunction('vtkBlockColors')

# Properties modified on vtkBlockColorsLUT
vtkBlockColorsLUT.IndexedColors = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.6299992370489051, 0.6299992370489051, 1.0, 0.6699931334401464, 0.5000076295109483, 0.3300068665598535, 1.0, 0.5000076295109483, 0.7499961852445258, 0.5300068665598535, 0.3499961852445258, 0.7000076295109483, 1.0, 0.7499961852445258, 0.5000076295109483]
vtkBlockColorsLUT.IndexedOpacities = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

# set scalar coloring
ColorBy(slice1Display, ('CELLS', 'proc'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(vtkBlockColorsLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
slice1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'proc'
procLUT = GetColorTransferFunction('proc')

# get opacity transfer function/opacity map for 'proc'
procPWF = GetOpacityTransferFunction('proc')

# create a new 'AMReX Particles Reader'
plt0 = AMReXParticlesReader(FileNames=AllPlotFiles)
plt0.PointArrayStatus = ['id', 'cpu', 'int_comp0', 'real_comp0', 'real_comp1', 'real_comp2']

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# show data in view
plt0Display = Show(plt0, renderView1)

# trace defaults for the display properties.
plt0Display.Representation = 'Surface'
plt0Display.ColorArrayName = [None, '']
plt0Display.OSPRayScaleArray = 'cpu'
plt0Display.OSPRayScaleFunction = 'PiecewiseFunction'
plt0Display.SelectOrientationVectors = 'None'
plt0Display.ScaleFactor = 0.11000000000000001
plt0Display.SelectScaleArray = 'None'
plt0Display.GlyphType = 'Arrow'
plt0Display.GlyphTableIndexArray = 'None'
plt0Display.GaussianRadius = 0.0055000000000000005
plt0Display.SetScaleArray = ['POINTS', 'cpu']
plt0Display.ScaleTransferFunction = 'PiecewiseFunction'
plt0Display.OpacityArray = ['POINTS', 'cpu']
plt0Display.OpacityTransferFunction = 'PiecewiseFunction'
plt0Display.DataAxesGrid = 'GridAxesRepresentation'
plt0Display.SelectionCellLabelFontFile = ''
plt0Display.SelectionPointLabelFontFile = ''
plt0Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
plt0Display.DataAxesGrid.XTitleFontFile = ''
plt0Display.DataAxesGrid.YTitleFontFile = ''
plt0Display.DataAxesGrid.ZTitleFontFile = ''
plt0Display.DataAxesGrid.XLabelFontFile = ''
plt0Display.DataAxesGrid.YLabelFontFile = ''
plt0Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
plt0Display.PolarAxes.PolarAxisTitleFontFile = ''
plt0Display.PolarAxes.PolarAxisLabelFontFile = ''
plt0Display.PolarAxes.LastRadialAxisTextFontFile = ''
plt0Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(plt0Display, ('FIELD', 'vtkBlockColors'))

# show color bar/color legend
plt0Display.SetScalarBarVisibility(renderView1, True)

# create a new 'Glyph'
glyph1 = Glyph(Input=plt0,
    GlyphType='Arrow')
glyph1.OrientationArray = ['POINTS', 'No orientation array']
glyph1.ScaleArray = ['POINTS', 'No scale array']
glyph1.ScaleFactor = 0.11000000000000001
glyph1.GlyphTransform = 'Transform2'

# Properties modified on glyph1
glyph1.GlyphType = 'Sphere'
glyph1.ScaleFactor = 0.02
glyph1.GlyphMode = 'All Points'

# show data in view
glyph1Display = Show(glyph1, renderView1)

# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.ColorArrayName = [None, '']
glyph1Display.OSPRayScaleArray = 'Normals'
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.SelectOrientationVectors = 'None'
glyph1Display.ScaleFactor = 0.11194985434412957
glyph1Display.SelectScaleArray = 'None'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GlyphTableIndexArray = 'None'
glyph1Display.GaussianRadius = 0.005597492717206478
glyph1Display.SetScaleArray = ['POINTS', 'Normals']
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityArray = ['POINTS', 'Normals']
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
glyph1Display.SelectionCellLabelFontFile = ''
glyph1Display.SelectionPointLabelFontFile = ''
glyph1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
glyph1Display.DataAxesGrid.XTitleFontFile = ''
glyph1Display.DataAxesGrid.YTitleFontFile = ''
glyph1Display.DataAxesGrid.ZTitleFontFile = ''
glyph1Display.DataAxesGrid.XLabelFontFile = ''
glyph1Display.DataAxesGrid.YLabelFontFile = ''
glyph1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
glyph1Display.PolarAxes.PolarAxisTitleFontFile = ''
glyph1Display.PolarAxes.PolarAxisLabelFontFile = ''
glyph1Display.PolarAxes.LastRadialAxisTextFontFile = ''
glyph1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(glyph1Display, ('FIELD', 'vtkBlockColors'))

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# reset view to fit data bounds
renderView1.ResetCamera(0.0452507212758, 1.16474926472, 1.88819074631, 1.90768933296, 0.239999994636, 0.259999990463)

# reset view to fit data
renderView1.ResetCamera()

animationScene1.Play()

# current camera placement for renderView1
renderView1.CameraPosition = [0.625, 1.0, 4.907525664877661]
renderView1.CameraFocalPoint = [0.625, 1.0, 0.25]
renderView1.CameraParallelScale = 1.2054563451241194

# save animation
SaveAnimation('pachinko_game.jpeg', renderView1, ImageResolution=[1400, 800],
    FrameWindow=[0, len(AllPlotFiles)-1], 
    # JPEG options
    Quality=85)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [0.625, 1.0, 4.907525664877661]
renderView1.CameraFocalPoint = [0.625, 1.0, 0.25]
renderView1.CameraParallelScale = 1.2054563451241194

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
