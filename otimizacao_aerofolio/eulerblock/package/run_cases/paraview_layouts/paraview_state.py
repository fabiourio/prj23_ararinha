# state file generated using paraview version 5.9.1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# Create a new 'Line Chart View'
lineChartView1 = CreateView('XYChartView')
lineChartView1.ViewSize = [733, 624]
lineChartView1.LegendPosition = [667, 578]
lineChartView1.LeftAxisRangeMinimum = 0.1
lineChartView1.LeftAxisRangeMaximum = 1.9000000000000001
lineChartView1.RightAxisRangeMaximum = 6.66
lineChartView1.TopAxisRangeMaximum = 6.66

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [734, 624]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [0.26144782352209983, 0.006431280975200515, 0.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [0.510598438985746, 0.0021546582800198834, 79.68305380698253]
renderView1.CameraFocalPoint = [0.510598438985746, 0.0021546582800198834, 0.0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.6671531404481784
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.SplitHorizontal(0, 0.500000)
layout1.AssignView(1, renderView1)
layout1.AssignView(2, lineChartView1)
layout1.SetSize(1468, 624)

# ----------------------------------------------------------------
# restore active view
SetActiveView(lineChartView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Legacy VTK Reader'
solutionvtk = LegacyVTKReader(registrationName='solution.vtk', FileNames=['../02_single_run/solution.vtk'])

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(registrationName='CellDatatoPointData1', Input=solutionvtk)
cellDatatoPointData1.CellDataArraytoprocess = ['Mach', 'p', 'resrho', 'resrhoe', 'resrhou', 'resrhov', 'rho', 'rhoe', 'rhou', 'rhov']

# create a new 'Extract Subset'
extractSubset1 = ExtractSubset(registrationName='ExtractSubset1', Input=cellDatatoPointData1)
extractSubset1.VOI = [0, 60, 0, 0, 0, 0]

# create a new 'Plot Data'
plotData1 = PlotData(registrationName='PlotData1', Input=extractSubset1)

# ----------------------------------------------------------------
# setup the visualization in view 'lineChartView1'
# ----------------------------------------------------------------

# show data from plotData1
plotData1Display = Show(plotData1, lineChartView1, 'XYChartRepresentation')

# trace defaults for the display properties.
plotData1Display.CompositeDataSetIndex = [0]
plotData1Display.UseIndexForXAxis = 0
plotData1Display.XArrayName = 'Points_X'
plotData1Display.SeriesVisibility = ['Mach', 'p']
plotData1Display.SeriesLabel = ['Mach', 'Mach', 'p', 'p', 'resrho', 'resrho', 'resrhoe', 'resrhoe', 'resrhou', 'resrhou', 'resrhov', 'resrhov', 'rho', 'rho', 'rhoe', 'rhoe', 'rhou', 'rhou', 'rhov', 'rhov', 'Points_X', 'Points_X', 'Points_Y', 'Points_Y', 'Points_Z', 'Points_Z', 'Points_Magnitude', 'Points_Magnitude']
plotData1Display.SeriesColor = ['Mach', '0', '0', '0', 'p', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845', 'resrho', '0.220004577706569', '0.4899977111467155', '0.7199969481956207', 'resrhoe', '0.30000762951094834', '0.6899977111467155', '0.2899977111467155', 'resrhou', '0.6', '0.3100022888532845', '0.6399938963912413', 'resrhov', '1', '0.5000076295109483', '0', 'rho', '0.6500038147554742', '0.3400015259021897', '0.16000610360875867', 'rhoe', '0', '0', '0', 'rhou', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845', 'rhov', '0.220004577706569', '0.4899977111467155', '0.7199969481956207', 'Points_X', '0.30000762951094834', '0.6899977111467155', '0.2899977111467155', 'Points_Y', '0.6', '0.3100022888532845', '0.6399938963912413', 'Points_Z', '1', '0.5000076295109483', '0', 'Points_Magnitude', '0.6500038147554742', '0.3400015259021897', '0.16000610360875867']
plotData1Display.SeriesPlotCorner = ['Mach', '0', 'Points_Magnitude', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'p', '0', 'resrho', '0', 'resrhoe', '0', 'resrhou', '0', 'resrhov', '0', 'rho', '0', 'rhoe', '0', 'rhou', '0', 'rhov', '0']
plotData1Display.SeriesLabelPrefix = ''
plotData1Display.SeriesLineStyle = ['Mach', '1', 'Points_Magnitude', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'p', '1', 'resrho', '1', 'resrhoe', '1', 'resrhou', '1', 'resrhov', '1', 'rho', '1', 'rhoe', '1', 'rhou', '1', 'rhov', '1']
plotData1Display.SeriesLineThickness = ['Mach', '2', 'Points_Magnitude', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'p', '2', 'resrho', '2', 'resrhoe', '2', 'resrhou', '2', 'resrhov', '2', 'rho', '2', 'rhoe', '2', 'rhou', '2', 'rhov', '2']
plotData1Display.SeriesMarkerStyle = ['Mach', '0', 'Points_Magnitude', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'p', '0', 'resrho', '0', 'resrhoe', '0', 'resrhou', '0', 'resrhov', '0', 'rho', '0', 'rhoe', '0', 'rhou', '0', 'rhov', '0']
plotData1Display.SeriesMarkerSize = ['Mach', '4', 'Points_Magnitude', '4', 'Points_X', '4', 'Points_Y', '4', 'Points_Z', '4', 'p', '4', 'resrho', '4', 'resrhoe', '4', 'resrhou', '4', 'resrhov', '4', 'rho', '4', 'rhoe', '4', 'rhou', '4', 'rhov', '4']

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from cellDatatoPointData1
cellDatatoPointData1Display = Show(cellDatatoPointData1, renderView1, 'StructuredGridRepresentation')

# get color transfer function/color map for 'Mach'
machLUT = GetColorTransferFunction('Mach')
machLUT.RGBPoints = [0.1709036683421313, 0.231373, 0.298039, 0.752941, 0.6598385447823422, 0.865003, 0.865003, 0.865003, 1.148773421222553, 0.705882, 0.0156863, 0.14902]
machLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'Mach'
machPWF = GetOpacityTransferFunction('Mach')
machPWF.Points = [0.1709036683421313, 0.0, 0.5, 0.0, 1.148773421222553, 1.0, 0.5, 0.0]
machPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
cellDatatoPointData1Display.Representation = 'Surface'
cellDatatoPointData1Display.ColorArrayName = ['POINTS', 'Mach']
cellDatatoPointData1Display.LookupTable = machLUT
cellDatatoPointData1Display.OSPRayScaleArray = 'rho'
cellDatatoPointData1Display.OSPRayScaleFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.SelectOrientationVectors = 'None'
cellDatatoPointData1Display.ScaleFactor = 2.9387606836259206
cellDatatoPointData1Display.SelectScaleArray = 'rho'
cellDatatoPointData1Display.GlyphType = 'Arrow'
cellDatatoPointData1Display.GlyphTableIndexArray = 'rho'
cellDatatoPointData1Display.GaussianRadius = 0.146938034181296
cellDatatoPointData1Display.SetScaleArray = ['POINTS', 'rho']
cellDatatoPointData1Display.ScaleTransferFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.OpacityArray = ['POINTS', 'rho']
cellDatatoPointData1Display.OpacityTransferFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.DataAxesGrid = 'GridAxesRepresentation'
cellDatatoPointData1Display.PolarAxes = 'PolarAxesRepresentation'
cellDatatoPointData1Display.ScalarOpacityFunction = machPWF
cellDatatoPointData1Display.ScalarOpacityUnitDistance = 2.8990877897773637

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
cellDatatoPointData1Display.ScaleTransferFunction.Points = [0.7313427679653585, 0.0, 0.5, 0.0, 1.2524234347973364, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
cellDatatoPointData1Display.OpacityTransferFunction.Points = [0.7313427679653585, 0.0, 0.5, 0.0, 1.2524234347973364, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for machLUT in view renderView1
machLUTColorBar = GetScalarBar(machLUT, renderView1)
machLUTColorBar.Title = 'Mach'
machLUTColorBar.ComponentTitle = ''

# set color bar visibility
machLUTColorBar.Visibility = 1

# show color legend
cellDatatoPointData1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# restore active source
SetActiveSource(plotData1)
# ----------------------------------------------------------------


if __name__ == '__main__':
    # generate extracts
    SaveExtracts(ExtractsOutputDirectory='extracts')