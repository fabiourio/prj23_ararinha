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
lineChartView1.ViewSize = [458, 410]
lineChartView1.LegendPosition = [353, 215]
lineChartView1.LeftAxisUseCustomRange = 1
lineChartView1.LeftAxisRangeMinimum = -0.2
lineChartView1.LeftAxisRangeMaximum = 2.0
lineChartView1.BottomAxisUseCustomRange = 1
lineChartView1.RightAxisUseCustomRange = 1
lineChartView1.RightAxisRangeMaximum = 6.66
lineChartView1.TopAxisUseCustomRange = 1
lineChartView1.TopAxisRangeMaximum = 6.66

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [458, 410]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [0.2568037407742496, -0.0036104127491496385, 0.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [0.5264421139937535, 0.010579183380217062, 10000.0]
renderView1.CameraFocalPoint = [0.5264421139937535, 0.010579183380217062, 0.0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.5511323419857597
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
layout1.SetSize(917, 410)

# ----------------------------------------------------------------
# restore active view
SetActiveView(lineChartView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Legacy VTK Reader'
import os
dirname = os.path.dirname(__file__)
solutionvtk = LegacyVTKReader(registrationName='solution.vtk', FileNames=[os.path.join(dirname,'..\\02_single_run\\solution.vtk')])

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
plotData1Display.SeriesVisibility = ['Mach', 'p', 'pressure']
plotData1Display.SeriesLabel = ['Mach', 'Mach', 'p', 'p', 'resrho', 'resrho', 'resrhoe', 'resrhoe', 'resrhou', 'resrhou', 'resrhov', 'resrhov', 'rho', 'rho', 'rhoe', 'rhoe', 'rhou', 'rhou', 'rhov', 'rhov', 'Points_X', 'Points_X', 'Points_Y', 'Points_Y', 'Points_Z', 'Points_Z', 'Points_Magnitude', 'Points_Magnitude', 'pressure', 'pressure', 'psi_rho', 'psi_rho', 'psi_rhoe', 'psi_rhoe', 'psi_rhou', 'psi_rhou', 'psi_rhov', 'psi_rhov', 'res_rho', 'res_rho', 'res_rhoe', 'res_rhoe', 'res_rhou', 'res_rhou', 'res_rhov', 'res_rhov']
plotData1Display.SeriesColor = ['Mach', '0', '0', '0', 'p', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845', 'resrho', '0.220004577706569', '0.4899977111467155', '0.7199969481956207', 'resrhoe', '0.30000762951094834', '0.6899977111467155', '0.2899977111467155', 'resrhou', '0.6', '0.3100022888532845', '0.6399938963912413', 'resrhov', '1', '0.5000076295109483', '0', 'rho', '0.6500038147554742', '0.3400015259021897', '0.16000610360875867', 'rhoe', '0', '0', '0', 'rhou', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845', 'rhov', '0.220004577706569', '0.4899977111467155', '0.7199969481956207', 'Points_X', '0.30000762951094834', '0.6899977111467155', '0.2899977111467155', 'Points_Y', '0.6', '0.3100022888532845', '0.6399938963912413', 'Points_Z', '1', '0.5000076295109483', '0', 'Points_Magnitude', '0.6500038147554742', '0.3400015259021897', '0.16000610360875867', 'pressure', '1', '0', '0', 'psi_rho', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845', 'psi_rhoe', '0.220004577706569', '0.4899977111467155', '0.7199969481956207', 'psi_rhou', '0.30000762951094834', '0.6899977111467155', '0.2899977111467155', 'psi_rhov', '0.6', '0.3100022888532845', '0.6399938963912413', 'res_rho', '1', '0.5000076295109483', '0', 'res_rhoe', '0.6500038147554742', '0.3400015259021897', '0.16000610360875867', 'res_rhou', '0', '0', '0', 'res_rhov', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845']
plotData1Display.SeriesPlotCorner = ['Mach', '0', 'Points_Magnitude', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'p', '0', 'pressure', '0', 'psi_rho', '0', 'psi_rhoe', '0', 'psi_rhou', '0', 'psi_rhov', '0', 'res_rho', '0', 'res_rhoe', '0', 'res_rhou', '0', 'res_rhov', '0', 'resrho', '0', 'resrhoe', '0', 'resrhou', '0', 'resrhov', '0', 'rho', '0', 'rhoe', '0', 'rhou', '0', 'rhov', '0']
plotData1Display.SeriesLabelPrefix = ''
plotData1Display.SeriesLineStyle = ['Mach', '1', 'Points_Magnitude', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'p', '1', 'pressure', '1', 'psi_rho', '1', 'psi_rhoe', '1', 'psi_rhou', '1', 'psi_rhov', '1', 'res_rho', '1', 'res_rhoe', '1', 'res_rhou', '1', 'res_rhov', '1', 'resrho', '1', 'resrhoe', '1', 'resrhou', '1', 'resrhov', '1', 'rho', '1', 'rhoe', '1', 'rhou', '1', 'rhov', '1']
plotData1Display.SeriesLineThickness = ['Mach', '2', 'Points_Magnitude', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'p', '2', 'pressure', '2', 'psi_rho', '2', 'psi_rhoe', '2', 'psi_rhou', '2', 'psi_rhov', '2', 'res_rho', '2', 'res_rhoe', '2', 'res_rhou', '2', 'res_rhov', '2', 'resrho', '2', 'resrhoe', '2', 'resrhou', '2', 'resrhov', '2', 'rho', '2', 'rhoe', '2', 'rhou', '2', 'rhov', '2']
plotData1Display.SeriesMarkerStyle = ['Mach', '0', 'Points_Magnitude', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'p', '0', 'pressure', '0', 'psi_rho', '0', 'psi_rhoe', '0', 'psi_rhou', '0', 'psi_rhov', '0', 'res_rho', '0', 'res_rhoe', '0', 'res_rhou', '0', 'res_rhov', '0', 'resrho', '0', 'resrhoe', '0', 'resrhou', '0', 'resrhov', '0', 'rho', '0', 'rhoe', '0', 'rhou', '0', 'rhov', '0']
plotData1Display.SeriesMarkerSize = ['Mach', '4', 'Points_Magnitude', '4', 'Points_X', '4', 'Points_Y', '4', 'Points_Z', '4', 'p', '4', 'pressure', '4', 'psi_rho', '4', 'psi_rhoe', '4', 'psi_rhou', '4', 'psi_rhov', '4', 'res_rho', '4', 'res_rhoe', '4', 'res_rhou', '4', 'res_rhov', '4', 'resrho', '4', 'resrhoe', '4', 'resrhou', '4', 'resrhov', '4', 'rho', '4', 'rhoe', '4', 'rhou', '4', 'rhov', '4']

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from cellDatatoPointData1
cellDatatoPointData1Display = Show(cellDatatoPointData1, renderView1, 'StructuredGridRepresentation')

# get color transfer function/color map for 'Mach'
machLUT = GetColorTransferFunction('Mach')
machLUT.RGBPoints = [0.27665520344357664, 0.231373, 0.298039, 0.752941, 0.7162229884367939, 0.865003, 0.865003, 0.865003, 1.155790773430011, 0.705882, 0.0156863, 0.14902]
machLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'Mach'
machPWF = GetOpacityTransferFunction('Mach')
machPWF.Points = [0.27665520344357664, 0.0, 0.5, 0.0, 1.155790773430011, 1.0, 0.5, 0.0]
machPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
cellDatatoPointData1Display.Representation = 'Surface'
cellDatatoPointData1Display.ColorArrayName = ['POINTS', 'Mach']
cellDatatoPointData1Display.LookupTable = machLUT
cellDatatoPointData1Display.SelectTCoordArray = 'None'
cellDatatoPointData1Display.SelectNormalArray = 'None'
cellDatatoPointData1Display.SelectTangentArray = 'None'
cellDatatoPointData1Display.OSPRayScaleArray = 'rho'
cellDatatoPointData1Display.OSPRayScaleFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.SelectOrientationVectors = 'None'
cellDatatoPointData1Display.ScaleFactor = 2.9377462401438104
cellDatatoPointData1Display.SelectScaleArray = 'rho'
cellDatatoPointData1Display.GlyphType = 'Arrow'
cellDatatoPointData1Display.GlyphTableIndexArray = 'rho'
cellDatatoPointData1Display.GaussianRadius = 0.14688731200719052
cellDatatoPointData1Display.SetScaleArray = ['POINTS', 'rho']
cellDatatoPointData1Display.ScaleTransferFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.OpacityArray = ['POINTS', 'rho']
cellDatatoPointData1Display.OpacityTransferFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.DataAxesGrid = 'GridAxesRepresentation'
cellDatatoPointData1Display.PolarAxes = 'PolarAxesRepresentation'
cellDatatoPointData1Display.ScalarOpacityFunction = machPWF
cellDatatoPointData1Display.ScalarOpacityUnitDistance = 2.8978579955061714

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
cellDatatoPointData1Display.ScaleTransferFunction.Points = [0.6792717756956008, 0.0, 0.5, 0.0, 1.232404246110169, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
cellDatatoPointData1Display.OpacityTransferFunction.Points = [0.6792717756956008, 0.0, 0.5, 0.0, 1.232404246110169, 1.0, 0.5, 0.0]

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