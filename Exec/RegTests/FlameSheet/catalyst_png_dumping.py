from paraview import catalyst
from paraview.simple import *  # noqa: F403



aMRGaussianPulseSource1 = GetActiveSource()


# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.008, 0.016, 0.1072]
renderView1.CameraFocalPoint = [0.008, 0.016, 0.0]
renderView1.CameraParallelScale = 0.017888543819998316

# get display properties
aMRGaussianPulseSource1Display = GetDisplayProperties(aMRGaussianPulseSource1, view=renderView1)

# set scalar coloring
ColorBy(aMRGaussianPulseSource1Display, ('CELLS', 'x_velocity'))

# rescale color and/or opacity maps used to include current data range
aMRGaussianPulseSource1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
aMRGaussianPulseSource1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'x_velocity'
centroidLUT = GetColorTransferFunction('x_velocity')

# get opacity transfer function/opacity map for 'x_velocity'
centroidPWF = GetOpacityTransferFunction('x_velocity')

# get 2D transfer function for 'x_velocity'
centroidTF2D = GetTransferFunction2D('x_velocity')

# change representation type
aMRGaussianPulseSource1Display.SetRepresentationType('Surface')

# create extractor
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
# trace defaults for the extractor.
pNG1.Trigger = 'Time Step'

# init the 'PNG' selected for 'Writer'
pNG1.Writer.FileName = 'RenderView1_{timestep:06d}{camera}.png'
pNG1.Writer.ImageResolution = [1528, 746]
pNG1.Writer.Format = 'PNG'

options = catalyst.Options()
options.ExtractsOutputDirectory = "datasets"