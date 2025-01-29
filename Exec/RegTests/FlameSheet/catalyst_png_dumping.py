import os 
from paraview import catalyst
from paraview.simple import *  # noqa: F403


variable = "Y(CH3)"

# set a new colormap
activeSource1 = GetActiveSource()
current_file_path = os.path.dirname(os.path.abspath(__file__))
json_file_path = os.path.join(current_file_path, "AMR_Colormap.json")
ImportPresets(filename=json_file_path, location=16)

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.009308476712947374, 0.015579097547887256, 0.1072]
renderView1.CameraFocalPoint = [0.009308476712947374, 0.015579097547887256, 0.0]
renderView1.CameraParallelScale = 0.021645138022197957

# get display properties
activeSource1Display = GetDisplayProperties(activeSource1, view=renderView1)

# set scalar coloring
ColorBy(activeSource1Display, ('CELLS', variable))

# show color bar/color legend
activeSource1Display.SetScalarBarVisibility(renderView1, True)

# change representation type
activeSource1Display.SetRepresentationType('Surface')

# get color transfer function/color map for 'source'
sourceLUT = GetColorTransferFunction(variable)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
sourceLUT.ApplyPreset('AMR_Colormap', True)

# get color legend/bar for sourceLUT in view renderView1
sourceLUTColorBar = GetScalarBar(sourceLUT, renderView1)

# change scalar bar placement
sourceLUTColorBar.WindowLocation = 'Any Location'
sourceLUTColorBar.Position = [0.6713147410358565, 0.18417462482946795]
sourceLUTColorBar.ScalarBarLength = 0.3299999999999999

# Properties modified on sourceLUTColorBar
sourceLUTColorBar.LabelFontSize = 10
sourceLUTColorBar.ScalarBarLength = 0.5
sourceLUTColorBar.TitleFontSize = 10

# create extractor
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
# trace defaults for the extractor.
pNG1.Trigger = 'Time Step'

# init the 'PNG' selected for 'Writer'
pNG1.Writer.FileName = '{timestep:06d}{camera}.png'
pNG1.Writer.ImageResolution = [1528, 746]
pNG1.Writer.Format = 'PNG'

options = catalyst.Options()
options.ExtractsOutputDirectory = "CH3"