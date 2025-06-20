import os 
from paraview import catalyst
from paraview.simple import *  # noqa: F403


variable = "density"

# set a new colormap
activeSource1 = GetActiveSource()
current_file_path = os.path.dirname(os.path.abspath(__file__))
json_file_path = os.path.join(current_file_path, "AMR_Colormap.json")
ImportPresets(filename=json_file_path, location=16)

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

renderView1.InteractionMode = '3D'
# current camera placement for renderView1
renderView1.CameraPosition = [-0.04904332546407211, 0.012375295615402978, -0.041652301422069396]
renderView1.CameraFocalPoint = [0.00800000037997961, 0.01600000075995922, 0.00800000037997961]
renderView1.CameraViewUp = [-0.030664588682070313, 0.9988190134439574, -0.037686355406156155]
renderView1.CameraParallelScale = 0.019595918873021582


# get display properties
# activeSource1Display = GetDisplayProperties(activeSource1, view=renderView1)

# set scalar coloring
# ColorBy(activeSource1Display, ('CELLS', variable))

# show color bar/color legend
#activeSource1Display.SetScalarBarVisibility(renderView1, True)

# change representation type
#activeSource1Display.SetRepresentationType('Surface')

# create a new 'Clip'
clip1 = Clip(registrationName='Clip1', Input=activeSource1)

# Properties modified on clip1.ClipType
clip1.ClipType.Origin = [0.005609997317906828, 0.020109299961136952 , 0.006674611818954681]


# Properties modified on clip1.ClipType
clip1.ClipType.Normal = [-0.5246100083723001, 0.819989536944208, -0.2289137357557183]

# show data in view
clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip1Display.Representation = 'Surface'

ColorBy(clip1Display, ('CELLS', variable))



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

SetActiveSource(pNG1)

options = catalyst.Options()
options.ExtractsOutputDirectory = "In-SITU/1Node1Task"