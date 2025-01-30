import os 
from paraview import catalyst
from paraview.simple import *  # noqa: F403

# The list of varriables which need to be plotted
variables = ["density","x_velocity","y_velocity"]


# set a new colormap
activeSource1 = GetActiveSource()
current_file_path = os.path.dirname(os.path.abspath(__file__))
json_file_path = os.path.join(current_file_path, "AMR_Colormap.json")
ImportPresets(filename=json_file_path, location=16)


# Create a single global Catalyst options object
options = catalyst.Options()
options.ExtractsOutputDirectory = "In-SITU"

for i,variable in enumerate(variables, start=1):
    if not os.path.exists(f"In-SITU/{variable}"):
        os.makedirs(f"In-SITU/{variable}")
    renderView1 = CreateView("RenderView", registrationName=f"renderView{i}")
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [0.009158251038704926, 0.015785509066906494, 0.1072]   
    renderView1.CameraFocalPoint = [0.009158251038704926, 0.015785509066906494, 0.0]
    renderView1.CameraParallelScale = 0.017888543819998316

    # get display properties
    activeSource1Display = GetDisplayProperties(activeSource1, view=renderView1)

    # change representation type
    activeSource1Display.SetRepresentationType('Surface')

    # set scalar coloring
    ColorBy(activeSource1Display, ('CELLS', variable))

    # show color bar/color legend
    activeSource1Display.SetScalarBarVisibility(renderView1, True)


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
    pNG1.Writer.FileName = os.path.join(variable, '{timestep:06d}{camera}.png')
    pNG1.Writer.ImageResolution = [1528, 746]
    pNG1.Writer.Format = 'PNG'

