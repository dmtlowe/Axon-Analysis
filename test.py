from panel_splitter import split_panels

image_path = r"Z:\Uni\PhD\Projects\Axon Analysis\Images\24-08-12_OreR_e3DIV_W_2-36_8bit_panel.tif"
output_folder = r"Z:\Uni\PhD\Projects\Axon Analysis\Images\24-08-12_OreR_e3DIV_W_2-36_8bit_panel_split"

panels = split_panels(image_path, output_dir=output_folder)

