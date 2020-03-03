import os
# import arcpy

data_dir = os.path.dirname(os.path.abspath(__file__))
wbt_dir = os.path.dirname(data_dir)
out_dir = os.path.join(os.path.expanduser("~"), "temp")

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

tbx = os.path.join(wbt_dir, "WhiteboxTools.pyt")
arcpy.ImportToolbox(tbx)

input = os.path.join(data_dir, "Wetlands\CLSA_Wetland_Polygons.shp")
output = os.path.join(out_dir, "test.tif")
arcpy.VectorPolygonsToRaster_WBT(input, "FID", output, cell_size=10)
print("Results saved at: {}".format(output))