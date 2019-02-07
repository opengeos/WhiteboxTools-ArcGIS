import arcpy
import os
from os import path

exe_name = "whitebox_tools.exe"
exe_path = path.join(path.dirname(path.abspath(__file__)), "WBT", exe_name)

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Toolbox"
        self.alias = ""

        # List of tool classes associated with this toolbox
        self.tools = [AbsoluteValue]


class AbsoluteValue(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Absolute Value"
        self.description = "Calculates the absolute value of every cell in a raster."
        self.category = "Math and Stats Analysis"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = []
        # First parameter
        param0 = arcpy.Parameter(
            displayName="Input raster file",
            name="input",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")
        params.append(param0)

        # Second parameter
        param1 = arcpy.Parameter(
            displayName="Output raster file",
            name="output",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output")
        params.append(param1)

        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        if not path.exists(exe_path):
            messages.addErrorMessage("WhiteboxTools executable could not be found!")
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText

        cmd = exe_path + " --run=AbsoluteValue" + " --input=" + input + " --output=" + output
        os.popen(cmd)

        return
