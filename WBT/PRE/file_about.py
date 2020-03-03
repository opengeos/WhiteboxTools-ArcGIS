class Help(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Help"
        self.description = "Help description for WhiteboxTools"
        self.category = "About WhiteboxTools"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
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
        messages.addMessage(wbt.help())
        return


class License(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "License"
        self.description = "License information for WhiteboxTools."
        self.category = "About WhiteboxTools"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
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
        messages.addMessage(wbt.license())
        return


class Version(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Version"
        self.description = "Version information for WhiteboxTools."
        self.category = "About WhiteboxTools"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
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
        messages.addMessage(wbt.version())
        return


class ListTools(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "List Tools"
        self.description = "All available tools in WhiteboxTools."
        self.category = "About WhiteboxTools"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""

        # First parameter
        param0 = arcpy.Parameter(
            displayName="Keywords",
            name="keywords",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        # param0.multiValue = True
        param0.value = "lidar"
        params = [param0]
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
        param0 = parameters[0].valueAsText

        if param0 is None:
            tools = wbt.list_tools()
        else:
            tools = wbt.list_tools([param0])
            
        for index, tool in enumerate(sorted(tools)):
            messages.addMessage("{}. {}: {}".format(index + 1, tool, tools[tool]))
        return


class ToolHelp(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Tool Help"
        self.description = "Help description for a specific tool."
        self.category = "About WhiteboxTools"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        tool_name = arcpy.Parameter(
            displayName="Select a tool",
            name="tool_name",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        tool_name.value = "Lidar Info"
        tool_name.filter.type = "ValueList"
        tool_name.filter.list = tool_labels

        params = [tool_name]
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
        param0 = parameters[0].valueAsText
        tool_name = param0.replace(" ", "").strip()
        messages.addMessage(wbt.tool_help(tool_name))
        return


class ToolParameters(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Tool Parameters"
        self.description = "Tool parameter descriptions for a specific tool."
        self.category = "About WhiteboxTools"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        tool_name = arcpy.Parameter(
            displayName="Select a tool",
            name="tool_name",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        tool_name.value = "Lidar Info"
        tool_name.filter.type = "ValueList"
        tool_name.filter.list = tool_labels

        params = [tool_name]
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
        param0 = parameters[0].valueAsText
        tool_name = param0.replace(" ", "").strip()
        messages.addMessage(wbt.tool_parameters(tool_name))
        return


class ViewCode(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "View Code"
        self.description = "Source code for a specific tool."
        self.category = "About WhiteboxTools"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        tool_name = arcpy.Parameter(
            displayName="Select a tool",
            name="tool_name",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        tool_name.value = "Lidar Info"
        tool_name.filter.type = "ValueList"
        tool_name.filter.list = tool_labels

        params = [tool_name]
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
        param0 = parameters[0].valueAsText
        tool_name = param0.replace(" ", "").strip()        
        messages.addMessage("Opening default browser...")
        webbrowser.get('windows-default').open(wbt.view_code(tool_name))
        messages.addMessage(wbt.view_code(tool_name))
        return


class RunTool(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Run Tool"
        self.description = "Runs a tool and specifies tool arguments."
        self.category = "About WhiteboxTools"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        tool_name = arcpy.Parameter(
            displayName="Select a tool",
            name="tool_name",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        tool_name.value = "Breach Depressions"
        tool_name.filter.type = "ValueList"
        tool_name.filter.list = tool_labels

        args = arcpy.Parameter(
            displayName="Arguments",
            name="agrs",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        args.value = '--dem="/path/to/DEM.tif"  --output="/path/to/output.tif"'

        params = [tool_name, args]
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
        param0 = parameters[0].valueAsText
        args = parameters[1].valueAsText
        tool_name = param0.replace(" ", "").strip()
        dir_path = os.path.dirname(os.path.realpath(__file__))    
        exe_path = os.path.join(dir_path, "WBT/whitebox_tools.exe")
        cmd = '{} --run={} {}'.format(exe_path, tool_name, args)
        if "-v" not in cmd:
            cmd = cmd + ' -v'  
        messages.addMessage(cmd)  
        messages.addMessage(os.popen(cmd).read().rstrip())
        return


