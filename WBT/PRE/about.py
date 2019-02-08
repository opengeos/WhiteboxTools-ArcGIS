class Help(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Help"
        self.description = "Help description for WhiteboxTools"
        self.category = "About WhiteboxTools"

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
        return


class License(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "License"
        self.description = "License information for WhiteboxTools."
        self.category = "About WhiteboxTools"

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
        return


class Version(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Version"
        self.description = "Version information for WhiteboxTools."
        self.category = "About WhiteboxTools"

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
        return



class ListTools(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "List Tools"
        self.description = "All available tools in WhiteboxTools."
        self.category = "About WhiteboxTools"

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
        return

## name conflict with ArcGIS Python Toolbox
# class Toolbox(object):
#     def __init__(self):
#         """Define the tool (tool name is the name of the class)."""
#         self.label = "Toolbox"
#         self.description = "The toolbox for a specific tool."
#         self.category = "About WhiteboxTools"

#     def getParameterInfo(self):
#         """Define parameter definitions"""
#         params = None
#         return params

#     def isLicensed(self):
#         """Set whether tool is licensed to execute."""
#         return True

#     def updateParameters(self, parameters):
#         """Modify the values and properties of parameters before internal
#         validation is performed.  This method is called whenever a parameter
#         has been changed."""
#         return

#     def updateMessages(self, parameters):
#         """Modify the messages created by internal validation for each tool
#         parameter.  This method is called after internal validation."""
#         return

#     def execute(self, parameters, messages):
#         """The source code of the tool."""
#         return


class ToolHelp(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Tool Help"
        self.description = "Help description for a specific tool."
        self.category = "About WhiteboxTools"

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
        return


class ToolParameters(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Tool Parameters"
        self.description = "Tool parameter descriptions for a specific tool."
        self.category = "About WhiteboxTools"

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
        return



class ViewCode(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "View Code"
        self.description = "Source code for a specific tool."
        self.category = "About WhiteboxTools"

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
        return


class RunTool(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Run Tool"
        self.description = "Runs a tool and specifies tool arguments."
        self.category = "About WhiteboxTools"

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
        return


class Update(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Update"
        self.description = "Download the latest version of WhiteboxTools"
        self.category = "About WhiteboxTools"

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
        return