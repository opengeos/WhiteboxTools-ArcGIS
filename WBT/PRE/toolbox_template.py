import arcpy


class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Toolbox"
        self.alias = ""

        # List of tool classes associated with this toolbox
        tools = []        
        tools.append(Help)
        tools.append(License)
        tools.append(Version)
        tools.append(ListTools)
        tools.append(ToolHelp)
        tools.append(ToolParameters)
        tools.append(ViewCode)
        tools.append(RunTool)
        tools.append(Update)