##################################################################
# Steps for updating WhiteboxTools-ArcGIS
# Step 0: - Update the whitebox_tools.json file in the whiteboxgui repository
# Step 1 - Delete the existing develop branch: git branch -D develop
# Step 2 - Create a new develop branch: git checkout -b develop
# Step 3 - Delete the old WhiteboxTools_win_amd64.zip in the root folder if needed
# Step 4 - Run automation.py
# Step 5 - Commit and push changes
# Step 6 - Merge pull request on GitHub
# Step 7 - Switch to master branch and pull updates: git checkout master | git pull
##################################################################

import json
import os
import re
import shutil
import sys
import whitebox
import urllib.request
from zipfile import ZipFile
from urllib.request import urlopen

wbt = whitebox.WhiteboxTools()


def to_camelcase(name):
    """
    Convert snake_case name to CamelCase name
    """
    return "".join(x.title() for x in name.split("_"))


def to_label(name):
    """
    Convert snake_case name to Title case label
    """
    return " ".join(x.title() for x in name.split("_"))


def to_snakecase(name):
    """
    Convert CamelCase name to snake_case name
    """
    s1 = re.sub("(.)([A-Z][a-z]+)", r"\1_\2", name)
    return re.sub("([a-z0-9])([A-Z])", r"\1_\2", s1).lower()


def write_header(file_path, tool_list):
    """
    Generate Python script header for ArcGIS Python Toolbox
    """
    f_header = open(file_path, "w")
    f_header.write("import arcpy\n")
    f_header.write("import os\n")
    f_header.write("import webbrowser\n")
    f_header.write("from WBT.whitebox_tools import WhiteboxTools\n")
    f_header.write("if sys.version_info < (3, 0):\n")
    f_header.write("    from StringIO import StringIO\n")
    f_header.write("else:\n")
    f_header.write("    from io import StringIO\n\n")
    f_header.write("wbt = WhiteboxTools()\n")
    # ValueList (Dropdown List) for About WhiteboxTools functions (e.g., Run Tool, View Code)
    f_header.write("tool_labels = []\n\n")
    tool_list.sort()
    for tool in tool_list:
        f_header.write('tool_labels.append("{}")\n'.format(tool))
    f_header.write("\n\n")
    f_header.close()


def get_wbt_dict():
    """Generate a dictionary containing information for all tools.

    Returns:
        dict: The dictionary containing information for all tools.
    """

    url = "https://github.com/opengeos/whiteboxgui/raw/master/whiteboxgui/data/whitebox_tools.json"
    response = urlopen(url)
    wbt_dict = json.loads(response.read())

    return wbt_dict


def generate_tool_template(tool):
    """
    Generate function block of each tool for the toolbox
    """
    tool_params = []
    for index, item in enumerate(tool["parameters"]):
        if index < 0:
            tool_params.append(item)
        else:
            item_dup = "{}={}".format(item, item)
            tool_params.append(item_dup)

    lines = []
    lines.append("class {}(object):\n".format(tool["name"]))
    lines.append("    def __init__(self):\n")
    lines.append('        self.label = "{}"\n'.format(tool["label"]))
    lines.append('        self.description = "{}"\n'.format(
        tool["description"]))
    lines.append('        self.category = "{}"\n\n'.format(tool["category"]))
    lines.append("    def getParameterInfo(self):\n")
    # Loop through parameters
    lines.append(define_tool_params(tool["parameters"]))
    lines.append("        return params\n\n")
    lines.append("    def updateParameters(self, parameters):\n")
    lines.append("        return\n\n")
    lines.append("    def updateMessages(self, parameters):\n")
    lines.append("        for param in parameters:\n")
    lines.append("            param_str = param.valueAsText\n")
    lines.append("            if param_str is not None:\n")
    lines.append("                try:\n")
    lines.append("                    desc = arcpy.Describe(param_str)\n")
    lines.append(
        '                    if (".gdb\\\\" in desc.catalogPath) or (".mdb\\\\" in desc.catalogPath):\n'
    )
    lines.append(
        '                        param.setErrorMessage("Datasets stored in a Geodatabase are not supported.")\n'
    )
    lines.append("                except:\n")
    lines.append("                    param.clearMessage()\n")
    lines.append("        return\n\n")
    lines.append("    def execute(self, parameters, messages):\n")
    # Access parameters through parameters[x].valueAsText
    lines.append(define_execute(tool["parameters"]))
    # redirect standard output to tool dialogue
    lines.append("        old_stdout = sys.stdout\n")
    lines.append("        result = StringIO()\n")
    lines.append("        sys.stdout = result\n")

    # line = '        wbt.{}({})\n'.format(to_snakecase(tool['name']), ', '.join(tool['parameters']).replace(", class,", ", cls,"))
    line = "        wbt.{}({})\n".format(
        to_snakecase(tool["name"]),
        ", ".join(tool_params).replace(", class=class,",
                                       ", cls=cls,").replace("input=input", "i=i"),
    )

    # Deal with name conflict with reserved Python functions (and, or, not)
    if tool["name"] == "And":
        line = line.replace("and", "And")
    elif tool["name"] == "Or":
        line = line.replace("or", "Or")
    elif tool["name"] == "Not":
        line = line.replace("not", "Not")
    lines.append(line)
    lines.append("        sys.stdout = old_stdout\n")
    lines.append("        result_string = result.getvalue()\n")
    lines.append("        messages.addMessage(result_string)\n")
    lines.append("        return\n\n\n")
    return lines


def define_tool_params(params):
    """
    Generate function block for each tool parameter
    """
    lines = []
    for param in params:
        items = params[param]
        if items["optional"]:
            parameter_type = "Optional"
        else:
            parameter_type = "Required"

        if "NewFile" in items["parameter_type"]:
            direction = "Output"
        else:
            direction = "Input"

        if param == "class":  # parameter cannot use Python reserved keyword
            param = "cls"

        data_type = get_data_type(items["parameter_type"])

        # if data_type['multi_value'] and param == "i":
        #     param = "inputs"

        if data_type["data_type"] == '"DERasterDataset"' and direction == "Output":
            data_type["data_type"] = '"DEFile"'
            data_type["data_filter"] = '["tif"]'
            # if a filter is used, the parameter must be changed to required.
            parameter_type = "Required"
        elif data_type["data_type"] == '"DERasterDataset"' and direction == "Input":
            data_type["data_type"] = '"GPRasterLayer"'
        elif data_type["data_type"] == '"DEShapefile"' and direction == "Input":
            data_type["data_type"] = '"GPFeatureLayer"'
        elif (
            data_type["data_type"] == '["DERasterDataset", "GPDouble"]'
            and direction == "Input"
        ):
            data_type["data_type"] = '["GPRasterLayer", "GPDouble"]'

        if data_type["data_filter"] == '["html"]':
            parameter_type = "Required"

        lines.append("        {} = arcpy.Parameter(\n".format(param))
        lines.append('            displayName="{}",\n'.format(items["name"]))
        lines.append('            name="{}",\n'.format(param))
        lines.append("            datatype={},\n".format(
            data_type["data_type"]))
        lines.append(
            '            parameterType="{}",\n'.format(parameter_type))
        lines.append('            direction="{}")\n'.format(direction))

        if data_type["multi_value"]:
            lines.append("        {}.multiValue = True\n".format(param))

        if len(data_type["dependency_field"]) > 0:
            if data_type["dependency_field"] == "input":
                data_type["dependency_field"] = "i"
            lines.append(
                "        {}.parameterDependencies = [{}.name]\n".format(
                    param, data_type["dependency_field"]
                )
            )

        if data_type["data_filter"] != "[]":
            if data_type["filter_type"] == '"ValueList"':
                lines.append(
                    '        {}.filter.type = "ValueList"\n'.format(param))

            if (
                data_type["data_filter"]
                in ['["Point"]', '["Polyline"]', '["Polygon"]', '["las", "zip"]']
            ) and direction == "Output":
                pass
            else:
                lines.append(
                    "        {}.filter.list = {}\n".format(
                        param, data_type["data_filter"]
                    )
                )

        if items["default_value"] is not None:
            if (items["default_value"] != "null") and (len(items["default_value"]) > 0):
                if "false" in items["default_value"]:
                    items["default_value"] = False
                elif "true" in items["default_value"]:
                    items["default_value"] = True

                lines.append(
                    "        {}.value = '{}'\n\n".format(
                        param, items["default_value"])
                )
        else:
            lines.append("\n")

    line = "        params = [{}]\n\n".format(", ".join(params))
    if "class" in line:
        line = line.replace(", class,", ", cls,")

    lines.append(line)
    lines = "".join(lines)
    return lines


def define_execute(params):
    """
    Accessing tool parameters
    """
    lines = []
    for index, param in enumerate(params):
        # get the full path to a input raster or vector layer
        param_type = params[param]["parameter_type"]
        inputRasVec = []
        inputRasVec.append({"ExistingFile": "Raster"})
        # inputRasVec.append({'ExistingFileOrFloat': 'Raster'})
        inputRasVec.append({"ExistingFile": {"Vector": "Point"}})
        inputRasVec.append({"ExistingFile": {"Vector": "Line"}})
        inputRasVec.append({"ExistingFile": {"Vector": "Polygon"}})
        inputRasVec.append({"ExistingFile": {"Vector": "LineOrPolygon"}})
        inputRasVec.append({"ExistingFile": {"Vector": "Any"}})

        optional = False

        if "optional" in params[param].keys():
            if params[param]["optional"] == "true":
                optional = True

        # deal with multi-value input
        items = params[param]
        data_type = get_data_type(items["parameter_type"])
        # if data_type['multi_value'] and param == 'i':
        #     param = "inputs"

        if param == "class":
            param = "cls"

        if param == "input":
            param = "i"

        # deal with multi-value inputs
        lines.append(
            "        {} = parameters[{}].valueAsText\n".format(param, index))
        if data_type["multi_value"]:
            lines.append("        if {} is not None:\n".format(param))
            lines.append('            items = {}.split(";")\n'.format(param))
            lines.append("            items_path = []\n")
            lines.append("            for item in items:\n")
            lines.append(
                "                items_path.append(arcpy.Describe(item).catalogPath)\n"
            )
            lines.append(
                '            {} = ";".join(items_path)\n'.format(param))

        if param_type in inputRasVec:
            #     lines.append('        desc = arcpy.Describe({})\n'.format(param))
            #     lines.append('        {} = desc.catalogPath\n'.format(param))
            # if param_type == "Optional":
            lines.append("        if {} is not None:\n".format(param))
            lines.append(
                "            desc = arcpy.Describe({})\n".format(param))
            lines.append("            {} = desc.catalogPath\n".format(param))
        elif param_type == {"ExistingFileOrFloat": "Raster"}:
            lines.append("        if {} is not None:\n".format(param))
            lines.append("            try:\n")
            lines.append(
                "                {} = str(float({}))\n".format(param, param))
            lines.append("            except:\n")
            lines.append(
                "                desc = arcpy.Describe({})\n".format(param))
            lines.append(
                "                {} = desc.catalogPath\n".format(param))

            # lines.append('        if ({} is not None) and {}.isnumeric() == False:\n'.format(param, param))

        # if param == "cell_size":
        #     print(param)

        # if param_type in inputRasVec:
        #     lines.append('        if {} is not None:\n'.format(param))
        #     lines.append('            desc = arcpy.Describe({})\n'.format(param))
        #     lines.append('            {} = desc.catalogPath\n'.format(param))
        #     lines.append('            if (".gdb\\\\" in desc.catalogPath) or (".mdb\\\\" in desc.catalogPath):\n')
        #     lines.append('                 arcpy.AddError("Datasets stored in a Geodatabase are not supported.")\n')
        # elif optional:
        #     lines.append('        if {} is None:\n'.format(param))
        #     lines.append('            {} = None\n'.format(param))

    lines = "".join(lines)
    return lines


def get_data_type(param):
    """
    Convert WhiteboxTools data types to ArcGIS data types
    """
    data_type = '"GPString"'  # default data type
    data_filter = "[]"  # https://goo.gl/EaVNzg
    filter_type = '""'
    multi_value = False
    dependency_field = ""

    # ArcGIS data types: https://goo.gl/95JtFu
    data_types = {
        "Boolean": '"GPBoolean"',
        "Integer": '"GPLong"',
        "Float": '"GPDouble"',
        "String": '"GPString"',
        "StringOrNumber": '["GPString", "GPDouble"]',
        "Directory": '"DEFolder"',
        "Raster": '"DERasterDataset"',
        "Csv": '"DEFile"',
        "Text": '"DEFile"',
        "Html": '"DEFile"',
        "Lidar": '"DEFile"',
        "Vector": '"DEShapefile"',
        "RasterAndVector": '["DERasterDataset", "DEShapefile"]',
        "ExistingFileOrFloat": '["DERasterDataset", "GPDouble"]',
        "ExistingFile": '["DERasterDataset"]',
    }

    vector_filters = {
        "Point": '["Point"]',
        "Points": '["Point"]',
        "Line": '["Polyline"]',
        "Lines": '["Polyline"]',
        "Polygon": '["Polygon"]',
        "Polygons": '["Polygon"]',
        "LineOrPolygon": '["Polyline", "Polygon"]',
        "Any": "[]",
    }

    if type(param) is str:
        data_type = data_types[param]

    else:
        for item in param:
            if item == "FileList":
                multi_value = True
            elif item == "OptionList":
                filter_type = '"ValueList"'
                data_filter = param[item]

            if param[item] == "Csv":
                data_filter = '["csv"]'
            elif param[item] == "Lidar":
                data_filter = '["las", "zip"]'
            elif param[item] == "Html":
                data_filter = '["html"]'

            if type(param[item]) is str:
                data_type = data_types[param[item]]
            elif type(param[item]) is dict:
                sub_item = param[item]
                for sub_sub_item in sub_item:
                    data_type = data_types[sub_sub_item]
                    if data_type == '"DEShapefile"':
                        data_filter = vector_filters[sub_item[sub_sub_item]]
            elif item == "VectorAttributeField":
                data_type = '"Field"'
                dependency_field = param[item][1].replace("--", "")
            else:
                data_type = '"GPString"'

            if param == {"ExistingFileOrFloat": "Raster"}:
                data_type = '["DERasterDataset", "GPDouble"]'

    ret = {}
    ret["data_type"] = data_type
    ret["data_filter"] = data_filter
    ret["filter_type"] = filter_type
    ret["multi_value"] = multi_value
    ret["dependency_field"] = dependency_field

    return ret


def get_github_url(tool_name, category):
    """
    Generate source code link on Github
    """
    # prefix = "https://github.com/jblindsay/whitebox-tools/blob/master/src/tools"
    url = wbt.view_code(tool_name).strip()
    # url = "{}/{}/{}.rs".format(prefix, category, tool_name)
    return url


def get_github_tag(tool_name, category):
    """
    Get GitHub HTML tag
    """
    # prefix = "https://github.com/jblindsay/whitebox-tools/blob/master/src/tools"
    # url = "{}/{}/{}.rs".format(prefix, category, tool_name)
    url = wbt.view_code(tool_name).strip()
    # print(tool_name)
    # if tool_name == "split_vector_lines":
    #     print(url)
    if "RUST_BACKTRACE" in url:
        url = "https://github.com/jblindsay/whitebox-tools/tree/master/whitebox-tools-app/src/tools"
    html_tag = "<a href='{}' target='_blank'>GitHub</a>".format(url)
    return html_tag


def get_book_url(tool_name, category):
    """
    Get link to WhiteboxTools User Manual
    """
    prefix = "https://www.whiteboxgeo.com/manual/wbt_book/available_tools"
    url = "{}/{}.html#{}".format(prefix, category, tool_name)
    return url


def get_book_tag(tool_name, category):
    """
    Get User Manual HTML tag
    """
    prefix = "https://www.whiteboxgeo.com/manual/wbt_book/available_tools"
    url = "{}/{}.html#{}".format(prefix, category, tool_name)
    html_tag = "<a href='{}' target='_blank'>WhiteboxTools User Manual</a>".format(
        url)
    return html_tag


dir_path = os.path.dirname(os.path.realpath(__file__))
wbt_dir = os.path.dirname(dir_path)
root_dir = os.path.dirname(wbt_dir)
wbt_win_zip = os.path.join(root_dir, "WhiteboxTools_win_amd64.zip")

# wbt_py = os.path.join(dir_path, "whitebox_tools.py")
wbt_py = os.path.join(wbt_dir, "whitebox_tools.py")

file_header_py = os.path.join(dir_path, "file_header.py")
file_toolbox_py = os.path.join(dir_path, "file_toolbox.py")
file_tool_py = os.path.join(dir_path, "file_tool.py")
file_about_py = os.path.join(dir_path, "file_about.py")

file_wbt_py = os.path.join(dir_path, "WhiteboxTools.py")
file_wbt_pyt = os.path.join(
    os.path.dirname(os.path.dirname(dir_path)), "WhiteboxTools.pyt"
)

if not os.path.exists(wbt_win_zip):
    print("Downloading WhiteboxTools binary ...")
    url = "https://www.whiteboxgeo.com/WBT_Windows/WhiteboxTools_win_amd64.zip"
    urllib.request.urlretrieve(url, wbt_win_zip)  # Download WhiteboxTools
else:
    print("WhiteboxTools binary already exists.")

print("Decompressing WhiteboxTools_win_amd64.zip ...")
with ZipFile(wbt_win_zip, "r") as zipObj:
    # Extract all the contents of zip file to the root directory
    zipObj.extractall(root_dir)

MACOSX = os.path.join(root_dir, "__MACOSX")
if os.path.exists(MACOSX):
    shutil.rmtree(MACOSX)


tools_dict = get_wbt_dict()
tool_labels = [tools_dict[tool]["label"] for tool in tools_dict.keys()]

write_header(file_header_py, tool_labels)

f_wbt = open(file_wbt_py, "w")

# write toolbox header
with open(file_header_py) as f:
    lines = f.readlines()
    f_wbt.writelines(lines)

# write toolbox class
with open(file_toolbox_py) as f:
    lines = f.readlines()
    f_wbt.writelines(lines)

# write tool class
for tool_name in tools_dict.keys():
    f_wbt.write("        tools.append({})\n".format(tool_name))
f_wbt.write("\n        self.tools = tools\n\n\n")

with open(file_about_py) as f:
    lines = f.readlines()
    f_wbt.writelines(lines)

for tool_name in tools_dict:
    print(tool_name)
    lines = generate_tool_template(tools_dict[tool_name])
    f_wbt.writelines(lines)

f_wbt.close()

# copy WhiteboxTools.py to WhiteboxTool.pyt (final ArcGIS Python Toolbox)
if os.path.exists(file_wbt_pyt):
    os.remove(file_wbt_pyt)
    shutil.copyfile(file_wbt_py, file_wbt_pyt)

outlines = []
with open(wbt_py) as f:
    lines = f.readlines()
    for line in lines:
        if line.strip() == "import urllib.request":
            line = ""
        if line.strip() == "from subprocess import CalledProcessError, Popen, PIPE, STDOUT":
            line = """
from subprocess import CalledProcessError, Popen, PIPE, STDOUT
if sys.version_info.major == 2:
    import urllib2 as urlopen
else:
    import urllib.request as urlopen
            """
        if 'f"={toolname}"' in line:
            line = '                args.append("={}".format(toolname))'
        if line.strip() == 'args2.append(f"--max_procs={val}")':
            line = '            args2.append("--max_procs={}".format(val))'
        if 'f"Warning: Unrecognized extension ext_name {ext_name}' in line:
            line = '                    print("Warning: Unrecognized extension ext_name {}. Installing the GTE instead...".format(ext_name))\n'
        if line.strip() == "for entry in os.scandir(f'./{unzipped_dir_name}'):":
            line = "            for entry in os.scandir('./{}'.format(unzipped_dir_name)):\n"
        if line.strip() == "new_path = entry.path.replace(f'{unzipped_dir_name}', 'plugins')":
            line = "                new_path = entry.path.replace('{}'.format(unzipped_dir_name), 'plugins')\n"
        if line.strip() == "if os.path.exists(f'./{unzipped_dir_name}'):":
            line = "            if os.path.exists('./{}'.format(unzipped_dir_name)):\n"
        if line.strip() == "shutil.rmtree(f'./{unzipped_dir_name}')":
            line = "                shutil.rmtree('./{}'.format(unzipped_dir_name))\n"
        if "urllib.request" in line and "import" not in line:
            line = line.replace("urllib.request", "urlopen")

        outlines.append(line)

with open(wbt_py, "w") as f:
    f.writelines(outlines)
