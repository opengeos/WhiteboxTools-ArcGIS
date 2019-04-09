import os 
import re
import shutil
import sys
import ast
import whitebox

wbt = whitebox.WhiteboxTools()

def to_camelcase(name):
    '''
    Convert snake_case name to CamelCase name 
    '''
    return ''.join(x.title() for x in name.split('_'))

def to_label(name):
    '''
    Convert snake_case name to Title case label 
    '''
    return ' '.join(x.title() for x in name.split('_'))

def to_snakecase(name):
    '''
    Convert CamelCase name to snake_case name 
    '''
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()


def write_header(file_path, tool_list):
    '''
    Generate Python script header for ArcGIS Python Toolbox  
    '''
    f_header = open(file_path, "w")
    f_header.write("import arcpy\n")
    f_header.write("import os\n")
    f_header.write("from WBT.whitebox_tools import WhiteboxTools\n")
    f_header.write('if sys.version_info < (3, 0):\n')
    f_header.write('    from StringIO import StringIO\n')
    f_header.write('else:\n')
    f_header.write('    from io import StringIO\n\n')
    f_header.write("wbt = WhiteboxTools()\n")
    # ValueList (Dropdown List) for About WhiteboxTools functions (e.g., Run Tool, View Code)
    f_header.write("tool_labels = []\n\n")      
    tool_list.sort()
    for tool in tool_list:
        f_header.write('tool_labels.append("{}")\n'.format(tool))
    f_header.write("\n\n")
    f_header.close()


def get_tool_params(tool_name):
    '''
    Convert tool parameters output string to a dictionary  
    '''
    out_str = wbt.tool_parameters(tool_name)
    start_index = out_str.index('[') + 1
    end_index = len(out_str.strip()) - 2
    params = out_str[start_index : end_index]

    sub_params = params.split('{"name"')
    param_list = []

    for param in sub_params:
        param = param.strip()
        if len(param) > 0:
            item = '"name"' + param
            item = item[ : item.rfind("}")].strip()
            param_list.append(item)

    params_dict = {}
    for item in param_list:
        param_dict = {}
        index_name = item.find("name")
        index_flags = item.find("flags")
        index_description = item.find("description")
        index_parameter_type = item.find("parameter_type")
        index_default_value = item.find("default_value")
        index_optional = item.find("optional")

        name = item[index_name - 1 : index_flags - 2].replace('"name":', '')
        name = name.replace('"', '')
        param_dict['name'] = name

        flags = item[index_flags - 1 : index_description -2].replace('"flags":', '')
        if flags.count("--") == 1 :
            flags = flags.split('--')[1][: -2]
        elif flags.count("--") == 2:
            flags = flags.split('--')[2][: -2]
        else:
            flags = flags.split('-')[1][: -2]
        param_dict['flags'] = flags

        desc = item[index_description - 1 : index_parameter_type - 2].replace('"description":', '')
        desc = desc.replace('"', '')
        param_dict['description'] = desc

        param_type = item[index_parameter_type - 1 : index_default_value - 2].replace('"parameter_type":', '')
        param_type = ast.literal_eval(param_type)
        param_dict['parameter_type'] = param_type

        default_value = item[index_default_value - 1 : index_optional - 2].replace('"default_value":', '')
        param_dict['default_value'] = default_value

        optional = item[index_optional - 1 :].replace('"optional":', '')
        param_dict['optional'] = optional

        params_dict[flags] = param_dict

    return params_dict


# def get_param_types(tools):
#     '''
#     Get unique parameter types  
#     '''
#     parameter_types = []
#     for tool in tools:
#         params = tools[tool]['parameters']
#         for param in params:
#             param_type = params[param]['parameter_type']
#             if param_type not in parameter_types:
#                 parameter_types.append(param_type)
#     return parameter_types


def generate_tool_template(tool):
    '''
    Generate function block of each tool for the toolbox  
    '''
    tool_params = []
    for item in tool['parameters']:
        item_dup = "{}={}".format(item, item)
        tool_params.append(item_dup)

    lines = []
    lines.append('class {}(object):\n'.format(tool['name']))
    lines.append('    def __init__(self):\n')
    lines.append('        self.label = "{}"\n'.format(tool['label']))
    lines.append('        self.description = "{}"\n'.format(tool['description']))
    lines.append('        self.category = "{}"\n\n'.format(tool['category']))
    lines.append('    def getParameterInfo(self):\n')
    # Loop through parameters 
    lines.append(define_tool_params(tool['parameters']))
    lines.append('        return params\n\n')
    lines.append('    def updateParameters(self, parameters):\n')
    lines.append('        return\n\n')
    lines.append('    def updateMessages(self, parameters):\n')
    lines.append('        return\n\n')
    lines.append('    def execute(self, parameters, messages):\n')
    # Access parameters throught parameters[x].valueAsText
    lines.append(define_execute(tool['parameters']))
    # redirect standard output to tool dialogue
    lines.append('        old_stdout = sys.stdout\n')
    lines.append('        result = StringIO()\n')
    lines.append('        sys.stdout = result\n')
    # line = '        wbt.{}({})\n'.format(to_snakecase(tool['name']), ', '.join(tool['parameters']).replace(", class,", ", cls,"))
    line = '        wbt.{}({})\n'.format(to_snakecase(tool['name']), ', '.join(tool_params).replace(", class=class,", ", cls=cls,"))

    # Deal with name conflict with reserved Python functions (and, or, not)
    if tool['name'] == "And":
        line = line.replace("and", "And")
    elif tool['name'] == "Or":
        line = line.replace("or", "Or")
    elif tool['name'] == "Not":
        line = line.replace("not", "Not")
    lines.append(line)
    lines.append('        sys.stdout = old_stdout\n')
    lines.append('        result_string = result.getvalue()\n')
    lines.append('        messages.addMessage(result_string)\n')
    lines.append('        return\n\n\n')
    return lines


def define_tool_params(params):
    '''
    Generate function block for each tool parameter 
    '''
    lines = []
    for param in params:
        items = params[param]
        if items['optional'] == 'false':
            parameter_type="Required"
        else:
            parameter_type="Optional"
        
        if 'NewFile' in items['parameter_type']:
            direction="Output"
        else:
            direction="Input"

        if param == "class":   # parameter cannot use Python reserved keyword
            param = "cls"

        data_type = get_data_type(items['parameter_type'])

        if data_type['data_type'] == '"DERasterDataset"' and direction == "Output":
            data_type['data_type'] = '"DEFile"'
            data_type['data_filter'] = '["tif"]'
            parameter_type = "Required"   # if a filter is used, the parameter must be changed to required.

        if data_type['data_filter'] == '["html"]':
            parameter_type = "Required"

        lines.append('        {} = arcpy.Parameter(\n'.format(param))
        lines.append('            displayName="{}",\n'.format(items['name']))
        lines.append('            name="{}",\n'.format(param))
        lines.append('            datatype={},\n'.format(data_type['data_type']))
        lines.append('            parameterType="{}",\n'.format(parameter_type))
        lines.append('            direction="{}")\n'.format(direction))

        if data_type['multi_value']:
            lines.append('        {}.multiValue = True\n'.format(param))

        if len(data_type['dependency_field']) > 0:
            lines.append('        {}.parameterDependencies = [{}.name]\n'.format(param, data_type['dependency_field']))

        if data_type['data_filter'] != '[]':
            if data_type['filter_type'] == '"ValueList"':
                lines.append('        {}.filter.type = "ValueList"\n'.format(param))
            lines.append('        {}.filter.list = {}\n'.format(param, data_type['data_filter']))

        if (items['default_value'] != 'null') and (len(items['default_value']) > 0):
            lines.append('\n        {}.value = {}\n\n'.format(param, items['default_value']))
        else:
            lines.append('\n')
        
    line = '        params = [{}]\n\n'.format(', '.join(params))    
    if "class" in line:
        line = line.replace(", class,", ", cls,")
    
    lines.append(line)
    lines = ''.join(lines)
    return lines


def define_execute(params):
    '''
    Accessing tool parameters
    '''
    lines = []
    for index, param in enumerate(params):
        if param == 'class':
            param = "cls"
        lines.append('        {} = parameters[{}].valueAsText\n'.format(param, index))
    lines = ''.join(lines)    
    return lines


def get_data_type(param):
    '''
    Convert WhiteboxTools data types to ArcGIS data types
    '''
    data_type = '"GPString"'  # default data type
    data_filter = '[]'   # https://goo.gl/EaVNzg
    filter_type = '""'
    multi_value = False
    dependency_field = ''

    # ArcGIS data types: https://goo.gl/95JtFu
    data_types = {
        'Boolean': '"GPBoolean"',
        'Integer': '"GPLong"',
        'Float': '"GPDouble"',
        'String': '"GPString"',
        'StringOrNumber': '["GPString", "GPDouble"]',
        'Directory': '"DEFolder"',
        'Raster': '"DERasterDataset"',
        'Csv': '"DEFile"',
        'Text': '"DEFile"',
        'Html': '"DEFile"',
        'Lidar': '"DEFile"',
        'Vector': '"DEShapefile"',
        'RasterAndVector': '["DERasterDataset", "DEShapefile"]'
    }

    vector_filters = {
        'Point': '["Point"]',
        'Line': '["Polyline"]',
        'Polygon': '["Polygon"]',
        'LineOrPolygon': '["Polyline", "Polygon"]',
        'Any': '[]'
    }

    if type(param) is str:
        data_type = data_types[param]

    else:
        for item in param:
            if item == 'FileList':
                multi_value = True
            elif item == 'OptionList':
                filter_type = '"ValueList"'
                data_filter = param[item]

            if param[item] == 'Csv':
                data_filter = '["csv"]'
            elif param[item] == 'Lidar':
                data_filter = '["las", "zip"]'
            elif param[item] == 'Html':
                data_filter = '["html"]'

            if type(param[item]) is str:
                data_type = data_types[param[item]]
            elif type(param[item]) is dict:
                sub_item = param[item]
                for sub_sub_item in sub_item:
                    data_type = data_types[sub_sub_item]
                    if data_type == '"DEShapefile"':
                        data_filter = vector_filters[sub_item[sub_sub_item]] 
            elif item == 'VectorAttributeField':
                data_type = '"Field"'
                dependency_field = param[item][1].replace('--', '')                
            else:
                data_type = '"GPString"'
    ret = {}
    ret['data_type'] = data_type
    ret['data_filter'] = data_filter
    ret['filter_type'] = filter_type
    ret['multi_value'] = multi_value
    ret['dependency_field'] = dependency_field

    return ret


def get_github_url(tool_name, category):
    '''
    Generate source code link on Github 
    '''    
    # prefix = "https://github.com/jblindsay/whitebox-tools/blob/master/src/tools"
    url = wbt.view_code(tool_name).strip()
    # url = "{}/{}/{}.rs".format(prefix, category, tool_name)
    return url


def get_github_tag(tool_name, category):
    '''
    Get GitHub HTML tag
    '''    
    # prefix = "https://github.com/jblindsay/whitebox-tools/blob/master/src/tools"
    # url = "{}/{}/{}.rs".format(prefix, category, tool_name)
    url = wbt.view_code(tool_name).strip()
    html_tag = "<a href='{}' target='_blank'>GitHub</a>".format(url)
    return html_tag


def get_book_url(tool_name, category):
    '''
    Get link to WhiteboxTools User Mannual 
    '''    
    prefix = "https://jblindsay.github.io/wbt_book/available_tools"
    url = "{}/{}.html#{}".format(prefix, category, tool_name)
    return url
    


def get_book_tag(tool_name, category):
    '''
    Get User Manual HTML tag
    '''    
    prefix = "https://jblindsay.github.io/wbt_book/available_tools"
    url = "{}/{}.html#{}".format(prefix, category, tool_name)
    html_tag = "<a href='{}' target='_blank'>WhiteboxTools User Manual</a>".format(url)
    return html_tag


dir_path = os.path.dirname(os.path.realpath(__file__))
wbt_py = os.path.join(dir_path, "whitebox_tools.py")

file_header_py = os.path.join(dir_path, "file_header.py")
file_toolbox_py = os.path.join(dir_path, "file_toolbox.py")
file_tool_py = os.path.join(dir_path, "file_tool.py")
file_about_py = os.path.join(dir_path, "file_about.py")

file_wbt_py = os.path.join(dir_path, "WhiteboxTools.py")
file_wbt_pyt = os.path.join(os.path.dirname(os.path.dirname(dir_path)), "WhiteboxTools.pyt")

toolboxes = {
    "# Data Tools #": "Data Tools",
    "# GIS Analysis #": "GIS Analysis",
    "# Geomorphometric Analysis #": "Geomorphometric Analysis",
    "# Hydrological Analysis #": "Hydrological Analysis",
    "# Image Processing Tools #": "Image Processing Tools",
    "# LiDAR Tools #": "LiDAR Tools",
    "# Math and Stats Tools #": "Math and Stats Tools",
    "# Stream Network Analysis #": "Stream Network Analysis"
}

github_cls = {
    "Data Tools": "data_tools",
    "GIS Analysis": "gis_analysis",
    "Geomorphometric Analysis": "terrain_analysis",
    "Hydrological Analysis": "hydro_analysis",
    "Image Processing Tools": "image_analysis",
    "LiDAR Tools": "lidar_analysis",
    "Math and Stats Tools": "math_stat_analysis",
    "Stream Network Analysis": "stream_network_analysis"
}

book_cls = {
    "Data Tools": "data_tools",
    "GIS Analysis": "gis_analysis",
    "Geomorphometric Analysis": "geomorphometric_analysis",
    "Hydrological Analysis": "hydrological_analysis",
    "Image Processing Tools": "image_processing_tools",
    "LiDAR Tools": "lidar_tools",
    "Math and Stats Tools": "mathand_stats_tools",
    "Stream Network Analysis": "stream_network_analysis"
}





tools_dict = {}
tool_labels = []
category = ''

tool_index = 1

with open(wbt_py) as f:
    lines = f.readlines()

    for index, line in enumerate(lines):
        if index > 360:            
            line = line.strip()

            if line in toolboxes:
                category = toolboxes[line]
            
            if line.startswith("def"):
                func_title = line.replace("def", "", 1).strip().split("(")[0]
                func_name = to_camelcase(func_title)

                func_label = to_label(func_title)
                tool_labels.append(func_label)
                func_desc = lines[index+1].replace('"""', '').strip()

                github_tag = get_github_tag(func_title, github_cls[category])
                book_tag = get_book_tag(func_name, book_cls[category])
                full_desc = "{} View detailed help documentation on {} and source code on {}.".format(func_desc, book_tag, github_tag)

                func_dict = {}
                func_dict['name'] = func_name
                func_dict["category"] = category
                func_dict["label"] = func_label
                func_dict["description"] = full_desc
                
                tool_index = tool_index + 1
                func_params = get_tool_params(func_name)
                func_dict["parameters"] = func_params
                tools_dict[func_name] = func_dict

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
    lines = generate_tool_template(tools_dict[tool_name])
    f_wbt.writelines(lines)

f_wbt.close()

# copy WhiteboxTools.py to WhiteboxTool.pyt (final ArcGIS Python Toolbox)
if os.path.exists(file_wbt_pyt):
    os.remove(file_wbt_pyt)
    shutil.copyfile(file_wbt_py, file_wbt_pyt)