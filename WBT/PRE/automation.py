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
    f_header = open(file_path, "w")
    f_header.write("import arcpy\n")
    f_header.write("from WBT.whitebox_tools import WhiteboxTools\n")
    f_header.write("wbt = WhiteboxTools()\n")
    f_header.write("wbt.set_verbose_mode(True)\n\n")
    f_header.write("tool_labels = []\n")
    tool_list.sort()
    for tool in tool_list:
        f_header.write('tool_labels.append("{}")\n'.format(tool))
    f_header.write("\n\n")
    f_header.close()


def extract_func_params(line):
    if line.startswith("def"):
        line = line.replace("self, i,", "input,")
        line = line.replace("self, i=None,", "input,")
        line = line.replace("self, output, i=None,", "input, output,")
        line = line.replace("self, ", "")
        line = line.replace("callback=None", "")
        line = line.replace("def ", "")
        line = line.replace(":", "")
        line = line.replace(", )", ")")
        line = line[line.index("(")+1 : len(line) - 1]
        params_temp = line.split(",")
        params = []
        for param in params_temp:
            params.append(param.strip())

    return params


def get_tool_params(tool_name):

    out_str = wbt.tool_parameters(tool_name)
    start_index = out_str.index('[') + 1
    end_index = len(out_str.strip()) - 2
    params = out_str[start_index : end_index]
    # print(params)

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
        # print("{}\n".format(item))
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


def get_param_types(tools):

    parameter_types = []
    for tool in tools:
        params = tools[tool]['parameters']
        for param in params:
            param_type = params[param]['parameter_type']
            if param_type not in parameter_types:
                parameter_types.append(param_type)
    return parameter_types


def generate_tool_template(tool):
    lines = []
    lines.append('class {}(object):\n'.format(tool['name']))
    lines.append('    def __init__(self):\n')
    lines.append('        self.label = "{}"\n'.format(tool['label']))
    lines.append('        self.description = "{}"\n'.format(tool['description']))
    lines.append('        self.category = "{}"\n\n'.format(tool['category']))
    lines.append('    def getParameterInfo(self):\n')
    lines.append(define_tool_params(tool['parameters']))
    lines.append('        return params\n\n')
    lines.append('    def updateParameters(self, parameters):\n')
    lines.append('        return\n\n')
    lines.append('    def updateMessages(self, parameters):\n')
    lines.append('        return\n\n')
    lines.append('    def execute(self, parameters, messages):\n')
    lines.append(define_execute(tool['parameters']))
    line = '        messages.addMessage(wbt.{}({}))\n'.format(to_snakecase(tool['name']), ', '.join(tool['parameters']).replace(", class,", ", class1,"))
    if tool['name'] == "And":
        line = line.replace("and", "And")
    elif tool['name'] == "Or":
        line = line.replace("or", "Or")
    elif tool['name'] == "Not":
        line = line.replace("not", "Not")
    lines.append(line)
    lines.append('        return\n\n\n')
    return lines


def define_tool_params(params):
    lines = []
    # num_params = len(params)
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

        # data_type = '"GPString"'

        if param == "class":   # parameter cannot use Python reserved keyword
            param = "class1"

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
        # lines.append('            datatype="{}",\n'.format(data_type))
        lines.append('            parameterType="{}",\n'.format(parameter_type))
        lines.append('            direction="{}")\n'.format(direction))

        if len(data_type['dependency_field']) > 0:
            lines.append('        {}.parameterDependencies = [{}.name]\n'.format(param, data_type['dependency_field']))

        if data_type['data_filter'] != '[]':
            if data_type['filter_type'] == '"ValueList"':
                lines.append('        {}.filter.type = "ValueList"\n'.format(param))
            # lines.append('        {}.filter.type = {}\n'.format(param, data_type['filter_type']))
            lines.append('        {}.filter.list = {}\n'.format(param, data_type['data_filter']))

        if (items['default_value'] != 'null') and (len(items['default_value']) > 0):
            lines.append('\n        {}.value = {}\n\n'.format(param, items['default_value']))
        else:
            lines.append('\n')
        
    line = '        params = [{}]\n\n'.format(', '.join(params))    
    if "class" in line:
        line = line.replace(", class,", ", class1,")
    
    lines.append(line)
    lines = ''.join(lines)
    return lines


def define_execute(params):
    lines = []
    for index, param in enumerate(params):
        if param == 'class':
            param = "class1"
        lines.append('        {} = parameters[{}].valueAsText\n'.format(param, index))
    lines = ''.join(lines)    
    return lines


def get_data_type(param):
    data_type = '"GPString"'  # default data type
    data_filter = '[]'   # https://goo.gl/EaVNzg
    filter_type = '""'
    multi_value = False
    dependency_field = ''

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

            # elif item == 'NewFile' and param[item] == 'Raster':
            #     data_filter = '["tif"]'
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


dir_path = os.path.dirname(os.path.realpath(__file__))
wbt_py = os.path.join(dir_path, "whitebox_tools.py")

file_header_py = os.path.join(dir_path, "file_header.py")
file_toolbox_py = os.path.join(dir_path, "file_toolbox.py")
file_tool_py = os.path.join(dir_path, "file_tool.py")
file_about_py = os.path.join(dir_path, "file_about.py")

file_wbt_py = os.path.join(dir_path, "WhiteboxTools.py")
file_wbt_pyt = os.path.join(os.path.dirname(os.path.dirname(dir_path)), "WhiteboxTools.pyt")

tool_template_py = os.path.join(dir_path, "tool_template.py")           # code chuck for each tool
toolbox_template_py = os.path.join(dir_path, "toolbox_template.py")     # code chuck for toolbox header
about_py = os.path.join(dir_path, "about.py")  

tools_py = os.path.join(dir_path, "tools.py")
toolbox_py = os.path.join(dir_path, "toolbox.py")

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
                # print(category)
            
            if line.startswith("def"):
                func_title = line.replace("def", "", 1).strip().split("(")[0]
                func_name = to_camelcase(func_title)

                func_label = to_label(func_title)
                tool_labels.append(func_label)
                func_desc = lines[index+1].replace('"""', '').strip()

                func_dict = {}
                func_dict['name'] = func_name
                func_dict["category"] = category
                func_dict["label"] = func_label
                func_dict["description"] = func_desc
                
                # print("{}: {}".format(tool_index, func_name))
                tool_index = tool_index + 1
                func_params = get_tool_params(func_name)
                func_dict["parameters"] = func_params

                # func_params = extract_func_params(line)
                # print(func_params)


                # print("{}: {} - {} - {}".format(category, func_name, func_label, description))
                tools_dict[func_name] = func_dict


write_header(file_header_py, tool_labels)

f_wbt = open(file_wbt_py, "w")

with open(file_header_py) as f:
    lines = f.readlines()
    f_wbt.writelines(lines)

with open(file_toolbox_py) as f:
    lines = f.readlines()
    f_wbt.writelines(lines)

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

if os.path.exists(file_wbt_pyt):
    os.remove(file_wbt_pyt)
    shutil.copyfile(file_wbt_py, file_wbt_pyt)




# types = []
# param_types = get_param_types(tools_dict)
# for param in param_types:
#     param_type = type(param)
#     if param_type == str:
#         print("{} - String".format(param))
#         # if type(param_type) not in types:
#         #     types.append(param)
#     else:
#         print("{} - Dictionary".format(param))
#         # if type(param_type) not in types:
#         #     types.append(param.keys())
#         for item in param:
#             if (param[item] not in types) and type(param[item]) != str:
#                 types.append(param[item])


# print(len(param_types))
# for item in types:
#     print(item)

# for item in types:
#     print(get_data_type(item))

# tool_name = "LidarElevationSlice"
# params = tools_dict[tool_name]['parameters']
# print(params)
# print(len(params))

# print(wbt.tool_parameters(tool_name))
# lines = define_tool_params(tools_dict['AddPointCoordinatesToTable'])
# for line in lines:
#     print(line)

# lines = define_tool_params(params)
# for line in lines:
#     print(line)