import os 
import re

dir_path = os.path.dirname(os.path.realpath(__file__))
wbt_py = os.path.join(dir_path, "whitebox_tools.py")
tools_py = os.path.join(dir_path, "tools.py")
tool_template_py = os.path.join(dir_path, "tool_template.py")
toolbox_template_py = os.path.join(dir_path, "toolbox_template.py")

def to_camelcase(name):
    '''
    Convert snake_case name to CamelCase name 
    '''
    return ''.join(x.title() for x in name.split('_'))

def to_label(name):
    '''
    Convert snake_case name to CamelCase name 
    '''
    return ' '.join(x.title() for x in name.split('_'))

def to_snakecase(name):
    '''
    Convert CamelCase name to snake_case name 
    '''
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()


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

category = None


ff = open(tools_py, "w")
f_temp = open(tool_template_py)
tool_template_lines = f_temp.readlines()
f_temp.close()

f_toolbox = open(toolbox_template_py, "a")

with open(wbt_py) as f:
    lines = f.readlines()

    for index, line in enumerate(lines):
        if index > 360:            
            line = line.strip()

            # Create an R script for each toolbox
            if line in toolboxes:
                category = toolboxes[line]
                # print(category)
            
            if line.startswith("def"):
                title = line.replace("def", "").strip().split("(")[0]
                label = to_label(title)
                name = to_camelcase(title)
                description = lines[index+1].replace('"""', '').strip()
                print("{}: {} - {} - {}".format(category, name, label, description))

                f_toolbox.write("        tools.append({})\n".format(name))

                cur_line = ""
                for tool_line in tool_template_lines:
                    if tool_line.strip() == "class Tool(object):":
                        cur_line = tool_line.replace("Tool", name)
                    elif tool_line.strip() == 'self.label = "Tool"':
                        cur_line = tool_line.replace("Tool", label)
                    elif tool_line.strip() == 'self.description = ""':
                        cur_line = tool_line.replace('self.description = ""', 'self.description = "{}"'.format(description))
                    elif tool_line.strip() == 'self.canRunInBackground = False':
                        cur_line = tool_line.replace("canRunInBackground", "category")
                        cur_line = cur_line.replace("False", '"{}"'.format(category))
                    else:
                        cur_line = tool_line
                    
                    ff.write(cur_line)
                
                ff.write("\n\n")


f_toolbox.write("        self.tools = tools")

ff.close()
f_toolbox.close()
