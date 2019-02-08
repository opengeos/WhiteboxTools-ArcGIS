import os 
import re
import shutil

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


dir_path = os.path.dirname(os.path.realpath(__file__))
wbt_py = os.path.join(dir_path, "whitebox_tools.py")
wbt_pyt = os.path.join(os.path.dirname(os.path.dirname(dir_path)), "WhiteboxTools.pyt")


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


f_tool_template = open(tool_template_py)
f_toolbox_template = open(toolbox_template_py)
f_tools = open(tools_py, "w")
f_toolbox = open(toolbox_py, "w")

lines_toolbox_template = f_toolbox_template.readlines()
f_toolbox_template.close()

# write toolbox header to final toolbox Python script
for line in lines_toolbox_template:
    f_toolbox.write(line)
f_toolbox.write("\n\n")


lines_tool_template = f_tool_template.readlines()
f_tool_template.close()

category = None  #default tool category


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
                for tool_line in lines_tool_template:
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
                    
                    f_tools.write(cur_line)
                
                f_tools.write("\n\n\n")


f_toolbox.write("\n        self.tools = tools\n\n")

f_tools.close()

f_about = open(about_py)
for line in f_about.readlines():
    f_toolbox.write(line)
f_toolbox.write('\n')
f_about.close()

f_tools = open(tools_py)
for line in f_tools.readlines():
    f_toolbox.write(line)

f_tools.close()
f_toolbox.close()

if os.path.exists(wbt_pyt):
    os.remove(wbt_pyt)
    shutil.copyfile(toolbox_py, wbt_pyt)