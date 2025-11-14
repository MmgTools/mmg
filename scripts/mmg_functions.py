from mmg_classes import *
import os, re

def setAPIFunctions(api : pythonAPI):

    type_dict = {
        "void"    : "None",
        "int"     : "ctypes.c_int",
        "int8_t"  : "ctypes.c_int8",
        "double"  : "ctypes.c_double",
        "MMG5_int": "ctypes.c_int",
        "char*"   : "ctypes.c_char_p"
    }

    header_file_path = os.getenv("HEADER_FILE")
    header   = open(header_file_path,"r")
    # read each prototype found in libmmg3d.h
    # manipulate strings to extract each piece of information: return type,
    # function name and list of arguments (type and name of each argument)
    line_buf = ""
    for line in header:
        if (line.find("LIBMMG3D_EXPORT") != -1 or line_buf != ""):
            if (line.find("...") == -1): # ignore variadic functions
                if (line.find(";") == -1): # functions may be declared on several lines in headers
                    line_buf = line
                    continue

                if (line_buf != ""):
                    line = line_buf.strip(" \n") + line.strip(" ")
                    line_buf = ""

                line_split = line.split(maxsplit=2)
                # line_split[0] : "LIBMMG3D_EXPORT"
                # line_split[1] : return type
                # line_split[2] : function name and argument list
                if (line_split[1] != "extern" ):
                    prototype = line_split[2].split("(",maxsplit=1)
                    # prototype[0] : function name
                    # prototype[1] : argument list (types and names)
                    name = prototype[0]
                    restype = type_dict[line_split[1]]
                    arglist_str = prototype[1].strip(" ) \n;").split(",")
                    arglist = []

                    if (name.find("parsar") != -1):
                        continue
                    if (name.find("commonFunc") != -1):
                        continue
                    for item in arglist_str:
                        item_split = item.split()
                        offset = 0
                        if (len(item_split)>=3): # const keyword must be ignored
                            offset = 1
                        var_type = item_split[0+offset]
                        if (var_type != "void"):
                            var_name = item_split[1+offset]
                        if ((var_type.find("_p") != -1) or (var_type.find("*") != -1) and (var_type.find("char") == -1)):
                            ptr = 1
                            var_type = var_type.replace("_p","_")
                            var_type = var_type.replace("*","")
                        elif ((var_name.find("*") != -1 ) and (var_type.find("char") != -1 )):
                            ptr = 0
                            var_name = var_name.replace("*","")
                            var_type = "char*"
                        elif ((var_name.find("*") != -1)):
                            ptr = 1
                            var_name = var_name.replace("*","")
                        else:
                            ptr = 0
                        if (var_type.find("MMG5_") == -1):
                            var_type = type_dict[var_type]
                        if (var_type.find("MMG5_int") != -1):
                            var_type = type_dict[var_type]
                        array_size_brackets = re.search(r"[\d+]", var_name)
                        if (array_size_brackets):
                            array_size = array_size_brackets.group(0).strip("[]")
                            var_name_tmp = var_name.split("[",maxsplit=1)
                            var_name = var_name_tmp[0]
                            var_type = var_type + "*" + array_size
                        if (var_name.find("lambda") != -1):
                            var_name = var_name + "0"
                        arglist.append(arg(var_name,var_type,ptr))

                    func = mmgFunction(name,restype,arglist)
                    api.addFunction(func)
