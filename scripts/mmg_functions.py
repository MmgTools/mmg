from mmg_classes import *
import os

def setAPIFunctions(api : pythonAPI):

    type_dict = {
        "void"    : "None",
        "int"     : "ctypes.c_int",
        "int8_t"  : "ctypes.c_int8",
        "double"  : "ctypes.c_double",
        "MMG5_int": "MMG5_int",
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
            if (line.find("...") == -1):
                if (line.find(";") == -1):
                    line_buf = line
                    continue

                if (line_buf != ""):
                    line = line_buf.strip(" \n") + line.strip(" ")
                    line_buf = ""

                line_split = line.split(maxsplit=2)
                if (line_split[1] != "extern" ):
                    prototype = line_split[2].split("(",maxsplit=1)
                    name = prototype[0]
                    restype = type_dict[line_split[1]]
                    arglist_str = prototype[1].strip(" ) \n;").split(",")
                    arglist = []
                    for item in arglist_str:
                        item_split = item.split()
                        offset = 0
                        if (len(item_split)>=3):
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
                        else:
                            ptr = 0
                        if (var_type.find("MMG5_") == -1):
                            var_type = type_dict[var_type]
                        arglist.append(arg(var_name,var_type,ptr))

                    func = mmgFunction(name,restype,arglist)
                    api.addFunction(func)
