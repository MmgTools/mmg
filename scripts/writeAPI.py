import os

class mmgFunction():
    def __init__(self,name,rtype,args,encode):
        self.name = name
        self.return_type = rtype
        self.args = args
        self.str_encode = encode

class arg():
    def __init__(self,arg_name,arg_type,arg_pointer=0):
        self.name    = arg_name
        self.type    = arg_type
        self.pointer = arg_pointer

class mmgClass():
    def __init__(self,name):
        self.name = name
        self.args = []

    def addArgs(self,args):
        for f in args:
            self.args.append(f)

class pythonAPI:

    def __init__(self, namespace="mmg"):
        self.ns = namespace
        self.classes   = []
        self.functions = []
        self.typenames = {}

    def addClass(self,cl):
        self.classes.append(cl)

    def addFunction(self,fn):
        self.functions.append(fn)

    def addType(self,key,value):
        self.typenames[key] = value

    def writeAPI(self):
        def writeClass(f, cl):
            f.write("class " + cl.name + "(ctypes.Structure):\n")
            indentcl = "    "
            indentfl = 4*indentcl
            f.write(indentcl + "_fields_ = [")
            for a in cl.args:
                if  (a == cl.args[0]):
                    f.write("(\""+a.name+"\","+a.type+"),\n")
                elif(a == cl.args[-1]):
                    f.write(indentfl+"(\""+a.name+"\","+a.type+")")
                else:
                    f.write(indentfl+"(\""+a.name+"\","+a.type+"),\n")
            f.write("]\n\n")

        def writeFunctionResArgs(f,fn):
            if (not fn.str_encode):
                f.write("lib." + fn.name + ".argtypes = (")
                for a in fn.args:
                    if (a.pointer):
                        f.write("ctypes.POINTER(" + a.type + "),")
                    else:
                        f.write(a.type + ",")
                f.write(")\n")
            f.write("lib." + fn.name + ".restype  = ")
            f.write(fn.return_type)
            f.write("\n\n")

        def writeFunction(f, fn):
            indentfn = "    "
            f.write("def " + fn.name + "(")
            for a in fn.args:
                if (a == "str"):
                    if (not (a == fn.args[-1])):
                        f.write(a.name + ",")
                    else:
                        f.write(a.name)
                else:
                    if (not (a == fn.args[-1])):
                        f.write(a.name + ": " + a.type + ",")
                    else:
                        f.write(a.name + ": " + a.type)
            f.write("):\n")
            f.write(indentfn)
            if (fn.str_encode):
                f.write("if (isinstance(name,str)):\n")
                f.write(indentfn + indentfn + "name = ctypes.c_char_p(name.encode('utf-8'))\n")
                f.write(indentfn)
            if (not (fn.return_type == "None")):
                f.write("ier = ")
            f.write("lib." + fn.name + "(")
            for a in fn.args:
                if (not (a == fn.args[-1])):
                    if (a.pointer):
                        f.write("ctypes.byref(" + a.name + "),")
                    else:
                        f.write(a.name + ",")
                else:
                    if (a.pointer):
                        f.write("ctypes.byref(" + a.name + ")")
                    else:
                        f.write(a.name)
            f.write(")")
            f.write("\n")
            if (not (fn.return_type == "None")):
                f.write(indentfn + "return ier\n")
            f.write("\n")


        with open(self.ns + ".py", "w") as f:
            f.write(python_header)
            for cl in self.classes:
                writeClass(f,cl)
            for fn in self.functions:
                writeFunctionResArgs(f,fn)
            for fn in self.functions:
                writeFunction(f,fn)

python_header = """\
import ctypes
import os

lib = ctypes.CDLL("{libpath}")

MMG5_int = "ctypes.c_int"

""".format(libpath=os.getenv("SHARED_LIB_FILE"))
