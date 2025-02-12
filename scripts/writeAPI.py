class mmgFunction():
    def __init__(self,name,rtype,args):
        self.name = name
        self.return_type = rtype
        self.args = args

class arg():
    def __init__(self,arg_name,arg_type):
        self.name = arg_name
        self.type = arg_type

class mmgClass():
    def __init__(self,name,):
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

    def addClass(self,cl):
        self.classes.append(cl)

    def addFunction(self,fn):
        self.functions.append(fn)

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
            f.write("lib." + fn.name + ".argtypes = (")
            for a in fn.args:
                f.write("ctypes.POINTER(" + a.type + "),")
            f.write(")\n")
            f.write("lib." + fn.name + ".restype = ")
            f.write(fn.return_type)
            f.write("\n\n")

        def writeFunction(f, fn):
            indentfn = "    "
            f.write("def " + fn.name + "(")
            for a in fn.args:
                if (not (a == fn.args[-1])):
                    f.write(a.name + ": " + a.type + ",")
                else:
                    f.write(a.name + ": " + a.type)
            f.write("):\n")
            f.write(indentfn + "lib." + fn.name + "(")
            for a in fn.args:
                if (not (a == fn.args[-1])):
                    f.write("ctypes.byref(" + a.name + "),")
                else:
                    f.write("ctypes.byref(" + a.name + ")")
            f.write(")")


        with open(self.ns + ".py", "w") as f:
            f.write(python_header)
            for cl in self.classes:
                writeClass(f,cl)
            for fn in self.functions:
                writeFunctionResArgs(f,fn)
            for fn in self.functions:
                writeFunction(f,fn)

python_header = """import ctypes
import os

libdir  = os.path.dirname(os.path.realpath(__file__))
libpath = os.path.join(libdir,"libtest.dylib")

lib = ctypes.CDLL(libpath)

"""