from writeAPI import *
from mmg_classes import *

point_tmp = mmgClass("Point")
point.addArgs([arg("a","ctypes.c_int"),arg("x","ctypes.c_double"),arg("y","ctypes.c_double")])

printPoint = mmgFunction("printPoint","ctypes.c_int",[arg("p","Point")])

api = pythonAPI()

api.addClass(point_tmp)
api.addFunction(printPoint)

api.addClass(point)
api.addClass(mesh)



api.writeAPI()