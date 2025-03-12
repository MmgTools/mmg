from writeAPI import *
from mmg_classes import *

api = pythonAPI()

#printPoint = mmgFunction("printPoint","ctypes.c_int",[arg("p","Point")])
#api.addFunction(printPoint)

for cl in classes:
    api.addClass(cl)

api.writeAPI()