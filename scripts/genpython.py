from writeAPI import *
from mmg_classes import *
from mmg_functions import *

api = pythonAPI()
setAPIClasses(api)
setAPIFunctions(api)
api.write_api()

#printPoint = mmgFunction("printPoint","ctypes.c_int",[arg("p","Point")])
#api.addFunction(printPoint)
