from mmg_classes import *

def setAPIFunctions(api : pythonAPI):

    mmg3d_init_mesh = mmgFunction("MMG3D_Init_mesh","ctypes.c_int",[arg("*starter","ctypes.c_int")])
    api.addFunction(mmg3d_init_mesh)
    mmg3D_init_filenames = mmgFunction("MMG3D_Init_fileNames","None",[arg("mesh","MMG5_Mesh"),arg("sol","MMG5_Sol")])
    api.addFunction(mmg3D_init_filenames)