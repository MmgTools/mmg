from mmg_classes import *

def setAPIFunctions(api : pythonAPI):

    mmg3d_init_parameters = mmgFunction("MMG3D_Init_parameters","None",
                                        [arg("mesh","MMG5_Mesh",1)])
    api.addFunction(mmg3d_init_parameters)

    mmg3d_load_mesh = mmgFunction("MMG3D_loadMesh","ctypes.c_int",
                                  [arg("mesh","MMG5_Mesh",1),
                                   arg("name","ctypes.c_char_p")])
    api.addFunction(mmg3d_load_mesh)

    mmg3d_libmmg3d = mmgFunction("MMG3D_mmg3dlib","ctypes.c_int",
                                 [arg("mesh","MMG5_Mesh",1),
                                  arg("sol","MMG5_Sol",1)])
    api.addFunction(mmg3d_libmmg3d)
