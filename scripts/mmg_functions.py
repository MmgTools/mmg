from mmg_classes import *

def setAPIFunctions(api : pythonAPI):

    mmg3d_init_mesh = mmgFunction("MMG3D_Init_mesh","ctypes.c_int",[arg("*starter","ctypes.c_int")])
    api.addFunction(mmg3d_init_mesh)
    mmg3d_init_filenames = mmgFunction("MMG3D_Init_fileNames","None",[arg("mesh","MMG5_Mesh"),arg("sol","MMG5_Sol")])
    api.addFunction(mmg3d_init_filenames)
    mmg3d_init_parameters = mmgFunction("MMG3D_Init_parameters","None",[arg("mesh","MMG5_Mesh")])

    mmg3d_load_mesh = mmgFunction("MMG3D_loadMesh","ctypes.c_int",[arg("mesh","MMG5_Mesh"),arg("name","ctypes.c_char_p")])
    api.addFunction(mmg3d_load_mesh)

    mmg3d_libmmg3d = mmgFunction("MMG3D_mmg3dlib","ctypes.c_int",[arg("mesh","MMG5_Mesh"),arg("sol","MMG5_Sol")])
    api.addFunction(mmg3d_libmmg3d)

    mmg3d_free_mesh = mmgFunction("MMG3D_Free_all","ctypes.c_int",[arg("*starter","ctypes.c_int")])
    api.addFunction(mmg3d_free_mesh)
