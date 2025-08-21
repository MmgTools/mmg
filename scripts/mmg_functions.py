from mmg_classes import *

def setAPIFunctions(api : pythonAPI):

    mmg3d_init_parameters = mmgFunction("MMG3D_Init_parameters","None",
                                        [arg("mesh","MMG5_Mesh",1)])
    api.addFunction(mmg3d_init_parameters)

    mmg3d_set_inputmeshname = mmgFunction("MMG3D_Set_inputMeshName",
                                          "ctypes.c_int",
                                          [arg("mesh","MMG5_Mesh",1),
                                           arg("name","ctypes.c_char_p")])
    api.addFunction(mmg3d_set_inputmeshname)

    mmg3d_set_outputmeshname = mmgFunction("MMG3D_Set_outputMeshName",
                                          "ctypes.c_int",
                                          [arg("mesh","MMG5_Mesh",1),
                                           arg("name","ctypes.c_char_p")])
    api.addFunction(mmg3d_set_outputmeshname)

    mmg3d_set_inputsolname = mmgFunction("MMG3D_Set_inputSolName",
                                          "ctypes.c_int",
                                          [arg("sol","MMG5_Sol",1),
                                           arg("name","ctypes.c_char_p")])
    api.addFunction(mmg3d_set_inputsolname)

    mmg3d_set_outputsolname = mmgFunction("MMG3D_Set_outputSolName",
                                          "ctypes.c_int",
                                          [arg("sol","MMG5_Sol",1),
                                           arg("name","ctypes.c_char_p")])
    api.addFunction(mmg3d_set_outputsolname)


    mmg3d_set_inputparamname = mmgFunction("MMG3D_Set_inputParamName",
                                          "ctypes.c_int",
                                          [arg("mesh","MMG5_Mesh",1),
                                           arg("name","ctypes.c_char_p")])
    api.addFunction(mmg3d_set_inputparamname)

    mmg3d_set_solsize = mmgFunction("MMG3D_Set_solSize",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("sol","MMG5_pSol",1),
                                     arg("typEntity","ctypes.c_int"),
                                     arg("np","MMG5_int"),
                                     arg("typSol","ctypes.c_int")])
    api.addFunction(mmg3d_set_solsize)

    mmg3d_set_solsatverticessize = mmgFunction("MMG3D_Set_solsAtVerticesSize",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("sol","MMG5_pSol",1),
                                     arg("nsols","ctypes.c_int"),
                                     arg("nentities","MMG5_int"),
                                     arg("typSol","ctypes.c_int",1)])
    api.addFunction(mmg3d_set_solsatverticessize)

    mmg3d_set_meshsize = mmgFunction("MMG3D_Set_meshSize",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("np","MMG5_int"),
                                     arg("ne","MMG5_int"),
                                     arg("nprism","MMG5_int"),
                                     arg("nt","MMG5_int"),
                                     arg("nquad","MMG5_int"),
                                     arg("na","MMG5_int")])
    api.addFunction(mmg3d_set_meshsize)

    mmg3d_set_vertex = mmgFunction("MMG3D_Set_vertex",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("c0","ctypes.c_double"),
                                     arg("c1","ctypes.c_double"),
                                     arg("c2","ctypes.c_double"),
                                     arg("ref","MMG5_int"),
                                     arg("pos","MMG5_int")])
    api.addFunction(mmg3d_set_vertex)

    mmg3d_set_vertices = mmgFunction("MMG3D_Set_vertices",
                                     "ctypes.c_int",
                                     [arg("mesh","MMG5_Mesh",1),
                                      arg("vertices","ctypes.c_double",1),
                                      arg("refs","MMG5_int",1)])
    api.addFunction(mmg3d_set_vertices)

    mmg3d_set_tetrahedron = mmgFunction("MMG3D_Set_tetrahedron",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("v0","MMG5_int"),
                                     arg("v1","MMG5_int"),
                                     arg("v2","MMG5_int"),
                                     arg("v3","MMG5_int"),
                                     arg("ref","MMG5_int"),
                                     arg("pos","MMG5_int")])
    api.addFunction(mmg3d_set_tetrahedron)

    mmg3d_set_tetrahedra = mmgFunction("MMG3D_Set_tetrahedra",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("tetra","MMG5_int",1),
                                     arg("refs","MMG5_int",1)])
    api.addFunction(mmg3d_set_tetrahedra)







    mmg3d_load_mesh = mmgFunction("MMG3D_loadMesh","ctypes.c_int",
                                  [arg("mesh","MMG5_Mesh",1),
                                   arg("name","ctypes.c_char_p")])
    api.addFunction(mmg3d_load_mesh)

    mmg3d_libmmg3d = mmgFunction("MMG3D_mmg3dlib","ctypes.c_int",
                                 [arg("mesh","MMG5_Mesh",1),
                                  arg("sol","MMG5_Sol",1)])
    api.addFunction(mmg3d_libmmg3d)
