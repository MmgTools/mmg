from mmg_classes import *

def setAPIFunctions(api : pythonAPI):

    type_dict = {
        "void"    : "None",
        "int"     : "ctypes.c_int",
        "double"  : "ctypes.c_double",
        "MMG5_int": "MMG5_int"
    }

    # read each prototype found in libmmg3d.h
    # manipulate strings to extract each piece of information: return type,
    # function name and list of arguments (type and name of each argument)
    for line in header:
        if (line.find("LIBMMG3D_EXPORT") != -1):
            line_split = line.split(maxsplit=2)
            if (line_split[1] != "extern" ):
                prototype = line_split[2].split("(",maxsplit=1)
                arglist_str = prototype[1].strip(" )\n ;").split(",")
                
                name     = prototype[0]
                restype  = type_dict[line_split[1]]
                arglist = []
                print(arglist_str)
                for item in arglist_str:
                    item_split = item.split()
                    offset = 0
                    if (len(item_split)>=3):
                        offset = 1
                    print(item_split)
                    var_type = item_split[0+offset]
                    var_name = item_split[1+offset]
                    if ((var_type.find("_p") != -1) or (var_type.find("*") != -1)):
                        ptr = 1
                        var_type = var_type.replace("_p","_")
                        var_type = var_type.replace("*"," ")
                    else:
                        ptr = 0
                    arglist.append(arg(var_name,var_type,ptr))
            
            func = mmgFunction(name,restype,argtypes)
            api.addFunction(func)

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

    mmg3d_set_prism = mmgFunction("MMG3D_Set_prism",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("v0","MMG5_int"),
                                     arg("v1","MMG5_int"),
                                     arg("v2","MMG5_int"),
                                     arg("v3","MMG5_int"),
                                     arg("v4","MMG5_int"),
                                     arg("v5","MMG5_int"),
                                     arg("ref","MMG5_int"),
                                     arg("pos","MMG5_int")])
    api.addFunction(mmg3d_set_prism)

    mmg3d_set_prisms = mmgFunction("MMG3D_Set_prisms",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("prisms","MMG5_int",1),
                                     arg("refs","MMG5_int",1)])
    api.addFunction(mmg3d_set_prisms)

    mmg3d_set_triangle = mmgFunction("MMG3D_Set_triangle",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("v0","MMG5_int"),
                                     arg("v1","MMG5_int"),
                                     arg("v2","MMG5_int"),
                                     arg("ref","MMG5_int"),
                                     arg("pos","MMG5_int")])
    api.addFunction(mmg3d_set_triangle)

    mmg3d_set_triangles = mmgFunction("MMG3D_Set_triangles",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("tria","MMG5_int",1),
                                     arg("refs","MMG5_int",1)])
    api.addFunction(mmg3d_set_triangles)

    mmg3d_set_quadrilateral = mmgFunction("MMG3D_Set_quadrilateral",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("v0","MMG5_int"),
                                     arg("v1","MMG5_int"),
                                     arg("v2","MMG5_int"),
                                     arg("v3","MMG5_int"),
                                     arg("ref","MMG5_int"),
                                     arg("pos","MMG5_int")])
    api.addFunction(mmg3d_set_quadrilateral)

    mmg3d_set_quadrilaterals = mmgFunction("MMG3D_Set_quadrilaterals",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("quads","MMG5_int",1),
                                     arg("refs","MMG5_int",1)])
    api.addFunction(mmg3d_set_quadrilaterals)

    mmg3d_set_edge = mmgFunction("MMG3D_Set_edge",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("v0","MMG5_int"),
                                     arg("v1","MMG5_int"),
                                     arg("ref","MMG5_int"),
                                     arg("pos","MMG5_int")])
    api.addFunction(mmg3d_set_edge)

    mmg3d_set_corner = mmgFunction("MMG3D_Set_corner",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("k","MMG5_int")])
    api.addFunction(mmg3d_set_corner)

    mmg3d_unset_corner = mmgFunction("MMG3D_Unset_corner",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("k","MMG5_int")])
    api.addFunction(mmg3d_unset_corner)
 
    mmg3d_set_requiredvertex = mmgFunction("MMG3D_Set_requiredVertex",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("k","MMG5_int")])
    api.addFunction(mmg3d_set_requiredvertex)

    mmg3d_unset_requiredvertex = mmgFunction("MMG3D_Unset_requiredVertex",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("k","MMG5_int")])
    api.addFunction(mmg3d_unset_requiredvertex)

    mmg3d_set_requiredtetrahedron = mmgFunction("MMG3D_Set_requiredTetrahedron",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("k","MMG5_int")])
    api.addFunction(mmg3d_set_requiredtetrahedron)

    mmg3d_unset_requiredtetrahedron = mmgFunction("MMG3D_Unset_requiredTetrahedron",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("k","MMG5_int")])
    api.addFunction(mmg3d_unset_requiredtetrahedron)

    mmg3d_set_requiredtetrahedra = mmgFunction("MMG3D_Set_requiredTetrahedra",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("reqIdx","MMG5_int",1),
                                     arg("nreq","MMG5_int")])
    api.addFunction(mmg3d_set_requiredtetrahedra)

    mmg3d_unset_requiredtetrahedra = mmgFunction("MMG3D_unset_requiredTetrahedra",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("reqIdx","MMG5_int",1),
                                     arg("nreq","MMG5_int")])
    api.addFunction(mmg3d_unset_requiredtetrahedra)

    mmg3d_set_requiredtriangle = mmgFunction("MMG3D_Set_requiredTriangle",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("k","MMG5_int")])
    api.addFunction(mmg3d_set_requiredtriangle)

    mmg3d_unset_requiredtriangle = mmgFunction("MMG3D_Unset_requiredTriangle",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("k","MMG5_int")])
    api.addFunction(mmg3d_unset_requiredtriangle)

    mmg3d_set_requiredtriangles = mmgFunction("MMG3D_Set_requiredTriangles",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("reqIdx","MMG5_int",1),
                                     arg("nreq","MMG5_int")])
    api.addFunction(mmg3d_set_requiredtriangles)

    mmg3d_unset_requiredtriangles = mmgFunction("MMG3D_Unset_requiredTriangles",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("reqIdx","MMG5_int",1),
                                     arg("nreq","MMG5_int")])
    api.addFunction(mmg3d_unset_requiredtriangles)

    mmg3d_set_paralleltriangle = mmgFunction("MMG3D_Set_parallelTriangle",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("k","MMG5_int")])
    api.addFunction(mmg3d_set_paralleltriangle)

    mmg3d_unset_paralleltriangle = mmgFunction("MMG3D_Unset_parallelTriangle",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("k","MMG5_int")])
    api.addFunction(mmg3d_unset_paralleltriangle)

    mmg3d_set_paralleltriangles = mmgFunction("MMG3D_Set_parallelTriangles",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("parIdx","MMG5_int",1),
                                     arg("npar","MMG5_int")])
    api.addFunction(mmg3d_set_paralleltriangles)
 
    mmg3d_unset_paralleltriangles = mmgFunction("MMG3D_Unset_parallelTriangles",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("parIdx","MMG5_int",1),
                                     arg("npar","MMG5_int")])
    api.addFunction(mmg3d_unset_paralleltriangles)

    mmg3d_set_ridge = mmgFunction("MMG3D_Set_ridge",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("k","MMG5_int")])
    api.addFunction(mmg3d_set_ridge)

    mmg3d_unset_ridge = mmgFunction("MMG3D_Unset_ridge",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("k","MMG5_int")])
    api.addFunction(mmg3d_unset_ridge)

    mmg3d_set_requirededge = mmgFunction("MMG3D_Set_requiredEdge",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("k","MMG5_int")])
    api.addFunction(mmg3d_set_requirededge)

    mmg3d_unset_requirededge = mmgFunction("MMG3D_Unset_requiredEdge",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("k","MMG5_int")])
    api.addFunction(mmg3d_unset_requirededge)

    mmg3d_set_normalatvertex = mmgFunction("MMG3D_Set_normalAtVertex",
                                    "ctypes.c_int",
                                    [arg("mesh","MMG5_Mesh",1),
                                     arg("k","MMG5_int"),
                                     arg("n0","ctypes.c_double"),
                                     arg("n1","ctypes.c_double"),
                                     arg("n2","ctypes.c_double")])
    api.addFunction(mmg3d_set_normalatvertex)

    mmg3d_set_scalarsol = mmgFunction("MMG3D_Set_scalarSol",
                                    "ctypes.c_int",
                                    [arg("sol","MMG5_Sol",1),
                                     arg("s","ctypes.c_double"),
                                     arg("pos","MMG5_int")])
    api.addFunction(mmg3d_set_scalarsol)

    mmg3d_load_mesh = mmgFunction("MMG3D_loadMesh","ctypes.c_int",
                                  [arg("mesh","MMG5_Mesh",1),
                                   arg("name","ctypes.c_char_p")])
    api.addFunction(mmg3d_load_mesh)

    mmg3d_libmmg3d = mmgFunction("MMG3D_mmg3dlib","ctypes.c_int",
                                 [arg("mesh","MMG5_Mesh",1),
                                  arg("sol","MMG5_Sol",1)])
    api.addFunction(mmg3d_libmmg3d)
