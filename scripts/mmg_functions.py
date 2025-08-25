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

LIBMMG3D_EXPORT int  MMG3D_Set_scalarSol(MMG5_pSol met, double s,MMG5_int pos);
    mmg3d_set_scalarsol = mmgFunction("MMG3D_Set_scalarSol",
                                    "ctypes.c_int",
                                    [arg("sol","MMG5_Sol",1),
                                     arg("s","ctypes.c_double"),
                                     arg("pos","MMG5_int")])
    api.addFunction(mmg3d_set_scalarsol)
  LIBMMG3D_EXPORT int  MMG3D_Set_scalarSols(MMG5_pSol met, double *s);
  LIBMMG3D_EXPORT int MMG3D_Set_vectorSol(MMG5_pSol met, double vx,double vy, double vz,
 LIBMMG3D_EXPORT int MMG3D_Set_vectorSols(MMG5_pSol met, double *sols);
 LIBMMG3D_EXPORT int MMG3D_Set_tensorSol(MMG5_pSol met, double m11,double m12, double m13,
  LIBMMG3D_EXPORT int MMG3D_Set_tensorSols(MMG5_pSol met, double *sols);
  LIBMMG3D_EXPORT int  MMG3D_Set_ithSol_inSolsAtVertices(MMG5_pSol sol,int i, double* s,MMG5_int pos);
  LIBMMG3D_EXPORT int  MMG3D_Set_ithSols_inSolsAtVertices(MMG5_pSol sol,int i, double* s);
 LIBMMG3D_EXPORT void MMG3D_Set_handGivenMesh(MMG5_pMesh mesh);
 LIBMMG3D_EXPORT int MMG3D_Chk_meshData(MMG5_pMesh mesh, MMG5_pSol met);
 LIBMMG3D_EXPORT int  MMG3D_Set_iparameter(MMG5_pMesh mesh,MMG5_pSol sol, int iparam,
                                           MMG5_int val);
 LIBMMG3D_EXPORT int  MMG3D_Set_dparameter(MMG5_pMesh mesh,MMG5_pSol sol, int dparam,
                                           double val);
  LIBMMG3D_EXPORT int  MMG3D_Set_localParameter(MMG5_pMesh mesh, MMG5_pSol sol, int typ,
                                                MMG5_int ref,double hmin,double hmax,double hausd);
 LIBMMG3D_EXPORT int  MMG3D_Set_multiMat(MMG5_pMesh mesh, MMG5_pSol sol,MMG5_int ref,int split,
                                         MMG5_int rmin, MMG5_int rplus);
LIBMMG3D_EXPORT int  MMG3D_Set_lsBaseReference(MMG5_pMesh mesh, MMG5_pSol sol,MMG5_int br);
  LIBMMG3D_EXPORT int  MMG3D_Get_meshSize(MMG5_pMesh mesh, MMG5_int* np, MMG5_int* ne,MMG5_int *nprism, MMG5_int* nt,
                                          MMG5_int* nquad, MMG5_int* na);
  LIBMMG3D_EXPORT int  MMG3D_Get_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int* typEntity,
                                         MMG5_int* np,int* typSol);
  LIBMMG3D_EXPORT int  MMG3D_Get_solsAtVerticesSize(MMG5_pMesh mesh, MMG5_pSol* sol,int *nsols,
                                                    MMG5_int* nentities,int* typSol);
  LIBMMG3D_EXPORT int  MMG3D_Get_vertex(MMG5_pMesh mesh, double* c0, double* c1, double* c2,
                                        MMG5_int* ref,int* isCorner, int* isRequired);
 LIBMMG3D_EXPORT int  MMG3D_GetByIdx_vertex(MMG5_pMesh mesh, double* c0, double* c1, double* c2, MMG5_int* ref,
                                            int* isCorner, int* isRequired,MMG5_int idx);
  LIBMMG3D_EXPORT int  MMG3D_Get_vertices(MMG5_pMesh mesh, double* vertices, MMG5_int* refs,
                                          int* areCorners, int* areRequired);
  LIBMMG3D_EXPORT int  MMG3D_Get_tetrahedron(MMG5_pMesh mesh, MMG5_int* v0, MMG5_int* v1, MMG5_int* v2,
                                             MMG5_int* v3,MMG5_int* ref, int* isRequired);
  LIBMMG3D_EXPORT int  MMG3D_Get_tetrahedra(MMG5_pMesh mesh, MMG5_int* tetra,MMG5_int* refs,
                                            int* areRequired);
  LIBMMG3D_EXPORT int  MMG3D_Get_prism(MMG5_pMesh mesh, MMG5_int* v0, MMG5_int* v1, MMG5_int* v2,
                                       MMG5_int* v3,MMG5_int* v4,MMG5_int* v5,MMG5_int* ref, int* isRequired);
  LIBMMG3D_EXPORT int  MMG3D_Get_prisms(MMG5_pMesh mesh, MMG5_int* prisms,MMG5_int* refs,
                                        int* areRequired);
  LIBMMG3D_EXPORT int  MMG3D_Get_triangle(MMG5_pMesh mesh, MMG5_int* v0, MMG5_int* v1, MMG5_int* v2, MMG5_int* ref,
                                          int* isRequired);
  LIBMMG3D_EXPORT int  MMG3D_Get_triangles(MMG5_pMesh mesh, MMG5_int* tria, MMG5_int* refs,
                                           int* areRequired);
  LIBMMG3D_EXPORT int  MMG3D_Get_quadrilateral(MMG5_pMesh mesh, MMG5_int* v0, MMG5_int* v1, MMG5_int* v2,MMG5_int* v3,
                                               MMG5_int* ref, int* isRequired);
  LIBMMG3D_EXPORT int  MMG3D_Get_quadrilaterals(MMG5_pMesh mesh, MMG5_int* quads, MMG5_int* refs,
                                                int* areRequired);
  LIBMMG3D_EXPORT int  MMG3D_Get_edge(MMG5_pMesh mesh, MMG5_int* e0, MMG5_int* e1, MMG5_int* ref,
                                      int* isRidge, int* isRequired);
  LIBMMG3D_EXPORT int MMG3D_Set_edges(MMG5_pMesh mesh, MMG5_int *edges, MMG5_int* refs);
  LIBMMG3D_EXPORT int MMG3D_Get_edges(MMG5_pMesh mesh,MMG5_int *edges,MMG5_int* refs,
                                      int *areRidges,int *areRequired);
 LIBMMG3D_EXPORT int  MMG3D_Get_normalAtVertex(MMG5_pMesh mesh, MMG5_int k, double *n0, double *n1,
                                               double *n2) ;
  LIBMMG3D_EXPORT double MMG3D_Get_tetrahedronQuality(MMG5_pMesh mesh, MMG5_pSol met, MMG5_int k);
 LIBMMG3D_EXPORT int  MMG3D_Get_scalarSol(MMG5_pSol met, double* s);
 LIBMMG3D_EXPORT int  MMG3D_Get_scalarSols(MMG5_pSol met, double* s);
 LIBMMG3D_EXPORT int MMG3D_Get_vectorSol(MMG5_pSol met, double* vx, double* vy, double* vz);
 LIBMMG3D_EXPORT int MMG3D_Get_vectorSols(MMG5_pSol met, double* sols);
 LIBMMG3D_EXPORT int MMG3D_Get_tensorSol(MMG5_pSol met, double *m11,double *m12, double *m13,
                                         double *m22,double *m23, double *m33);
 LIBMMG3D_EXPORT int MMG3D_Get_tensorSols(MMG5_pSol met, double *sols);
  LIBMMG3D_EXPORT int  MMG3D_Get_ithSol_inSolsAtVertices(MMG5_pSol sol,int i, double* s,MMG5_int pos);
  LIBMMG3D_EXPORT int  MMG3D_Get_ithSols_inSolsAtVertices(MMG5_pSol sol,int i, double* s);
 LIBMMG3D_EXPORT int MMG3D_Get_iparameter(MMG5_pMesh mesh, MMG5_int iparam);
  LIBMMG3D_EXPORT int  MMG3D_Add_tetrahedron(MMG5_pMesh mesh, MMG5_int v0, MMG5_int v1,
                                             MMG5_int v2, MMG5_int v3, MMG5_int ref);
  LIBMMG3D_EXPORT MMG5_int  MMG3D_Add_vertex(MMG5_pMesh mesh, double c0, double c1,
                                             double c2, MMG5_int ref);



















    mmg3d_load_mesh = mmgFunction("MMG3D_loadMesh","ctypes.c_int",
                                  [arg("mesh","MMG5_Mesh",1),
                                   arg("name","ctypes.c_char_p")])
    api.addFunction(mmg3d_load_mesh)

    mmg3d_libmmg3d = mmgFunction("MMG3D_mmg3dlib","ctypes.c_int",
                                 [arg("mesh","MMG5_Mesh",1),
                                  arg("sol","MMG5_Sol",1)])
    api.addFunction(mmg3d_libmmg3d)
