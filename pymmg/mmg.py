import ctypes
import os
from enum import IntEnum

lib2d = ctypes.CDLL(os.getenv("SHARED_LIB_FILE2D"))
lib3d = ctypes.CDLL(os.getenv("SHARED_LIB_FILE3D"))
libs  = ctypes.CDLL(os.getenv("SHARED_LIB_FILES"))

MMG5_int = "ctypes.c_int"

MMG2D_LMAX = 1024

class MMG3D_Param(IntEnum):
    MMG3D_IPARAM_verbose = 0
    MMG3D_IPARAM_mem = 1
    MMG3D_IPARAM_debug = 2
    MMG3D_IPARAM_angle = 3
    MMG3D_IPARAM_iso = 4
    MMG3D_IPARAM_isosurf = 5
    MMG3D_IPARAM_nofem = 6
    MMG3D_IPARAM_opnbdy = 7
    MMG3D_IPARAM_lag = 8
    MMG3D_IPARAM_optim = 9
    MMG3D_IPARAM_optimLES = 10
    MMG3D_IPARAM_noinsert = 11
    MMG3D_IPARAM_noswap = 12
    MMG3D_IPARAM_nomove = 13
    MMG3D_IPARAM_nosurf = 14
    MMG3D_IPARAM_nreg = 15
    MMG3D_IPARAM_xreg = 16
    MMG3D_IPARAM_numberOfLocalParam = 17
    MMG3D_IPARAM_numberOfLSBaseReferences = 18
    MMG3D_IPARAM_numberOfMat = 19
    MMG3D_IPARAM_numsubdomain = 20
    MMG3D_IPARAM_renum = 21
    MMG3D_IPARAM_anisosize = 22
    MMG3D_IPARAM_octree = 23
    MMG3D_IPARAM_nosizreq = 24
    MMG3D_IPARAM_isoref = 25
    MMG3D_DPARAM_angleDetection = 26
    MMG3D_DPARAM_hmin = 27
    MMG3D_DPARAM_hmax = 28
    MMG3D_DPARAM_hsiz = 29
    MMG3D_DPARAM_hausd = 30
    MMG3D_DPARAM_hgrad = 31
    MMG3D_DPARAM_hgradreq = 32
    MMG3D_DPARAM_ls = 33
    MMG3D_DPARAM_xreg = 34
    MMG3D_DPARAM_rmc = 35
    MMG3D_PARAM_size = 36

class MMG2D_Param(IntEnum):
    MMG2D_IPARAM_verbose = 0
    MMG2D_IPARAM_mem = 1
    MMG2D_IPARAM_debug = 2
    MMG2D_IPARAM_angle = 3
    MMG2D_IPARAM_iso = 4
    MMG2D_IPARAM_isosurf = 5
    MMG2D_IPARAM_opnbdy = 6
    MMG2D_IPARAM_lag = 7
    MMG2D_IPARAM_3dMedit = 8
    MMG2D_IPARAM_optim = 9
    MMG2D_IPARAM_noinsert = 10
    MMG2D_IPARAM_noswap = 11
    MMG2D_IPARAM_nomove = 12
    MMG2D_IPARAM_nosurf = 13
    MMG2D_IPARAM_nreg = 14
    MMG2D_IPARAM_xreg = 15
    MMG2D_IPARAM_numsubdomain = 16
    MMG2D_IPARAM_numberOfLocalParam = 17
    MMG2D_IPARAM_numberOfLSBaseReferences = 18
    MMG2D_IPARAM_numberOfMat = 19
    MMG2D_IPARAM_anisosize = 20
    MMG2D_IPARAM_nosizreq = 21
    MMG2D_DPARAM_angleDetection = 22
    MMG2D_DPARAM_hmin = 23
    MMG2D_DPARAM_hmax = 24
    MMG2D_DPARAM_hsiz = 25
    MMG2D_DPARAM_hausd = 26
    MMG2D_DPARAM_hgrad = 27
    MMG2D_DPARAM_hgradreq = 28
    MMG2D_DPARAM_ls = 29
    MMG2D_DPARAM_xreg = 30
    MMG2D_DPARAM_rmc = 31
    MMG2D_IPARAM_nofem = 32
    MMG2D_IPARAM_isoref = 33

class MMGS_Param(IntEnum):
    MMGS_IPARAM_verbose = 0
    MMGS_IPARAM_mem = 1
    MMGS_IPARAM_debug = 2
    MMGS_IPARAM_angle = 3
    MMGS_IPARAM_iso = 4
    MMGS_IPARAM_isosurf = 5
    MMGS_IPARAM_isoref = 6
    MMGS_IPARAM_keepRef = 7
    MMGS_IPARAM_optim = 8
    MMGS_IPARAM_noinsert = 9
    MMGS_IPARAM_noswap = 10
    MMGS_IPARAM_nomove = 11
    MMGS_IPARAM_nreg = 12
    MMGS_IPARAM_xreg = 13
    MMGS_IPARAM_numberOfLocalParam = 14
    MMGS_IPARAM_numberOfLSBaseReferences = 15
    MMGS_IPARAM_numberOfMat = 16
    MMGS_IPARAM_numsubdomain = 17
    MMGS_IPARAM_renum = 18
    MMGS_IPARAM_anisosize = 19
    MMGS_IPARAM_nosizreq = 20
    MMGS_DPARAM_angleDetection = 21
    MMGS_DPARAM_hmin = 22
    MMGS_DPARAM_hmax = 23
    MMGS_DPARAM_hsiz = 24
    MMGS_DPARAM_hausd = 25
    MMGS_DPARAM_hgrad = 26
    MMGS_DPARAM_hgradreq = 27
    MMGS_DPARAM_ls = 28
    MMGS_DPARAM_xreg = 29
    MMGS_DPARAM_rmc = 30
    MMGS_PARAM_size = 31

class MMG5_Par(ctypes.Structure):
    _fields_ = [("hmin",ctypes.c_double),
                ("hmax",ctypes.c_double),
                ("hausd",ctypes.c_double),
                ("ref",ctypes.c_int),
                ("elt",ctypes.c_int8)]

class MMG5_Point(ctypes.Structure):
    _fields_ = [("c",ctypes.c_double * 3),
                ("n",ctypes.c_double * 3),
                ("ref",ctypes.c_int),
                ("xp",ctypes.c_int),
                ("tmp",ctypes.c_int),
                ("flag",ctypes.c_int),
                ("s",ctypes.c_int),
                ("tag",ctypes.c_uint16),
                ("tagdel",ctypes.c_int8)]

class MMG5_xPoint(ctypes.Structure):
    _fields_ = [("n1",ctypes.c_double * 3),
                ("n2",ctypes.c_double * 3),
                ("nnor",ctypes.c_int8)]

class MMG5_Edge(ctypes.Structure):
    _fields_ = [("a",ctypes.c_int),
                ("b",ctypes.c_int),
                ("ref",ctypes.c_int),
                ("base",ctypes.c_int),
                ("tag",ctypes.c_uint16)]

class MMG5_Tria(ctypes.Structure):
    _fields_ = [("qual",ctypes.c_double),
                ("v",ctypes.c_int*3),
                ("ref",ctypes.c_int),
                ("base",ctypes.c_int),
                ("cc",ctypes.c_int),
                ("edg",ctypes.c_int*3),
                ("flag",ctypes.c_int),
                ("tag",ctypes.c_uint16 * 3)]

class MMG5_Quad(ctypes.Structure):
    _fields_ = [("v",ctypes.c_int*4),
                ("ref",ctypes.c_int),
                ("base",ctypes.c_int),
                ("edg",ctypes.c_int*4),
                ("tag",ctypes.c_uint16 * 4)]

class MMG5_Tetra(ctypes.Structure):
    _fields_ = [("qual",ctypes.c_double),
                ("v",ctypes.c_int),
                ("ref",ctypes.c_int),
                ("base",ctypes.c_int),
                ("mark",ctypes.c_int),
                ("xt",ctypes.c_int),
                ("flag",ctypes.c_int),
                ("tag",ctypes.c_uint16)]

class MMG5_xTetra(ctypes.Structure):
    _fields_ = [("ref",ctypes.c_int*4),
                ("edg",ctypes.c_int*6),
                ("ftag",ctypes.c_uint16 * 4),
                ("tag",ctypes.c_uint16 * 6),
                ("ori",ctypes.c_int8)]

class MMG5_Prism(ctypes.Structure):
    _fields_ = [("v",ctypes.c_int*6),
                ("ref",ctypes.c_int),
                ("base",ctypes.c_int),
                ("flag",ctypes.c_int),
                ("xpr",ctypes.c_int),
                ("tag",ctypes.c_int8)]

class MMG5_xPrism(ctypes.Structure):
    _fields_ = [("ref",ctypes.c_int*5),
                ("edg",ctypes.c_int*9),
                ("ftag",ctypes.c_uint16 * 5),
                ("tag",ctypes.c_uint16 * 9)]

class MMG5_Mat(ctypes.Structure):
    _fields_ = [("dospl",ctypes.c_int8),
                ("ref",ctypes.c_int),
                ("rin",ctypes.c_int),
                ("rex",ctypes.c_int)]

class MMG5_InvMat(ctypes.Structure):
    _fields_ = [("offset",ctypes.c_int),
                ("size",ctypes.c_int),
                ("lookup",ctypes.POINTER(ctypes.c_int))]

class MMG5_Info(ctypes.Structure):
    _fields_ = [("par",ctypes.POINTER(MMG5_Par)),
                ("dhd",ctypes.c_double),
                ("hmin",ctypes.c_double),
                ("hmax",ctypes.c_double),
                ("hsiz",ctypes.c_double),
                ("hgrad",ctypes.c_double),
                ("hgradreq",ctypes.c_double),
                ("hausd",ctypes.c_double),
                ("min",ctypes.c_double * 3),
                ("max",ctypes.c_double * 3),
                ("delta",ctypes.c_double),
                ("ls",ctypes.c_double),
                ("lxreg",ctypes.c_double),
                ("rmc",ctypes.c_double),
                ("br",ctypes.POINTER(ctypes.c_int)),
                ("isoref",ctypes.c_int),
                ("nsd",ctypes.c_int),
                ("mem",ctypes.c_int),
                ("npar",ctypes.c_int),
                ("npari",ctypes.c_int),
                ("nbr",ctypes.c_int),
                ("nbri",ctypes.c_int),
                ("opnbdy",ctypes.c_int),
                ("renum",ctypes.c_int),
                ("PROctree",ctypes.c_int),
                ("nmati",ctypes.c_int),
                ("nmat",ctypes.c_int),
                ("imprim",ctypes.c_int),
                ("nreg",ctypes.c_int8),
                ("xreg",ctypes.c_int8),
                ("ddebug",ctypes.c_int8),
                ("badkal",ctypes.c_int8),
                ("iso",ctypes.c_int8),
                ("isosurf",ctypes.c_int8),
                ("setfem",ctypes.c_int8),
                ("fem",ctypes.c_int8),
                ("lag",ctypes.c_int8),
                ("parTyp",ctypes.c_int8),
                ("sethmin",ctypes.c_int8),
                ("sethmax",ctypes.c_int8),
                ("ani",ctypes.c_uint8),
                ("optim",ctypes.c_uint8),
                ("optimLES",ctypes.c_uint8),
                ("noinsert",ctypes.c_uint8),
                ("noswap",ctypes.c_uint8),
                ("nomove",ctypes.c_uint8),
                ("nosurf",ctypes.c_uint8),
                ("nosizereq",ctypes.c_uint8),
                ("metRidTyp",ctypes.c_uint8),
                ("fparam",ctypes.c_char_p),
                ("mat",ctypes.POINTER(MMG5_Mat)),
                ("invmat",MMG5_InvMat)]

class MMG5_hgeom(ctypes.Structure):
    _fields_ = [("a",ctypes.c_int),
                ("b",ctypes.c_int),
                ("ref",ctypes.c_int),
                ("nxt",ctypes.c_int),
                ("tag",ctypes.c_uint16)]

class MMG5_HGeom(ctypes.Structure):
    _fields_ = [("geom",ctypes.POINTER(MMG5_hgeom)),
                ("siz",ctypes.c_int),
                ("max",ctypes.c_int),
                ("nxt",ctypes.c_int)]

class MMG5_hedge(ctypes.Structure):
    _fields_ = [("a",ctypes.c_int),
                ("b",ctypes.c_int),
                ("nxt",ctypes.c_int),
                ("k",ctypes.c_int),
                ("s",ctypes.c_int)]

class MMG5_Hash(ctypes.Structure):
    _fields_ = [("siz",ctypes.c_int),
                ("max",ctypes.c_int),
                ("nxt",ctypes.c_int),
                ("item",ctypes.POINTER(MMG5_hedge))]

class MMG5_Mesh(ctypes.Structure):
    _fields_ = [("memMax",ctypes.c_size_t),
                ("memCur",ctypes.c_size_t),
                ("gap",ctypes.c_double),
                ("ver",ctypes.c_int),
                ("dim",ctypes.c_int),
                ("type",ctypes.c_int),
                ("npi",ctypes.c_int),
                ("nti",ctypes.c_int),
                ("nai",ctypes.c_int),
                ("nei",ctypes.c_int),
                ("np",ctypes.c_int),
                ("na",ctypes.c_int),
                ("nt",ctypes.c_int),
                ("ne",ctypes.c_int),
                ("npmax",ctypes.c_int),
                ("namax",ctypes.c_int),
                ("ntmax",ctypes.c_int),
                ("nemax",ctypes.c_int),
                ("xpmax",ctypes.c_int),
                ("xtmax",ctypes.c_int),
                ("nquad",ctypes.c_int),
                ("nprism",ctypes.c_int),
                ("nsols",ctypes.c_int),
                ("nc1",ctypes.c_int),
                ("base",ctypes.c_int),
                ("mark",ctypes.c_int),
                ("xp",ctypes.c_int),
                ("xt",ctypes.c_int),
                ("xpr",ctypes.c_int),
                ("npnil",ctypes.c_int),
                ("nenil",ctypes.c_int),
                ("nanil",ctypes.c_int),
                ("adja",ctypes.POINTER(ctypes.c_int)),
                ("adjt",ctypes.POINTER(ctypes.c_int)),
                ("adjapr",ctypes.POINTER(ctypes.c_int)),
                ("adjq",ctypes.POINTER(ctypes.c_int)),
                ("ipar",ctypes.POINTER(ctypes.c_int)),
                ("point",ctypes.POINTER(MMG5_Point)),
                ("xpoint",ctypes.POINTER(MMG5_xPoint)),
                ("tetra",ctypes.POINTER(MMG5_Tetra)),
                ("xtetra",ctypes.POINTER(MMG5_xTetra)),
                ("prism",ctypes.POINTER(MMG5_Prism)),
                ("xprism",ctypes.POINTER(MMG5_xPrism)),
                ("tria",ctypes.POINTER(MMG5_Tria)),
                ("quadra",ctypes.POINTER(MMG5_Quad)),
                ("edge",ctypes.POINTER(MMG5_Edge)),
                ("htab",MMG5_HGeom),
                ("info",MMG5_Info),
                ("namein",ctypes.c_char_p),
                ("nameout",ctypes.c_char_p)]

class MMG5_Sol(ctypes.Structure):
    _fields_ = [("ver",ctypes.c_int),
                ("dim",ctypes.c_int),
                ("np",ctypes.c_int),
                ("npmax",ctypes.c_int),
                ("npi",ctypes.c_int),
                ("size",ctypes.c_int),
                ("type",ctypes.c_int),
                ("entities",ctypes.c_int),
                ("m",ctypes.POINTER(ctypes.c_double)),
                ("umin",ctypes.c_double),
                ("umax",ctypes.c_double),
                ("namein",ctypes.c_char_p),
                ("nameout",ctypes.c_char_p)]

    def __init__(self):
        self.size = 1

def MMG3D_Init_fileNames(mesh: MMG5_Mesh,sol: MMG5_Sol):
    lib3d.MMG3D_Init_fileNames(ctypes.byref(mesh),ctypes.byref(sol))

def MMG3D_Init_parameters(mesh: MMG5_Mesh):
    lib3d.MMG3D_Init_parameters(ctypes.byref(mesh))

def MMG3D_Set_inputMeshName(mesh: MMG5_Mesh,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib3d.MMG3D_Set_inputMeshName(ctypes.byref(mesh),name)
    return ier

def MMG3D_Set_outputMeshName(mesh: MMG5_Mesh,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib3d.MMG3D_Set_outputMeshName(ctypes.byref(mesh),name)
    return ier

def MMG3D_Set_inputSolName(mesh: MMG5_Mesh,sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib3d.MMG3D_Set_inputSolName(ctypes.byref(mesh),ctypes.byref(sol),name)
    return ier

def MMG3D_Set_outputSolName(mesh: MMG5_Mesh,sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib3d.MMG3D_Set_outputSolName(ctypes.byref(mesh),ctypes.byref(sol),name)
    return ier

def MMG3D_Set_inputParamName(mesh: MMG5_Mesh,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib3d.MMG3D_Set_inputParamName(ctypes.byref(mesh),name)
    return ier

def MMG3D_Set_solSize(mesh: MMG5_Mesh,sol: MMG5_Sol,typEntity: ctypes.c_int,np: ctypes.c_int,typSol: ctypes.c_int):
    ier = lib3d.MMG3D_Set_solSize(ctypes.byref(mesh),ctypes.byref(sol),typEntity,np,typSol)
    return ier

def MMG3D_Set_solsAtVerticesSize(mesh: MMG5_Mesh,*sol: MMG5_Sol,nsols: ctypes.c_int,nentities: ctypes.c_int,typSol):
    ier = lib3d.MMG3D_Set_solsAtVerticesSize(ctypes.byref(mesh),ctypes.byref(*sol),nsols,nentities,typSol.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier

def MMG3D_Set_meshSize(mesh: MMG5_Mesh,np: ctypes.c_int,ne: ctypes.c_int,nprism: ctypes.c_int,nt: ctypes.c_int,nquad: ctypes.c_int,na: ctypes.c_int):
    ier = lib3d.MMG3D_Set_meshSize(ctypes.byref(mesh),np,ne,nprism,nt,nquad,na)
    return ier

def MMG3D_Set_vertex(mesh: MMG5_Mesh,c0: ctypes.c_double,c1: ctypes.c_double,c2: ctypes.c_double,ref: ctypes.c_int,pos: ctypes.c_int):
    ier = lib3d.MMG3D_Set_vertex(ctypes.byref(mesh),c0,c1,c2,ref,pos)
    return ier

def MMG3D_Set_vertices(mesh: MMG5_Mesh,vertices,refs):
    ier = lib3d.MMG3D_Set_vertices(ctypes.byref(mesh),vertices.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),refs.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier

def MMG3D_Set_tetrahedron(mesh: MMG5_Mesh,v0: ctypes.c_int,v1: ctypes.c_int,v2: ctypes.c_int,v3: ctypes.c_int,ref: ctypes.c_int,pos: ctypes.c_int):
    ier = lib3d.MMG3D_Set_tetrahedron(ctypes.byref(mesh),v0,v1,v2,v3,ref,pos)
    return ier

def MMG3D_Set_tetrahedra(mesh: MMG5_Mesh,refs):
    ier = lib3d.MMG3D_Set_tetrahedra(ctypes.byref(mesh),tetra.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),refs.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier, tetra

def MMG3D_Set_prism(mesh: MMG5_Mesh,v0: ctypes.c_int,v1: ctypes.c_int,v2: ctypes.c_int,v3: ctypes.c_int,v4: ctypes.c_int,v5: ctypes.c_int,ref: ctypes.c_int,pos: ctypes.c_int):
    ier = lib3d.MMG3D_Set_prism(ctypes.byref(mesh),v0,v1,v2,v3,v4,v5,ref,pos)
    return ier

def MMG3D_Set_prisms(mesh: MMG5_Mesh,refs):
    ier = lib3d.MMG3D_Set_prisms(ctypes.byref(mesh),prisms.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),refs.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier, prisms

def MMG3D_Set_triangle(mesh: MMG5_Mesh,v0: ctypes.c_int,v1: ctypes.c_int,v2: ctypes.c_int,ref: ctypes.c_int,pos: ctypes.c_int):
    ier = lib3d.MMG3D_Set_triangle(ctypes.byref(mesh),v0,v1,v2,ref,pos)
    return ier

def MMG3D_Set_triangles(mesh: MMG5_Mesh,tria,refs):
    ier = lib3d.MMG3D_Set_triangles(ctypes.byref(mesh),tria.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),refs.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier

def MMG3D_Set_quadrilateral(mesh: MMG5_Mesh,v0: ctypes.c_int,v1: ctypes.c_int,v2: ctypes.c_int,v3: ctypes.c_int,ref: ctypes.c_int,pos: ctypes.c_int):
    ier = lib3d.MMG3D_Set_quadrilateral(ctypes.byref(mesh),v0,v1,v2,v3,ref,pos)
    return ier

def MMG3D_Set_quadrilaterals(mesh: MMG5_Mesh,refs):
    ier = lib3d.MMG3D_Set_quadrilaterals(ctypes.byref(mesh),quads.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),refs.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier, quads

def MMG3D_Set_edge(mesh: MMG5_Mesh,v0: ctypes.c_int,v1: ctypes.c_int,ref: ctypes.c_int,pos: ctypes.c_int):
    ier = lib3d.MMG3D_Set_edge(ctypes.byref(mesh),v0,v1,ref,pos)
    return ier

def MMG3D_Set_corner(mesh: MMG5_Mesh,k: ctypes.c_int):
    ier = lib3d.MMG3D_Set_corner(ctypes.byref(mesh),k)
    return ier

def MMG3D_Unset_corner(mesh: MMG5_Mesh,k: ctypes.c_int):
    ier = lib3d.MMG3D_Unset_corner(ctypes.byref(mesh),k)
    return ier

def MMG3D_Set_requiredVertex(mesh: MMG5_Mesh,k: ctypes.c_int):
    ier = lib3d.MMG3D_Set_requiredVertex(ctypes.byref(mesh),k)
    return ier

def MMG3D_Unset_requiredVertex(mesh: MMG5_Mesh,k: ctypes.c_int):
    ier = lib3d.MMG3D_Unset_requiredVertex(ctypes.byref(mesh),k)
    return ier

def MMG3D_Set_requiredTetrahedron(mesh: MMG5_Mesh,k: ctypes.c_int):
    ier = lib3d.MMG3D_Set_requiredTetrahedron(ctypes.byref(mesh),k)
    return ier

def MMG3D_Unset_requiredTetrahedron(mesh: MMG5_Mesh,k: ctypes.c_int):
    ier = lib3d.MMG3D_Unset_requiredTetrahedron(ctypes.byref(mesh),k)
    return ier

def MMG3D_Set_requiredTetrahedra(mesh: MMG5_Mesh,nreq: ctypes.c_int):
    ier = lib3d.MMG3D_Set_requiredTetrahedra(ctypes.byref(mesh),reqIdx.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),nreq)
    return ier, reqIdx

def MMG3D_Unset_requiredTetrahedra(mesh: MMG5_Mesh,nreq: ctypes.c_int):
    ier = lib3d.MMG3D_Unset_requiredTetrahedra(ctypes.byref(mesh),reqIdx.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),nreq)
    return ier, reqIdx

def MMG3D_Set_requiredTriangle(mesh: MMG5_Mesh,k: ctypes.c_int):
    ier = lib3d.MMG3D_Set_requiredTriangle(ctypes.byref(mesh),k)
    return ier

def MMG3D_Unset_requiredTriangle(mesh: MMG5_Mesh,k: ctypes.c_int):
    ier = lib3d.MMG3D_Unset_requiredTriangle(ctypes.byref(mesh),k)
    return ier

def MMG3D_Set_requiredTriangles(mesh: MMG5_Mesh,nreq: ctypes.c_int):
    ier = lib3d.MMG3D_Set_requiredTriangles(ctypes.byref(mesh),reqIdx.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),nreq)
    return ier, reqIdx

def MMG3D_Unset_requiredTriangles(mesh: MMG5_Mesh,nreq: ctypes.c_int):
    ier = lib3d.MMG3D_Unset_requiredTriangles(ctypes.byref(mesh),reqIdx.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),nreq)
    return ier, reqIdx

def MMG3D_Set_parallelTriangle(mesh: MMG5_Mesh,k: ctypes.c_int):
    ier = lib3d.MMG3D_Set_parallelTriangle(ctypes.byref(mesh),k)
    return ier

def MMG3D_Unset_parallelTriangle(mesh: MMG5_Mesh,k: ctypes.c_int):
    ier = lib3d.MMG3D_Unset_parallelTriangle(ctypes.byref(mesh),k)
    return ier

def MMG3D_Set_parallelTriangles(mesh: MMG5_Mesh,npar: ctypes.c_int):
    ier = lib3d.MMG3D_Set_parallelTriangles(ctypes.byref(mesh),parIdx.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),npar)
    return ier, parIdx

def MMG3D_Unset_parallelTriangles(mesh: MMG5_Mesh,npar: ctypes.c_int):
    ier = lib3d.MMG3D_Unset_parallelTriangles(ctypes.byref(mesh),parIdx.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),npar)
    return ier, parIdx

def MMG3D_Set_ridge(mesh: MMG5_Mesh,k: ctypes.c_int):
    ier = lib3d.MMG3D_Set_ridge(ctypes.byref(mesh),k)
    return ier

def MMG3D_Unset_ridge(mesh: MMG5_Mesh,k: ctypes.c_int):
    ier = lib3d.MMG3D_Unset_ridge(ctypes.byref(mesh),k)
    return ier

def MMG3D_Set_requiredEdge(mesh: MMG5_Mesh,k: ctypes.c_int):
    ier = lib3d.MMG3D_Set_requiredEdge(ctypes.byref(mesh),k)
    return ier

def MMG3D_Unset_requiredEdge(mesh: MMG5_Mesh,k: ctypes.c_int):
    ier = lib3d.MMG3D_Unset_requiredEdge(ctypes.byref(mesh),k)
    return ier

def MMG3D_Set_normalAtVertex(mesh: MMG5_Mesh,k: ctypes.c_int,n0: ctypes.c_double,n1: ctypes.c_double,n2: ctypes.c_double):
    ier = lib3d.MMG3D_Set_normalAtVertex(ctypes.byref(mesh),k,n0,n1,n2)
    return ier

def MMG3D_Set_scalarSol(met: MMG5_Sol,s: ctypes.c_double,pos: ctypes.c_int):
    ier = lib3d.MMG3D_Set_scalarSol(ctypes.byref(met),s,pos)
    return ier

def MMG3D_Set_scalarSols(met: MMG5_Sol,s):
    ier = lib3d.MMG3D_Set_scalarSols(ctypes.byref(met),s.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    return ier

def MMG3D_Set_vectorSol(met: MMG5_Sol,vx: ctypes.c_double,vy: ctypes.c_double,vz: ctypes.c_double,pos: ctypes.c_int):
    ier = lib3d.MMG3D_Set_vectorSol(ctypes.byref(met),vx,vy,vz,pos)
    return ier

def MMG3D_Set_vectorSols(met: MMG5_Sol,sols):
    ier = lib3d.MMG3D_Set_vectorSols(ctypes.byref(met),sols.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    return ier

def MMG3D_Set_tensorSol(met: MMG5_Sol,m11: ctypes.c_double,m12: ctypes.c_double,m13: ctypes.c_double,m22: ctypes.c_double,m23: ctypes.c_double,m33: ctypes.c_double,pos: ctypes.c_int):
    ier = lib3d.MMG3D_Set_tensorSol(ctypes.byref(met),m11,m12,m13,m22,m23,m33,pos)
    return ier

def MMG3D_Set_tensorSols(met: MMG5_Sol,sols):
    ier = lib3d.MMG3D_Set_tensorSols(ctypes.byref(met),sols.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    return ier

def MMG3D_Set_ithSol_inSolsAtVertices(sol: MMG5_Sol,i: ctypes.c_int,s,pos: ctypes.c_int):
    ier = lib3d.MMG3D_Set_ithSol_inSolsAtVertices(ctypes.byref(sol),i,s.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),pos)
    return ier

def MMG3D_Set_ithSols_inSolsAtVertices(sol: MMG5_Sol,i: ctypes.c_int,s):
    ier = lib3d.MMG3D_Set_ithSols_inSolsAtVertices(ctypes.byref(sol),i,s.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    return ier

def MMG3D_Set_handGivenMesh(mesh: MMG5_Mesh):
    lib3d.MMG3D_Set_handGivenMesh(ctypes.byref(mesh))

def MMG3D_Chk_meshData(mesh: MMG5_Mesh,met: MMG5_Sol):
    ier = lib3d.MMG3D_Chk_meshData(ctypes.byref(mesh),ctypes.byref(met))
    return ier

def MMG3D_Set_iparameter(mesh: MMG5_Mesh,sol: MMG5_Sol,iparam: ctypes.c_int,val: ctypes.c_int):
    ier = lib3d.MMG3D_Set_iparameter(ctypes.byref(mesh),ctypes.byref(sol),iparam,val)
    return ier

def MMG3D_Set_dparameter(mesh: MMG5_Mesh,sol: MMG5_Sol,dparam: ctypes.c_int,val: ctypes.c_double):
    ier = lib3d.MMG3D_Set_dparameter(ctypes.byref(mesh),ctypes.byref(sol),dparam,val)
    return ier

def MMG3D_Set_localParameter(mesh: MMG5_Mesh,sol: MMG5_Sol,typ: ctypes.c_int,ref: ctypes.c_int,hmin: ctypes.c_double,hmax: ctypes.c_double,hausd: ctypes.c_double):
    ier = lib3d.MMG3D_Set_localParameter(ctypes.byref(mesh),ctypes.byref(sol),typ,ref,hmin,hmax,hausd)
    return ier

def MMG3D_Set_multiMat(mesh: MMG5_Mesh,sol: MMG5_Sol,ref: ctypes.c_int,split: ctypes.c_int,rmin: ctypes.c_int,rplus: ctypes.c_int):
    ier = lib3d.MMG3D_Set_multiMat(ctypes.byref(mesh),ctypes.byref(sol),ref,split,rmin,rplus)
    return ier

def MMG3D_Set_lsBaseReference(mesh: MMG5_Mesh,sol: MMG5_Sol,br: ctypes.c_int):
    ier = lib3d.MMG3D_Set_lsBaseReference(ctypes.byref(mesh),ctypes.byref(sol),br)
    return ier

def MMG3D_Get_meshSize(mesh: MMG5_Mesh):
    ier = lib3d.MMG3D_Get_meshSize(ctypes.byref(mesh),np.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),ne.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),nprism.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),nt.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),nquad.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),na.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier, np, ne, nprism, nt, nquad, na

def MMG3D_Get_solSize(mesh: MMG5_Mesh,sol: MMG5_Sol,typSol):
    ier = lib3d.MMG3D_Get_solSize(ctypes.byref(mesh),ctypes.byref(sol),typEntity.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),np.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),typSol.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier, typEntity, np

def MMG3D_Get_solsAtVerticesSize(mesh: MMG5_Mesh,sol: MMG5_Sol,typSol):
    ier = lib3d.MMG3D_Get_solsAtVerticesSize(ctypes.byref(mesh),ctypes.byref(sol),nsols.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),nentities.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),typSol.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier, nsols, nentities

def MMG3D_Get_vertex(mesh: MMG5_Mesh):
    ier = lib3d.MMG3D_Get_vertex(ctypes.byref(mesh),c0.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),c1.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),c2.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),ref.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),isCorner.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),isRequired.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier, c0, c1, c2, ref, isCorner, isRequired

def MMG3D_GetByIdx_vertex(mesh: MMG5_Mesh,idx: ctypes.c_int):
    ier = lib3d.MMG3D_GetByIdx_vertex(ctypes.byref(mesh),c0.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),c1.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),c2.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),ref.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),isCorner.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),isRequired.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),idx)
    return ier, c0, c1, c2, ref, isCorner, isRequired

def MMG3D_Get_vertices(mesh: MMG5_Mesh,vertices,refs,areCorners,areRequired):
    ier = lib3d.MMG3D_Get_vertices(ctypes.byref(mesh),vertices.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),refs.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),areCorners.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),areRequired.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier

def MMG3D_Get_tetrahedron(mesh: MMG5_Mesh):
    ier = lib3d.MMG3D_Get_tetrahedron(ctypes.byref(mesh),v0.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),v1.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),v2.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),v3.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),ref.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),isRequired.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier, v0, v1, v2, v3, ref, isRequired

def MMG3D_Get_tetrahedra(mesh: MMG5_Mesh,refs,areRequired):
    ier = lib3d.MMG3D_Get_tetrahedra(ctypes.byref(mesh),tetra.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),refs.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),areRequired.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier, tetra

def MMG3D_Get_prism(mesh: MMG5_Mesh):
    ier = lib3d.MMG3D_Get_prism(ctypes.byref(mesh),v0.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),v1.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),v2.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),v3.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),v4.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),v5.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),ref.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),isRequired.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier, v0, v1, v2, v3, v4, v5, ref, isRequired

def MMG3D_Get_prisms(mesh: MMG5_Mesh,refs,areRequired):
    ier = lib3d.MMG3D_Get_prisms(ctypes.byref(mesh),prisms.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),refs.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),areRequired.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier, prisms

def MMG3D_Get_triangle(mesh: MMG5_Mesh):
    ier = lib3d.MMG3D_Get_triangle(ctypes.byref(mesh),v0.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),v1.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),v2.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),ref.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),isRequired.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier, v0, v1, v2, ref, isRequired

def MMG3D_Get_triangles(mesh: MMG5_Mesh,tria,refs,areRequired):
    ier = lib3d.MMG3D_Get_triangles(ctypes.byref(mesh),tria.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),refs.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),areRequired.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier

def MMG3D_Get_quadrilateral(mesh: MMG5_Mesh):
    ier = lib3d.MMG3D_Get_quadrilateral(ctypes.byref(mesh),v0.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),v1.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),v2.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),v3.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),ref.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),isRequired.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier, v0, v1, v2, v3, ref, isRequired

def MMG3D_Get_quadrilaterals(mesh: MMG5_Mesh,refs,areRequired):
    ier = lib3d.MMG3D_Get_quadrilaterals(ctypes.byref(mesh),quads.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),refs.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),areRequired.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier, quads

def MMG3D_Get_edge(mesh: MMG5_Mesh):
    ier = lib3d.MMG3D_Get_edge(ctypes.byref(mesh),e0.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),e1.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),ref.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),isRidge.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),isRequired.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier, e0, e1, ref, isRidge, isRequired

def MMG3D_Set_edges(mesh: MMG5_Mesh,edges,refs):
    ier = lib3d.MMG3D_Set_edges(ctypes.byref(mesh),edges.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),refs.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier

def MMG3D_Get_edges(mesh: MMG5_Mesh,edges,refs,areRidges,areRequired):
    ier = lib3d.MMG3D_Get_edges(ctypes.byref(mesh),edges.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),refs.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),areRidges.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),areRequired.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier

def MMG3D_Get_normalAtVertex(mesh: MMG5_Mesh,k: ctypes.c_int):
    ier = lib3d.MMG3D_Get_normalAtVertex(ctypes.byref(mesh),k,n0.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),n1.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),n2.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    return ier, n0, n1, n2

def MMG3D_Get_tetrahedronQuality(mesh: MMG5_Mesh,met: MMG5_Sol,k: ctypes.c_int):
    ier = lib3d.MMG3D_Get_tetrahedronQuality(ctypes.byref(mesh),ctypes.byref(met),k)
    return ier

def MMG3D_Get_scalarSol(met: MMG5_Sol,s):
    ier = lib3d.MMG3D_Get_scalarSol(ctypes.byref(met),s.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    return ier

def MMG3D_Get_scalarSols(met: MMG5_Sol,s):
    ier = lib3d.MMG3D_Get_scalarSols(ctypes.byref(met),s.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    return ier

def MMG3D_Get_vectorSol(met: MMG5_Sol):
    ier = lib3d.MMG3D_Get_vectorSol(ctypes.byref(met),vx.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),vy.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),vz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    return ier, vx, vy, vz

def MMG3D_Get_vectorSols(met: MMG5_Sol,sols):
    ier = lib3d.MMG3D_Get_vectorSols(ctypes.byref(met),sols.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    return ier

def MMG3D_Get_tensorSol(met: MMG5_Sol):
    ier = lib3d.MMG3D_Get_tensorSol(ctypes.byref(met),m11.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),m12.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),m13.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),m22.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),m23.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),m33.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    return ier, m11, m12, m13, m22, m23, m33

def MMG3D_Get_tensorSols(met: MMG5_Sol,sols):
    ier = lib3d.MMG3D_Get_tensorSols(ctypes.byref(met),sols.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    return ier

def MMG3D_Get_ithSol_inSolsAtVertices(sol: MMG5_Sol,i: ctypes.c_int,s,pos: ctypes.c_int):
    ier = lib3d.MMG3D_Get_ithSol_inSolsAtVertices(ctypes.byref(sol),i,s.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),pos)
    return ier

def MMG3D_Get_ithSols_inSolsAtVertices(sol: MMG5_Sol,i: ctypes.c_int,s):
    ier = lib3d.MMG3D_Get_ithSols_inSolsAtVertices(ctypes.byref(sol),i,s.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    return ier

def MMG3D_Get_iparameter(mesh: MMG5_Mesh,iparam: ctypes.c_int):
    ier = lib3d.MMG3D_Get_iparameter(ctypes.byref(mesh),iparam)
    return ier

def MMG3D_Add_tetrahedron(mesh: MMG5_Mesh,v0: ctypes.c_int,v1: ctypes.c_int,v2: ctypes.c_int,v3: ctypes.c_int,ref: ctypes.c_int):
    ier = lib3d.MMG3D_Add_tetrahedron(ctypes.byref(mesh),v0,v1,v2,v3,ref)
    return ier

def MMG3D_Add_vertex(mesh: MMG5_Mesh,c0: ctypes.c_double,c1: ctypes.c_double,c2: ctypes.c_double,ref: ctypes.c_int):
    ier = lib3d.MMG3D_Add_vertex(ctypes.byref(mesh),c0,c1,c2,ref)
    return ier

def MMG3D_loadMesh(mesh: MMG5_Mesh,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib3d.MMG3D_loadMesh(ctypes.byref(mesh),name)
    return ier

def MMG3D_loadMshMesh(mesh: MMG5_Mesh,sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib3d.MMG3D_loadMshMesh(ctypes.byref(mesh),ctypes.byref(sol),name)
    return ier

def MMG3D_loadVtuMesh(mesh: MMG5_Mesh,met: MMG5_Sol,sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib3d.MMG3D_loadVtuMesh(ctypes.byref(mesh),ctypes.byref(met),ctypes.byref(sol),name)
    return ier

def MMG3D_loadVtuMesh_and_allData(mesh: MMG5_Mesh,*sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib3d.MMG3D_loadVtuMesh_and_allData(ctypes.byref(mesh),ctypes.byref(*sol),name)
    return ier

def MMG3D_loadVtkMesh(mesh: MMG5_Mesh,met: MMG5_Sol,sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib3d.MMG3D_loadVtkMesh(ctypes.byref(mesh),ctypes.byref(met),ctypes.byref(sol),name)
    return ier

def MMG3D_loadVtkMesh_and_allData(mesh: MMG5_Mesh,*sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib3d.MMG3D_loadVtkMesh_and_allData(ctypes.byref(mesh),ctypes.byref(*sol),name)
    return ier

def MMG3D_loadMshMesh_and_allData(mesh: MMG5_Mesh,*sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib3d.MMG3D_loadMshMesh_and_allData(ctypes.byref(mesh),ctypes.byref(*sol),name)
    return ier

def MMG3D_loadGenericMesh(mesh: MMG5_Mesh,met: MMG5_Sol,sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib3d.MMG3D_loadGenericMesh(ctypes.byref(mesh),ctypes.byref(met),ctypes.byref(sol),name)
    return ier

def MMG3D_saveMesh(mesh: MMG5_Mesh,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib3d.MMG3D_saveMesh(ctypes.byref(mesh),name)
    return ier

def MMG3D_saveMshMesh(mesh: MMG5_Mesh,sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib3d.MMG3D_saveMshMesh(ctypes.byref(mesh),ctypes.byref(sol),name)
    return ier

def MMG3D_saveMshMesh_and_allData(mesh: MMG5_Mesh,*sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib3d.MMG3D_saveMshMesh_and_allData(ctypes.byref(mesh),ctypes.byref(*sol),name)
    return ier

def MMG3D_saveVtkMesh(mesh: MMG5_Mesh,sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib3d.MMG3D_saveVtkMesh(ctypes.byref(mesh),ctypes.byref(sol),name)
    return ier

def MMG3D_saveVtkMesh_and_allData(mesh: MMG5_Mesh,*sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib3d.MMG3D_saveVtkMesh_and_allData(ctypes.byref(mesh),ctypes.byref(*sol),name)
    return ier

def MMG3D_saveVtuMesh(mesh: MMG5_Mesh,sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib3d.MMG3D_saveVtuMesh(ctypes.byref(mesh),ctypes.byref(sol),name)
    return ier

def MMG3D_saveVtuMesh_and_allData(mesh: MMG5_Mesh,*sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib3d.MMG3D_saveVtuMesh_and_allData(ctypes.byref(mesh),ctypes.byref(*sol),name)
    return ier

def MMG3D_saveTetgenMesh(mesh: MMG5_Mesh,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib3d.MMG3D_saveTetgenMesh(ctypes.byref(mesh),name)
    return ier

def MMG3D_saveGenericMesh(mesh: MMG5_Mesh,sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib3d.MMG3D_saveGenericMesh(ctypes.byref(mesh),ctypes.byref(sol),name)
    return ier

def MMG3D_loadSol(mesh: MMG5_Mesh,met: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib3d.MMG3D_loadSol(ctypes.byref(mesh),ctypes.byref(met),name)
    return ier

def MMG3D_loadAllSols(mesh: MMG5_Mesh,*sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib3d.MMG3D_loadAllSols(ctypes.byref(mesh),ctypes.byref(*sol),name)
    return ier

def MMG3D_saveSol(mesh: MMG5_Mesh,met: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib3d.MMG3D_saveSol(ctypes.byref(mesh),ctypes.byref(met),name)
    return ier

def MMG3D_saveAllSols(mesh: MMG5_Mesh,*sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib3d.MMG3D_saveAllSols(ctypes.byref(mesh),ctypes.byref(*sol),name)
    return ier

def MMG3D_Free_allSols(mesh: MMG5_Mesh,*sol: MMG5_Sol):
    ier = lib3d.MMG3D_Free_allSols(ctypes.byref(mesh),ctypes.byref(*sol))
    return ier

def MMG3D_mmg3dlib(mesh: MMG5_Mesh,met: MMG5_Sol):
    ier = lib3d.MMG3D_mmg3dlib(ctypes.byref(mesh),ctypes.byref(met))
    return ier

def MMG3D_mmg3dls(mesh: MMG5_Mesh,sol: MMG5_Sol,met: MMG5_Sol):
    ier = lib3d.MMG3D_mmg3dls(ctypes.byref(mesh),ctypes.byref(sol),ctypes.byref(met))
    return ier

def MMG3D_mmg3dmov(mesh: MMG5_Mesh,met: MMG5_Sol,disp: MMG5_Sol):
    ier = lib3d.MMG3D_mmg3dmov(ctypes.byref(mesh),ctypes.byref(met),ctypes.byref(disp))
    return ier

def MMG3D_defaultValues(mesh: MMG5_Mesh):
    ier = lib3d.MMG3D_defaultValues(ctypes.byref(mesh))
    return ier

def MMG3D_parsop(mesh: MMG5_Mesh,met: MMG5_Sol):
    ier = lib3d.MMG3D_parsop(ctypes.byref(mesh),ctypes.byref(met))
    return ier

def MMG3D_usage(name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib3d.MMG3D_usage(name)
    return ier

def MMG3D_stockOptions(mesh: MMG5_Mesh):
    ier = lib3d.MMG3D_stockOptions(ctypes.byref(mesh),info.ctypes.data_as(ctypes.POINTER(MMG5_Info)))
    return ier, info

def MMG3D_destockOptions(mesh: MMG5_Mesh):
    lib3d.MMG3D_destockOptions(ctypes.byref(mesh),info.ctypes.data_as(ctypes.POINTER(MMG5_Info)))

def MMG3D_mmg3dcheck(mesh: MMG5_Mesh,met: MMG5_Sol,sol: MMG5_Sol,critmin: ctypes.c_double,lmin: ctypes.c_double,lmax: ctypes.c_double,metRidTyp: ctypes.c_int8):
    ier = lib3d.MMG3D_mmg3dcheck(ctypes.byref(mesh),ctypes.byref(met),ctypes.byref(sol),critmin,lmin,lmax,eltab.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),metRidTyp)
    return ier, eltab

def MMG3D_searchqua(mesh: MMG5_Mesh,met: MMG5_Sol,critmin: ctypes.c_double,metRidTyp: ctypes.c_int8):
    lib3d.MMG3D_searchqua(ctypes.byref(mesh),ctypes.byref(met),critmin,eltab.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),metRidTyp)

def MMG3D_searchlen(mesh: MMG5_Mesh,met: MMG5_Sol,lmin: ctypes.c_double,lmax: ctypes.c_double,metRidTyp: ctypes.c_int8):
    ier = lib3d.MMG3D_searchlen(ctypes.byref(mesh),ctypes.byref(met),lmin,lmax,eltab.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),metRidTyp)
    return ier, eltab

def MMG3D_Get_adjaTet(mesh: MMG5_Mesh,kel: ctypes.c_int,listet: ctypes.c_int*4):
    ier = lib3d.MMG3D_Get_adjaTet(ctypes.byref(mesh),kel,listet)
    return ier

def MMG3D_hashTetra(mesh: MMG5_Mesh,pack: ctypes.c_int):
    ier = lib3d.MMG3D_hashTetra(ctypes.byref(mesh),pack)
    return ier

def MMG3D_Set_constantSize(mesh: MMG5_Mesh,met: MMG5_Sol):
    ier = lib3d.MMG3D_Set_constantSize(ctypes.byref(mesh),ctypes.byref(met))
    return ier

def MMG3D_switch_metricStorage(mesh: MMG5_Mesh,met: MMG5_Sol):
    ier = lib3d.MMG3D_switch_metricStorage(ctypes.byref(mesh),ctypes.byref(met))
    return ier

def MMG3D_setfunc(mesh: MMG5_Mesh,met: MMG5_Sol):
    lib3d.MMG3D_setfunc(ctypes.byref(mesh),ctypes.byref(met))

def MMG3D_Get_numberOfNonBdyTriangles(mesh: MMG5_Mesh):
    ier = lib3d.MMG3D_Get_numberOfNonBdyTriangles(ctypes.byref(mesh),nb_tria.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier, nb_tria

def MMG3D_Get_nonBdyTriangle(mesh: MMG5_Mesh,idx: ctypes.c_int):
    ier = lib3d.MMG3D_Get_nonBdyTriangle(ctypes.byref(mesh),v0.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),v1.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),v2.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),ref.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),idx)
    return ier, v0, v1, v2, ref

def MMG3D_Get_tetFromTria(mesh: MMG5_Mesh,ktri: ctypes.c_int):
    ier = lib3d.MMG3D_Get_tetFromTria(ctypes.byref(mesh),ktri,ktet.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),iface.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier, ktet, iface

def MMG3D_Get_tetsFromTria(mesh: MMG5_Mesh,ktri: ctypes.c_int,ktet: ctypes.c_int*2,iface: ctypes.c_int*2):
    ier = lib3d.MMG3D_Get_tetsFromTria(ctypes.byref(mesh),ktri,ktet,iface)
    return ier

def MMG3D_Compute_eigenv(m: ctypes.c_double*6,lambda0: ctypes.c_double*3,vp: ctypes.c_double*3):
    ier = lib3d.MMG3D_Compute_eigenv(m,lambda0,vp)
    return ier

def MMG3D_Clean_isoSurf(mesh: MMG5_Mesh):
    ier = lib3d.MMG3D_Clean_isoSurf(ctypes.byref(mesh))
    return ier

def MMG3D_Free_solutions(mesh: MMG5_Mesh,sol: MMG5_Sol):
    lib3d.MMG3D_Free_solutions(ctypes.byref(mesh),ctypes.byref(sol))

def MMG2D_Init_fileNames(mesh: MMG5_Mesh,sol: MMG5_Sol):
    lib2d.MMG2D_Init_fileNames(ctypes.byref(mesh),ctypes.byref(sol))

def MMG2D_Init_parameters(mesh: MMG5_Mesh):
    lib2d.MMG2D_Init_parameters(ctypes.byref(mesh))

def MMG2D_Set_inputMeshName(mesh: MMG5_Mesh,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_Set_inputMeshName(ctypes.byref(mesh),name)
    return ier

def MMG2D_Set_outputMeshName(mesh: MMG5_Mesh,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_Set_outputMeshName(ctypes.byref(mesh),name)
    return ier

def MMG2D_Set_inputSolName(mesh: MMG5_Mesh,sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_Set_inputSolName(ctypes.byref(mesh),ctypes.byref(sol),name)
    return ier

def MMG2D_Set_outputSolName(mesh: MMG5_Mesh,sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_Set_outputSolName(ctypes.byref(mesh),ctypes.byref(sol),name)
    return ier

def MMG2D_Set_inputParamName(mesh: MMG5_Mesh,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_Set_inputParamName(ctypes.byref(mesh),name)
    return ier

def MMG2D_Set_iparameter(mesh: MMG5_Mesh,sol: MMG5_Sol,iparam: ctypes.c_int,val: ctypes.c_int):
    ier = lib2d.MMG2D_Set_iparameter(ctypes.byref(mesh),ctypes.byref(sol),iparam,val)
    return ier

def MMG2D_Set_dparameter(mesh: MMG5_Mesh,sol: MMG5_Sol,dparam: ctypes.c_int,val: ctypes.c_double):
    ier = lib2d.MMG2D_Set_dparameter(ctypes.byref(mesh),ctypes.byref(sol),dparam,val)
    return ier

def MMG2D_Set_localParameter(mesh: MMG5_Mesh,sol: MMG5_Sol,typ: ctypes.c_int,ref: ctypes.c_int,hmin: ctypes.c_double,hmax: ctypes.c_double,hausd: ctypes.c_double):
    ier = lib2d.MMG2D_Set_localParameter(ctypes.byref(mesh),ctypes.byref(sol),typ,ref,hmin,hmax,hausd)
    return ier

def MMG2D_Set_multiMat(mesh: MMG5_Mesh,sol: MMG5_Sol,ref: ctypes.c_int,split: ctypes.c_int,rmin: ctypes.c_int,rplus: ctypes.c_int):
    ier = lib2d.MMG2D_Set_multiMat(ctypes.byref(mesh),ctypes.byref(sol),ref,split,rmin,rplus)
    return ier

def MMG2D_Set_lsBaseReference(mesh: MMG5_Mesh,sol: MMG5_Sol,br: ctypes.c_int):
    ier = lib2d.MMG2D_Set_lsBaseReference(ctypes.byref(mesh),ctypes.byref(sol),br)
    return ier

def MMG2D_Set_meshSize(mesh: MMG5_Mesh,np: ctypes.c_int,nt: ctypes.c_int,nquad: ctypes.c_int,na: ctypes.c_int):
    ier = lib2d.MMG2D_Set_meshSize(ctypes.byref(mesh),np,nt,nquad,na)
    return ier

def MMG2D_Set_solSize(mesh: MMG5_Mesh,sol: MMG5_Sol,typEntity: ctypes.c_int,np: ctypes.c_int,typSol: ctypes.c_int):
    ier = lib2d.MMG2D_Set_solSize(ctypes.byref(mesh),ctypes.byref(sol),typEntity,np,typSol)
    return ier

def MMG2D_Set_solsAtVerticesSize(mesh: MMG5_Mesh,*sol: MMG5_Sol,nsols: ctypes.c_int,nentities: ctypes.c_int,typSol):
    ier = lib2d.MMG2D_Set_solsAtVerticesSize(ctypes.byref(mesh),ctypes.byref(*sol),nsols,nentities,typSol.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier

def MMG2D_Set_vertex(mesh: MMG5_Mesh,c0: ctypes.c_double,c1: ctypes.c_double,ref: ctypes.c_int,pos: ctypes.c_int):
    ier = lib2d.MMG2D_Set_vertex(ctypes.byref(mesh),c0,c1,ref,pos)
    return ier

def MMG2D_Set_vertices(mesh: MMG5_Mesh,vertices,refs):
    ier = lib2d.MMG2D_Set_vertices(ctypes.byref(mesh),vertices.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),refs.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier

def MMG2D_Set_corner(mesh: MMG5_Mesh,k: ctypes.c_int):
    ier = lib2d.MMG2D_Set_corner(ctypes.byref(mesh),k)
    return ier

def MMG2D_Unset_corner(mesh: MMG5_Mesh,k: ctypes.c_int):
    ier = lib2d.MMG2D_Unset_corner(ctypes.byref(mesh),k)
    return ier

def MMG2D_Set_requiredVertex(mesh: MMG5_Mesh,k: ctypes.c_int):
    ier = lib2d.MMG2D_Set_requiredVertex(ctypes.byref(mesh),k)
    return ier

def MMG2D_Unset_requiredVertex(mesh: MMG5_Mesh,k: ctypes.c_int):
    ier = lib2d.MMG2D_Unset_requiredVertex(ctypes.byref(mesh),k)
    return ier

def MMG2D_Set_triangle(mesh: MMG5_Mesh,v0: ctypes.c_int,v1: ctypes.c_int,v2: ctypes.c_int,ref: ctypes.c_int,pos: ctypes.c_int):
    ier = lib2d.MMG2D_Set_triangle(ctypes.byref(mesh),v0,v1,v2,ref,pos)
    return ier

def MMG2D_Set_triangles(mesh: MMG5_Mesh,tria,refs):
    ier = lib2d.MMG2D_Set_triangles(ctypes.byref(mesh),tria.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),refs.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier

def MMG2D_Set_requiredTriangle(mesh: MMG5_Mesh,k: ctypes.c_int):
    ier = lib2d.MMG2D_Set_requiredTriangle(ctypes.byref(mesh),k)
    return ier

def MMG2D_Unset_requiredTriangle(mesh: MMG5_Mesh,k: ctypes.c_int):
    ier = lib2d.MMG2D_Unset_requiredTriangle(ctypes.byref(mesh),k)
    return ier

def MMG2D_Set_quadrilateral(mesh: MMG5_Mesh,v0: ctypes.c_int,v1: ctypes.c_int,v2: ctypes.c_int,v3: ctypes.c_int,ref: ctypes.c_int,pos: ctypes.c_int):
    ier = lib2d.MMG2D_Set_quadrilateral(ctypes.byref(mesh),v0,v1,v2,v3,ref,pos)
    return ier

def MMG2D_Set_quadrilaterals(mesh: MMG5_Mesh,quadra,refs):
    ier = lib2d.MMG2D_Set_quadrilaterals(ctypes.byref(mesh),quadra.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),refs.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier

def MMG2D_Set_edge(mesh: MMG5_Mesh,v0: ctypes.c_int,v1: ctypes.c_int,ref: ctypes.c_int,pos: ctypes.c_int):
    ier = lib2d.MMG2D_Set_edge(ctypes.byref(mesh),v0,v1,ref,pos)
    return ier

def MMG2D_Set_edges(mesh: MMG5_Mesh,edges,refs):
    ier = lib2d.MMG2D_Set_edges(ctypes.byref(mesh),edges.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),refs.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier

def MMG2D_Set_requiredEdge(mesh: MMG5_Mesh,k: ctypes.c_int):
    ier = lib2d.MMG2D_Set_requiredEdge(ctypes.byref(mesh),k)
    return ier

def MMG2D_Unset_requiredEdge(mesh: MMG5_Mesh,k: ctypes.c_int):
    ier = lib2d.MMG2D_Unset_requiredEdge(ctypes.byref(mesh),k)
    return ier

def MMG2D_Set_parallelEdge(mesh: MMG5_Mesh,k: ctypes.c_int):
    ier = lib2d.MMG2D_Set_parallelEdge(ctypes.byref(mesh),k)
    return ier

def MMG2D_Set_scalarSol(met: MMG5_Sol,s: ctypes.c_double,pos: ctypes.c_int):
    ier = lib2d.MMG2D_Set_scalarSol(ctypes.byref(met),s,pos)
    return ier

def MMG2D_Set_scalarSols(met: MMG5_Sol,s):
    ier = lib2d.MMG2D_Set_scalarSols(ctypes.byref(met),s.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    return ier

def MMG2D_Set_vectorSol(met: MMG5_Sol,vx: ctypes.c_double,vy: ctypes.c_double,pos: ctypes.c_int):
    ier = lib2d.MMG2D_Set_vectorSol(ctypes.byref(met),vx,vy,pos)
    return ier

def MMG2D_Set_vectorSols(met: MMG5_Sol,sols):
    ier = lib2d.MMG2D_Set_vectorSols(ctypes.byref(met),sols.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    return ier

def MMG2D_Set_tensorSol(met: MMG5_Sol,m11: ctypes.c_double,m12: ctypes.c_double,m22: ctypes.c_double,pos: ctypes.c_int):
    ier = lib2d.MMG2D_Set_tensorSol(ctypes.byref(met),m11,m12,m22,pos)
    return ier

def MMG2D_Set_tensorSols(met: MMG5_Sol,sols):
    ier = lib2d.MMG2D_Set_tensorSols(ctypes.byref(met),sols.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    return ier

def MMG2D_Set_ithSol_inSolsAtVertices(sol: MMG5_Sol,i: ctypes.c_int,s,pos: ctypes.c_int):
    ier = lib2d.MMG2D_Set_ithSol_inSolsAtVertices(ctypes.byref(sol),i,s.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),pos)
    return ier

def MMG2D_Set_ithSols_inSolsAtVertices(sol: MMG5_Sol,i: ctypes.c_int,s):
    ier = lib2d.MMG2D_Set_ithSols_inSolsAtVertices(ctypes.byref(sol),i,s.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    return ier

def MMG2D_Get_meshSize(mesh: MMG5_Mesh):
    ier = lib2d.MMG2D_Get_meshSize(ctypes.byref(mesh),np.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),nt.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),nquad.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),na.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier, np, nt, nquad, na

def MMG2D_Get_solSize(mesh: MMG5_Mesh,sol: MMG5_Sol,typSol):
    ier = lib2d.MMG2D_Get_solSize(ctypes.byref(mesh),ctypes.byref(sol),typEntity.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),np.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),typSol.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier, typEntity, np

def MMG2D_Get_solsAtVerticesSize(mesh: MMG5_Mesh,sol: MMG5_Sol,typSol):
    ier = lib2d.MMG2D_Get_solsAtVerticesSize(ctypes.byref(mesh),ctypes.byref(sol),nsols.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),nentities.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),typSol.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier, nsols, nentities

def MMG2D_Get_vertex(mesh: MMG5_Mesh):
    ier = lib2d.MMG2D_Get_vertex(ctypes.byref(mesh),c0.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),c1.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),ref.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),isCorner.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),isRequired.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier, c0, c1, ref, isCorner, isRequired

def MMG2D_GetByIdx_vertex(mesh: MMG5_Mesh,idx: ctypes.c_int):
    ier = lib2d.MMG2D_GetByIdx_vertex(ctypes.byref(mesh),c0.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),c1.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),ref.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),isCorner.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),isRequired.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),idx)
    return ier, c0, c1, ref, isCorner, isRequired

def MMG2D_Get_vertices(mesh: MMG5_Mesh,vertices,refs,areCorners,areRequired):
    ier = lib2d.MMG2D_Get_vertices(ctypes.byref(mesh),vertices.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),refs.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),areCorners.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),areRequired.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier

def MMG2D_Get_triangle(mesh: MMG5_Mesh):
    ier = lib2d.MMG2D_Get_triangle(ctypes.byref(mesh),v0.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),v1.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),v2.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),ref.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),isRequired.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier, v0, v1, v2, ref, isRequired

def MMG2D_Get_triangles(mesh: MMG5_Mesh,tria,refs,areRequired):
    ier = lib2d.MMG2D_Get_triangles(ctypes.byref(mesh),tria.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),refs.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),areRequired.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier

def MMG2D_Get_quadrilateral(mesh: MMG5_Mesh):
    ier = lib2d.MMG2D_Get_quadrilateral(ctypes.byref(mesh),v0.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),v1.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),v2.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),v3.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),ref.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),isRequired.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier, v0, v1, v2, v3, ref, isRequired

def MMG2D_Get_quadrilaterals(mesh: MMG5_Mesh,quadra,refs,areRequired):
    ier = lib2d.MMG2D_Get_quadrilaterals(ctypes.byref(mesh),quadra.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),refs.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),areRequired.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier

def MMG2D_Get_edge(mesh: MMG5_Mesh):
    ier = lib2d.MMG2D_Get_edge(ctypes.byref(mesh),e0.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),e1.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),ref.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),isRidge.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),isRequired.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier, e0, e1, ref, isRidge, isRequired

def MMG2D_Get_edges(mesh: MMG5_Mesh,edges,refs,areRidges,areRequired):
    ier = lib2d.MMG2D_Get_edges(ctypes.byref(mesh),edges.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),refs.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),areRidges.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),areRequired.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier

def MMG2D_Get_triangleQuality(mesh: MMG5_Mesh,met: MMG5_Sol,k: ctypes.c_int):
    ier = lib2d.MMG2D_Get_triangleQuality(ctypes.byref(mesh),ctypes.byref(met),k)
    return ier

def MMG2D_Get_scalarSol(met: MMG5_Sol,s):
    ier = lib2d.MMG2D_Get_scalarSol(ctypes.byref(met),s.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    return ier

def MMG2D_Get_scalarSols(met: MMG5_Sol,s):
    ier = lib2d.MMG2D_Get_scalarSols(ctypes.byref(met),s.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    return ier

def MMG2D_Get_vectorSol(met: MMG5_Sol):
    ier = lib2d.MMG2D_Get_vectorSol(ctypes.byref(met),vx.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),vy.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    return ier, vx, vy

def MMG2D_Get_vectorSols(met: MMG5_Sol,sols):
    ier = lib2d.MMG2D_Get_vectorSols(ctypes.byref(met),sols.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    return ier

def MMG2D_Get_tensorSol(met: MMG5_Sol):
    ier = lib2d.MMG2D_Get_tensorSol(ctypes.byref(met),m11.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),m12.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),m22.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    return ier, m11, m12, m22

def MMG2D_Get_tensorSols(met: MMG5_Sol,sols):
    ier = lib2d.MMG2D_Get_tensorSols(ctypes.byref(met),sols.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    return ier

def MMG2D_Get_ithSol_inSolsAtVertices(sol: MMG5_Sol,i: ctypes.c_int,s,pos: ctypes.c_int):
    ier = lib2d.MMG2D_Get_ithSol_inSolsAtVertices(ctypes.byref(sol),i,s.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),pos)
    return ier

def MMG2D_Get_ithSols_inSolsAtVertices(sol: MMG5_Sol,i: ctypes.c_int,s):
    ier = lib2d.MMG2D_Get_ithSols_inSolsAtVertices(ctypes.byref(sol),i,s.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    return ier

def MMG2D_Chk_meshData(mesh: MMG5_Mesh,met: MMG5_Sol):
    ier = lib2d.MMG2D_Chk_meshData(ctypes.byref(mesh),ctypes.byref(met))
    return ier

def MMG2D_Free_allSols(mesh: MMG5_Mesh,*sol: MMG5_Sol):
    ier = lib2d.MMG2D_Free_allSols(ctypes.byref(mesh),ctypes.byref(*sol))
    return ier

def MMG2D_loadMesh(mesh: MMG5_Mesh,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_loadMesh(ctypes.byref(mesh),name)
    return ier

def MMG2D_loadVtpMesh(mesh: MMG5_Mesh,met: MMG5_Sol,sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_loadVtpMesh(ctypes.byref(mesh),ctypes.byref(met),ctypes.byref(sol),name)
    return ier

def MMG2D_loadVtpMesh_and_allData(mesh: MMG5_Mesh,*sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_loadVtpMesh_and_allData(ctypes.byref(mesh),ctypes.byref(*sol),name)
    return ier

def MMG2D_loadVtuMesh(mesh: MMG5_Mesh,met: MMG5_Sol,sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_loadVtuMesh(ctypes.byref(mesh),ctypes.byref(met),ctypes.byref(sol),name)
    return ier

def MMG2D_loadVtuMesh_and_allData(mesh: MMG5_Mesh,*sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_loadVtuMesh_and_allData(ctypes.byref(mesh),ctypes.byref(*sol),name)
    return ier

def MMG2D_loadVtkMesh(mesh: MMG5_Mesh,met: MMG5_Sol,sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_loadVtkMesh(ctypes.byref(mesh),ctypes.byref(met),ctypes.byref(sol),name)
    return ier

def MMG2D_loadVtkMesh_and_allData(mesh: MMG5_Mesh,*sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_loadVtkMesh_and_allData(ctypes.byref(mesh),ctypes.byref(*sol),name)
    return ier

def MMG2D_loadMshMesh(mesh: MMG5_Mesh,sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_loadMshMesh(ctypes.byref(mesh),ctypes.byref(sol),name)
    return ier

def MMG2D_loadMshMesh_and_allData(mesh: MMG5_Mesh,*sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_loadMshMesh_and_allData(ctypes.byref(mesh),ctypes.byref(*sol),name)
    return ier

def MMG2D_loadSol(mesh: MMG5_Mesh,sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_loadSol(ctypes.byref(mesh),ctypes.byref(sol),name)
    return ier

def MMG2D_loadGenericMesh(mesh: MMG5_Mesh,met: MMG5_Sol,sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_loadGenericMesh(ctypes.byref(mesh),ctypes.byref(met),ctypes.byref(sol),name)
    return ier

def MMG2D_loadAllSols(mesh: MMG5_Mesh,*sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_loadAllSols(ctypes.byref(mesh),ctypes.byref(*sol),name)
    return ier

def MMG2D_saveMesh(mesh: MMG5_Mesh,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_saveMesh(ctypes.byref(mesh),name)
    return ier

def MMG2D_saveMshMesh(mesh: MMG5_Mesh,sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_saveMshMesh(ctypes.byref(mesh),ctypes.byref(sol),name)
    return ier

def MMG2D_saveMshMesh_and_allData(mesh: MMG5_Mesh,*sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_saveMshMesh_and_allData(ctypes.byref(mesh),ctypes.byref(*sol),name)
    return ier

def MMG2D_saveVtkMesh(mesh: MMG5_Mesh,sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_saveVtkMesh(ctypes.byref(mesh),ctypes.byref(sol),name)
    return ier

def MMG2D_saveVtkMesh_and_allData(mesh: MMG5_Mesh,*sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_saveVtkMesh_and_allData(ctypes.byref(mesh),ctypes.byref(*sol),name)
    return ier

def MMG2D_saveVtuMesh(mesh: MMG5_Mesh,sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_saveVtuMesh(ctypes.byref(mesh),ctypes.byref(sol),name)
    return ier

def MMG2D_saveVtuMesh_and_allData(mesh: MMG5_Mesh,*sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_saveVtuMesh_and_allData(ctypes.byref(mesh),ctypes.byref(*sol),name)
    return ier

def MMG2D_saveVtpMesh(mesh: MMG5_Mesh,sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_saveVtpMesh(ctypes.byref(mesh),ctypes.byref(sol),name)
    return ier

def MMG2D_saveVtpMesh_and_allData(mesh: MMG5_Mesh,*sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_saveVtpMesh_and_allData(ctypes.byref(mesh),ctypes.byref(*sol),name)
    return ier

def MMG2D_saveTetgenMesh(mesh: MMG5_Mesh,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_saveTetgenMesh(ctypes.byref(mesh),name)
    return ier

def MMG2D_saveGenericMesh(mesh: MMG5_Mesh,sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_saveGenericMesh(ctypes.byref(mesh),ctypes.byref(sol),name)
    return ier

def MMG2D_saveSol(mesh: MMG5_Mesh,sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_saveSol(ctypes.byref(mesh),ctypes.byref(sol),name)
    return ier

def MMG2D_saveAllSols(mesh: MMG5_Mesh,*sol: MMG5_Sol,name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_saveAllSols(ctypes.byref(mesh),ctypes.byref(*sol),name)
    return ier

def MMG2D_mmg2dlib(mesh: MMG5_Mesh,sol: MMG5_Sol):
    ier = lib2d.MMG2D_mmg2dlib(ctypes.byref(mesh),ctypes.byref(sol))
    return ier

def MMG2D_mmg2dmesh(mesh: MMG5_Mesh,sol: MMG5_Sol):
    ier = lib2d.MMG2D_mmg2dmesh(ctypes.byref(mesh),ctypes.byref(sol))
    return ier

def MMG2D_mmg2dls(mesh: MMG5_Mesh,sol: MMG5_Sol,met: MMG5_Sol):
    ier = lib2d.MMG2D_mmg2dls(ctypes.byref(mesh),ctypes.byref(sol),ctypes.byref(met))
    return ier

def MMG2D_mmg2dmov(mesh: MMG5_Mesh,met: MMG5_Sol,disp: MMG5_Sol):
    ier = lib2d.MMG2D_mmg2dmov(ctypes.byref(mesh),ctypes.byref(met),ctypes.byref(disp))
    return ier

def MMG2D_defaultValues(mesh: MMG5_Mesh):
    ier = lib2d.MMG2D_defaultValues(ctypes.byref(mesh))
    return ier

def MMG2D_parsop(mesh: MMG5_Mesh,met: MMG5_Sol):
    ier = lib2d.MMG2D_parsop(ctypes.byref(mesh),ctypes.byref(met))
    return ier

def MMG2D_usage(name: str):
    if (isinstance(name,str)):
        name = ctypes.c_char_p(name.encode('utf-8'))
    ier = lib2d.MMG2D_usage(name)
    return ier

def MMG2D_Set_constantSize(mesh: MMG5_Mesh,met: MMG5_Sol):
    ier = lib2d.MMG2D_Set_constantSize(ctypes.byref(mesh),ctypes.byref(met))
    return ier

def MMG2D_setfunc(mesh: MMG5_Mesh,met: MMG5_Sol):
    lib2d.MMG2D_setfunc(ctypes.byref(mesh),ctypes.byref(met))

def MMG2D_Get_numberOfNonBdyEdges(mesh: MMG5_Mesh):
    ier = lib2d.MMG2D_Get_numberOfNonBdyEdges(ctypes.byref(mesh),nb_edges.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier, nb_edges

def MMG2D_Get_nonBdyEdge(mesh: MMG5_Mesh,idx: ctypes.c_int):
    ier = lib2d.MMG2D_Get_nonBdyEdge(ctypes.byref(mesh),e0.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),e1.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),ref.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),idx)
    return ier, e0, e1, ref

def MMG2D_Get_adjaTri(mesh: MMG5_Mesh,kel: ctypes.c_int,listri: ctypes.c_int*3):
    ier = lib2d.MMG2D_Get_adjaTri(ctypes.byref(mesh),kel,listri)
    return ier

def MMG2D_Get_adjaVertices(mesh: MMG5_Mesh,ip: ctypes.c_int,lispoi: ctypes.c_int*MMG2D_LMAX):
    ier = lib2d.MMG2D_Get_adjaVertices(ctypes.byref(mesh),ip,lispoi)
    return ier

def MMG2D_Get_adjaVerticesFast(mesh: MMG5_Mesh,ip: ctypes.c_int,start: ctypes.c_int,lispoi: ctypes.c_int*MMG2D_LMAX):
    ier = lib2d.MMG2D_Get_adjaVerticesFast(ctypes.byref(mesh),ip,start,lispoi)
    return ier

def MMG2D_Get_triFromEdge(mesh: MMG5_Mesh,ked: ctypes.c_int):
    ier = lib2d.MMG2D_Get_triFromEdge(ctypes.byref(mesh),ked,ktri.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),ied.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))
    return ier, ktri, ied

def MMG2D_Get_trisFromEdge(mesh: MMG5_Mesh,ked: ctypes.c_int,ktri: ctypes.c_int*2,ied: ctypes.c_int*2):
    ier = lib2d.MMG2D_Get_trisFromEdge(ctypes.byref(mesh),ked,ktri,ied)
    return ier

def MMG2D_Compute_eigenv(m: ctypes.c_double*3,lambda0: ctypes.c_double*2,vp: ctypes.c_double*2):
    ier = lib2d.MMG2D_Compute_eigenv(m,lambda0,vp)
    return ier

def MMG2D_Reset_verticestags(mesh: MMG5_Mesh):
    lib2d.MMG2D_Reset_verticestags(ctypes.byref(mesh))

def MMG2D_Free_triangles(mesh: MMG5_Mesh):
    lib2d.MMG2D_Free_triangles(ctypes.byref(mesh))

def MMG2D_Free_edges(mesh: MMG5_Mesh):
    lib2d.MMG2D_Free_edges(ctypes.byref(mesh))

def MMG2D_Free_solutions(mesh: MMG5_Mesh,sol: MMG5_Sol):
    lib2d.MMG2D_Free_solutions(ctypes.byref(mesh),ctypes.byref(sol))

