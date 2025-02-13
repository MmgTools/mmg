from writeAPI import *

MMG5_int = "ctypes.c_int"

# Parameters

par = mmgClass("MMG5_Par")
par.addArgs([arg("hmin","ctypes.c_double")])
par.addArgs([arg("hmax","ctypes.c_double")])
par.addArgs([arg("hausd","ctypes.c_double")])
par.addArgs([arg("ref",MMG5_int)])
par.addArgs([arg("elt","ctypes.c_int8")])

MMG5_pPar = "ctypes.POINTER("+par.name+")"

# Point

point = mmgClass("MMG5_Point")
point.addArgs([arg("c","ctypes.c_double * 3")])
point.addArgs([arg("n","ctypes.c_double * 3")])
#point.addArgs([arg("src","ctypes.c_double")])
point.addArgs([arg("ref",MMG5_int)])
point.addArgs([arg("xp",MMG5_int)])
point.addArgs([arg("tmp",MMG5_int)])
point.addArgs([arg("flag",MMG5_int)])
point.addArgs([arg("s",MMG5_int)])
point.addArgs([arg("tag","ctypes.c_uint16")])
point.addArgs([arg("tagdel","ctypes.c_int8")])

MMG5_pPoint = "ctypes.POINTER("+point.name+")"

# xPoint

xpoint = mmgClass("MMG5_xPoint")
xpoint.addArgs(["n1","ctypes.c_double * 3"])
xpoint.addArgs(["n2","ctypes.c_double * 3"])
xpoint.addArgs(["nnor","ctypes.c_int8"])

MMG5_pxPoint = "ctypes.POINTER("+xpoint.name+")"

# mesh

mesh = mmgClass("MMG5_mesh")
mesh.addArgs([arg("memMax","ctypes.c_size_t")])
mesh.addArgs([arg("memCur","ctypes.c_size_t")])
mesh.addArgs([arg("gap","ctypes.c_double")])
mesh.addArgs([arg("ver","ctypes.c_int")])
mesh.addArgs([arg("dim","ctypes.c_int")])
mesh.addArgs([arg("type","ctypes.c_int")])
mesh.addArgs([arg("npi",MMG5_int),arg("nti",MMG5_int),arg("nai",MMG5_int),arg("nei",MMG5_int)])
mesh.addArgs([arg("np",MMG5_int),arg("na",MMG5_int),arg("nt",MMG5_int),arg("ne",MMG5_int)])
mesh.addArgs([arg("npmax",MMG5_int),arg("namax",MMG5_int),arg("ntmax",MMG5_int),arg("nemax",MMG5_int)])
mesh.addArgs([arg("xpmax",MMG5_int),arg("xtmax",MMG5_int)])
mesh.addArgs([arg("nquad",MMG5_int),arg("nprism",MMG5_int)])
mesh.addArgs([arg("nsols","ctypes.c_int")])
mesh.addArgs([arg("nc1",MMG5_int)])
mesh.addArgs([arg("base",MMG5_int)])
mesh.addArgs([arg("mark",MMG5_int)])
mesh.addArgs([arg("xp",MMG5_int),arg("xt",MMG5_int),arg("xpr",MMG5_int)])
mesh.addArgs([arg("npnil",MMG5_int)])
mesh.addArgs([arg("nenil",MMG5_int)])
mesh.addArgs([arg("nanil",MMG5_int)])
mesh.addArgs([arg("adja","ctypes.POINTER("+MMG5_int+")")])
mesh.addArgs([arg("adjt","ctypes.POINTER("+MMG5_int+")")])
mesh.addArgs([arg("adjapr","ctypes.POINTER("+MMG5_int+")")])
mesh.addArgs([arg("adjq","ctypes.POINTER("+MMG5_int+")")])
mesh.addArgs([arg("ipar","ctypes.POINTER(ctypes.c_int)")])
mesh.addArgs([arg("point",MMG5_pPoint)])
mesh.addArgs([arg("xpoint","MMG5_pxPoint")])
mesh.addArgs([arg("tetra","MMG5_pTetra")])
mesh.addArgs([arg("xtetra","MMG5_pxTetra")])
mesh.addArgs([arg("prism","MMG5_pPrism")])
mesh.addArgs([arg("xprism","MMG5_pxPrism")])
mesh.addArgs([arg("tria","MMG5_pTria")])
mesh.addArgs([arg("quadra","MMG5_pQuad")])
mesh.addArgs([arg("edge","MMG5_pEdge")])
mesh.addArgs([arg("htab","MMG5_HGeom")])
mesh.addArgs([arg("info","MMG5_Info")])
mesh.addArgs([arg("namein","ctypes.c_char_p")])
mesh.addArgs([arg("nameout","ctypes.c_char_p")])



