from mmg_classes import *

def setAPIFunctions(api : pythonAPI):

    mmg3D_init_filenames = mmgFunction("MMG3D_init_filenames","None",[arg("mesh",api.typenames.get("pMesh")),arg("sol",api.typenames.get("pSol"))])
    api.addFunction(mmg3D_init_filenames)