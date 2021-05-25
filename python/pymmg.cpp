
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <stdio.h>
#include "mmg/mmg3d/libmmg3d.h"

namespace py = pybind11;

class TetMesh {
	//remeshable mmg tetmesh

private:

	MMG5_pMesh mmgMesh = NULL;
	MMG5_pSol mmgSol = NULL;

	double* m_verts;
    int *m_tets, *m_tetRefs, *m_faces, *m_faceRefs;
	int m_nverts, m_ntets, m_nfaces;


public:

  int m_nreg = 0;

  double m_hmax = -1.0;
  double m_hmin = -1.0;
  double m_hgrad = 1.3;
  double m_hausd = 0.01;

  TetMesh(
      py::array_t<double> verts, 
      py::array_t<int> tets,
      py::array_t<int> tetRefs,
      py::array_t<int> faces,
      py::array_t<int> faceRefs,
      py::array_t<double>vertSizes
  ) {
	  //TODO check correct ordering/packing of input arrays

    /*handle numpy arrays*/
    py::buffer_info verts_buf = verts.request();
    py::buffer_info tets_buf = tets.request();
    py::buffer_info tetRefs_buf = tetRefs.request();
    py::buffer_info faces_buf = faces.request();
    py::buffer_info faceRefs_buf = faceRefs.request();
    py::buffer_info vertSizes_buf = vertSizes.request();

    // check arrays
    if (verts_buf.size / 3 != vertSizes_buf.size) {
	    throw std::runtime_error("Input shapes must match");
    }
    if (tets_buf.size / 4 != tetRefs_buf.size) {
            throw std::runtime_error("Input shapes must match");
    }
    if (faces_buf.size / 3 != faceRefs_buf.size) {
        throw std::runtime_error("Input shapes must match");
    }

    /*init mesh and sol (vertex scalar field here)*/
    MMG3D_Init_mesh(MMG5_ARG_start,
	    MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet, &mmgSol,
	    MMG5_ARG_end);

    if (MMG3D_Set_meshSize(mmgMesh, verts_buf.size / 3, tets_buf.size / 4, 0, faces_buf.size / 3, 0, 0) != 1) {
	    throw std::runtime_error("MMG3D_Set_meshSize failed");
    }

    if (MMG3D_Set_solSize(mmgMesh, mmgSol, MMG5_Vertex, vertSizes_buf.size, MMG5_Scalar) != 1) {
	    throw std::runtime_error("MMG3D_Set_solSize failed");
    }

    if (MMG3D_Set_vertices(mmgMesh, (double*)verts_buf.ptr, NULL) != 1) {
	    throw std::runtime_error("MMG3D_Set_vertices failed");
    }

    if (MMG3D_Set_tetrahedra(mmgMesh, (int*)tets_buf.ptr, (int*)tetRefs_buf.ptr) != 1) {
	    throw std::runtime_error("MMG3D_Set_tetra failed");
    }

    if (MMG3D_Set_triangles(mmgMesh, (int*)faces_buf.ptr, (int*)faceRefs_buf.ptr) != 1) {
        throw std::runtime_error("MMG3D_Set_triangles failed");
    }

    if (MMG3D_Set_scalarSols(mmgSol, (double*)vertSizes_buf.ptr) != 1) {
	    throw std::runtime_error("MMG3D_Set_scalarSols failed");
    }

    // get the mesh sizes
    MMG3D_Get_meshSize(mmgMesh, &m_nverts, &m_ntets, NULL, &m_nfaces, NULL, NULL);

    // get vertices
    m_verts = new double[m_nverts * 3];
    MMG3D_Get_vertices(mmgMesh, m_verts, NULL, NULL, NULL);

    // get tets and their refs for keeping track of domains/materials
    m_tets = new int[m_ntets * 4];
    m_tetRefs = new int[m_ntets];
    MMG3D_Get_tetrahedra(mmgMesh, m_tets, m_tetRefs, NULL);

    // get faces (triangles) and their refs for keeping track of bnd-surfaces
    m_faces = new int[m_nfaces * 3];
    m_faceRefs = new int[m_nfaces];
    MMG3D_Get_triangles(mmgMesh, m_faces, m_faceRefs, NULL);

	}

	~TetMesh() {

		MMG3D_Free_all(MMG5_ARG_start,
                   MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet, &mmgSol,
                   MMG5_ARG_end);

		delete(m_verts);
		delete(m_tets);
    delete(m_tetRefs);
    delete(m_faces);
    delete(m_faceRefs);

	}

  void remesh(){

    // set the control parameters
    MMG3D_Set_iparameter(mmgMesh, mmgSol, MMG3D_IPARAM_nreg, m_nreg);
    MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hmin, m_hmin);
    MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hmax, m_hmax);
    MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hausd, m_hausd);
    MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hgrad, m_hgrad);

    // run the mesher
    int ierr = MMG3D_mmg3dlib(mmgMesh, mmgSol);
    if (ierr == MMG5_STRONGFAILURE) {
      throw std::runtime_error("MMG was not able to remesh the mesh.");
    }
    else if (ierr == MMG5_LOWFAILURE){
      printf("MMG completed remeshing with a warning");
    }

    // get the new mesh sizes
    MMG3D_Get_meshSize(mmgMesh, &m_nverts, &m_ntets, NULL, &m_nfaces, NULL, NULL);

    // update vertices
    delete(m_verts);
    m_verts = new double[m_nverts * 3];
    MMG3D_Get_vertices(mmgMesh, m_verts, NULL, NULL, NULL);

    // update tets
    delete(m_tets);
    m_tets = new int[m_ntets * 4];
    m_tetRefs = new int[m_ntets];
    MMG3D_Get_tetrahedra(mmgMesh, m_tets, m_tetRefs, NULL);

    // get faces (triangles) and their refs for keeping track of bnd-surfaces
    m_faces = new int[m_nfaces * 3];
    m_faceRefs = new int[m_nfaces];
    MMG3D_Get_triangles(mmgMesh, m_faces, m_faceRefs, NULL);

  }

  int getNumberOfVerts() {
    return m_nverts;
  }

  int getNumberOfTets() {
      return m_ntets;
  }

  int getNumberOfFaces() {
      return m_nfaces;
  }

  py::array_t<double> getVerts(){
    return py::array_t<double>({m_nverts, 3}, &m_verts[0]);
  }

  py::array_t<int> getTets(){
    return py::array_t<int>({m_ntets, 4}, &m_tets[0]);
  }

  py::array_t<int> getTetRefs() {
      return py::array_t<int>({ m_ntets }, &m_tetRefs[0]);
  }

  py::array_t<int> getFaces() {
      return py::array_t<int>({ m_nfaces, 3 }, &m_faces[0]);
  }


  py::array_t<int> getFaceRefs() {
      return py::array_t<int>({ m_nfaces }, &m_faceRefs[0]);
  }


};

PYBIND11_MODULE(pymmg, m) {
	py::class_<TetMesh>(m, "TetMesh")
		.def(py::init<
            py::array_t<double>, 
            py::array_t<int>,
            py::array_t<int>, 
            py::array_t<int>,
            py::array_t<int>,
            py::array_t<double>
        >())
    .def("remesh", &TetMesh::remesh)
    .def("getNumberOfVerts", &TetMesh::getNumberOfVerts)
    .def("getNumberOfTets", &TetMesh::getNumberOfTets)
	.def("getNumberOfFaces", &TetMesh::getNumberOfFaces)
    .def("getVerts", &TetMesh::getVerts)
    .def("getTets", &TetMesh::getTets)
	.def("getTetRefs", &TetMesh::getTetRefs)
	.def("getFaces", &TetMesh::getFaces)
	.def("getFaceRefs", &TetMesh::getFaceRefs)
	.def_readwrite("nreg", &TetMesh::m_nreg)
    .def_readwrite("hmin", &TetMesh::m_hmin)
    .def_readwrite("hmax", &TetMesh::m_hmax)
    .def_readwrite("hgrad", &TetMesh::m_hgrad)
    .def_readwrite("hausd", &TetMesh::m_hausd);
}
