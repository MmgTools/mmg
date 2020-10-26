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
	int* m_tets;
	int m_nverts, m_ntets, m_nfaces;


public:

  double m_hmin = -1.0;
  double m_hgrad = 1.3;
  double m_hausd = 0.01;

  TetMesh(py::array_t<double> verts, py::array_t<int> tets, py::array_t<double>scalars) {
	  //TODO check correct ordering/packing of input arrays

	  /*handle numpy arrays*/
	  py::buffer_info verts_buf = verts.request();
	  py::buffer_info tets_buf = tets.request();
	  py::buffer_info scalars_buf = scalars.request();
	  
	  // check arrays
	  if (verts_buf.size == scalars_buf.size * 3) {
		  throw std::runtime_error("Input shapes must match");
	  }
	  
	  /*init mesh and sol (vertex scalar field here)*/
	  MMG3D_Init_mesh(MMG5_ARG_start,
		  MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet, &mmgSol,
		  MMG5_ARG_end);
	  
	  if (MMG3D_Set_meshSize(mmgMesh, verts_buf.size, tets_buf.size / 4, 0, 0, 0, 0) != 1) {
		  throw std::runtime_error("MMG3D_Set_meshSize failed");
	  }

	  if (MMG3D_Set_solSize(mmgMesh, mmgSol, MMG5_Vertex, verts_buf.size / 4, MMG5_Scalar) != 1) {
		  throw std::runtime_error("MMG3D_Set_solSize failed");
	  }

	  if (MMG3D_Set_vertices(mmgMesh, (double*)verts_buf.ptr, NULL) != 1) {
		  throw std::runtime_error("MMG3D_Set_vertices failed");
	  }
	  
	  if (MMG3D_Set_tetrahedra(mmgMesh, (int*)tets_buf.ptr, NULL) != 1) {
		  throw std::runtime_error("MMG3D_Set_tetra failed");
	  }

	  if (MMG3D_Set_scalarSols(mmgSol, (double*)scalars_buf.ptr) != 1) {
		  throw std::runtime_error("MMG3D_Set_scalarSols failed");
	  }

	  // get the mesh sizes
	  MMG3D_Get_meshSize(mmgMesh, &m_nverts, &m_ntets, NULL, &m_nfaces, NULL, NULL);

	  // set vertices
	  m_verts = new double[m_nverts * 3];
	  MMG3D_Get_vertices(mmgMesh, m_verts, NULL, NULL, NULL);

	  // set tets
	  m_tets = new int[m_ntets * 4];
	  MMG3D_Get_tetrahedra(mmgMesh, m_tets, NULL, NULL);

	}

	~TetMesh() {

		printf("pymmg destructor called\n");

		MMG3D_Free_all(MMG5_ARG_start,
			MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet, &mmgSol,
			MMG5_ARG_end);

		delete(m_verts);
		delete(m_tets);

	}

  void remesh(){

    // set the control parameters
    MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hmin, m_hmin);
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
	MMG3D_Get_tetrahedra(mmgMesh, m_tets, NULL, NULL);

  }

  int getNumberOfTets() {
	  return m_ntets;
  }

  int getNumberOfVerts() {
	  return m_nverts;
  }

  py::array_t<double> getVerts(){
    return py::array_t<double>({m_nverts, 3}, &m_verts[0]);
  }

  py::array_t<int> getTets(){
    return py::array_t<int>({m_ntets, 4}, &m_tets[0]);
  }

};

PYBIND11_MODULE(pymmg, m) {
	py::class_<TetMesh>(m, "TetMesh")
		.def(py::init<py::array_t<double>, py::array_t<int>, py::array_t<double>>())
    .def("remesh", &TetMesh::remesh)
	.def("getNumberOfVerts", &TetMesh::getNumberOfVerts)
	.def("getNumberOfTets", &TetMesh::getNumberOfTets)
    .def("getVerts", &TetMesh::getVerts)
    .def("getTets", &TetMesh::getTets);
}
