import pymmg
import numpy as np
from vtk.util import numpy_support
import vtk
import vtk_helpers

# import a mesh which consists of tetrahedrons and triangles indicating surface
# selection for applying boundary conditions
# TODO also use domain ref/marks for the volume...
ug = vtk_helpers.readUnstructuredGrid("out_102_wBnd.vtk")
N = ug.GetNumberOfPoints()
points = vtk_helpers.getPoints(ug)
tets, tetIds = vtk_helpers.getCellIds(ug, filter=vtk.VTK_TETRA)
faces, faceIds = vtk_helpers.getCellIds(ug, filter=vtk.VTK_TRIANGLE)
refs = numpy_support.vtk_to_numpy(
    ug.GetCellData().GetArray("boundary")
)
face_refs = refs[faceIds]

scalars = np.ones(N)*10.0

# create pymmg mesh object
# remember 0-idx to 1-idx
mesh = pymmg.TetMesh(points.ravel(),
                     tets.ravel()+1,
                     faces.ravel()+1, face_refs.ravel(),
                     scalars.ravel())

# do the remeshing
mesh.hmin = 1.0
mesh.hmax = 10.0
mesh.hausd = 0.05
mesh.hgrad = 1.3
mesh.remesh()

# create new mesh
new_mesh = vtk.vtkUnstructuredGrid()

# set points
vtkPoints = vtk.vtkPoints()
vtkPoints.SetData(numpy_support.numpy_to_vtk(mesh.getVerts(), deep=True))
new_mesh.SetPoints(vtkPoints)

# "refs" array
refsArray = vtk.vtkIntArray()
refsArray.SetName("refs")

# insert tets, remember 1-idx to 0-idx
for tet in mesh.getTets()-1:

    # Create a triangle
    tetra = vtk.vtkTetra()
    tetra.GetPointIds().SetId(0, tet[0])
    tetra.GetPointIds().SetId(1, tet[1])
    tetra.GetPointIds().SetId(2, tet[2])
    tetra.GetPointIds().SetId(3, tet[3])

    # Insert
    new_mesh.InsertNextCell(tetra.GetCellType(), tetra.GetPointIds())

    # Set dummy ref
    refsArray.InsertNextValue(0)

# insert faces, remember 1-idx to 0-idx
for tri, ref in zip(mesh.getFaces()-1, mesh.getFaceRefs()):

    # Create a triangle
    triangle = vtk.vtkTriangle()
    triangle.GetPointIds().SetId(0, tri[0])
    triangle.GetPointIds().SetId(1, tri[1])
    triangle.GetPointIds().SetId(2, tri[2])

    # Insert
    new_mesh.InsertNextCell(triangle.GetCellType(), triangle.GetPointIds())
    refsArray.InsertNextValue(ref)

# insert the refs array
new_mesh.GetCellData().AddArray(refsArray)

vtk_helpers.writeUnstructuredGrid(new_mesh,"new.vtk")
