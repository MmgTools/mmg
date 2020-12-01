import numpy as np
import vtk
from vtk.util import numpy_support
import os.path


# get the polydata object
def getPolydata(filename):
    extension = os.path.splitext(filename)[1]

    if extension == ".vtk":
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(filename)
        reader.Update()
        return reader.GetOutput()
    elif extension == ".stl":
        reader = vtk.vtkSTLReader()
        reader.SetFileName(filename)
        reader.Update()
        return reader.GetOutput()
    else:
        print("Bad extension: %s", extension)
        return None

# get the polydata object
def getSTL(filename):
    reader = vtk.vtkSTLReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()

# get points as array (npoints, ndim)
def getPoints(data):
    vtkDataArray = data.GetPoints().GetData()
    return numpy_support.vtk_to_numpy(vtkDataArray)

# get triangle cells connectivity as array (ncells, 3)
def getCellIds(data, filter=None):
    cellIds = []  # global cell ids (not 0..Ncell-1, when filtering)
    vertIds = []  # vert ids of each cell
    idList = vtk.vtkIdList()
    for cellId in range(data.GetNumberOfCells()):
        cellVertIds = []  # list of vertices ids
        cellType = data.GetCellType(cellId)
        if filter == None or cellType == filter:
            data.GetCellPoints(cellId, idList)
            for i in range(0, idList.GetNumberOfIds()):
                vertId = idList.GetId(i)
                cellVertIds.append(vertId)
            vertIds.append(cellVertIds)
            cellIds.append(cellId)
    return np.array(vertIds), cellIds

def getPointsNormals(polyData):

    normalGenerator = vtk.vtkPolyDataNormals()
    normalGenerator.SetInputData(polyData)
    normalGenerator.ComputePointNormalsOn()
    normalGenerator.ComputeCellNormalsOn()
    normalGenerator.SetFeatureAngle(30.0)
    normalGenerator.SetSplitting(0)  # turnoff splitting of cells, this is turned on by default!
    normalGenerator.SetConsistency(1)
    normalGenerator.SetAutoOrientNormals(1)  # is only valid for closed surfaces
    normalGenerator.SetFlipNormals(0)
    normalGenerator.SetNonManifoldTraversal(1)
    normalGenerator.Update()

    pts = numpy_support.vtk_to_numpy(normalGenerator.GetOutput().GetPoints().GetData())
    nms = numpy_support.vtk_to_numpy(normalGenerator.GetOutput().GetPointData().GetNormals())
    return pts, nms  # also points because seems to be different


# write STLdata in file
def writeSTL(STLdata,filename):
    w = vtk.vtkSTLWriter()
    w.SetInputData(STLdata)
    w.SetFileName(filename)
    w.Write()

# write polydata in file
def writePolydata(polydata,filename):
    w = vtk.vtkPolyDataWriter()
    w.SetInputData(polydata)
    w.SetFileName(filename)
    w.Write()


def readUnstructuredGrid(filename, tets_only=False):
    extension = os.path.splitext(filename)[1]
    if extension == ".vtk":
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(filename)
        reader.Update()
        #return reader.GetOutput()
    elif extension == ".vtu":
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(filename)
        reader.Update()
        #return reader.GetOutput()
    else:
        print("Bad extension: %s", extension)
        return None

    if tets_only:  # may renumber points!

        print("reading unstructured grid keeping only tets -> this may result in renumbered points!")

        # select tets only
        ug = reader.GetOutput()
        ids = vtk.vtkIdTypeArray()
        ids.SetNumberOfComponents(1)
        for i in range(ug.GetNumberOfCells()):
            if ug.GetCellType(i) == vtk.VTK_TETRA:
                ids.InsertNextValue(i)

        selectionNode = vtk.vtkSelectionNode()
        selectionNode.SetFieldType(vtk.vtkSelectionNode.CELL)
        selectionNode.SetContentType(vtk.vtkSelectionNode.INDICES)
        selectionNode.SetSelectionList(ids)

        selection = vtk.vtkSelection()
        selection.AddNode(selectionNode)

        extractSelection = vtk.vtkExtractSelection()
        extractSelection.SetInputData(0, ug)
        extractSelection.SetInputData(1, selection)
        extractSelection.Update()

        return extractSelection.GetOutput()
    else:
        return reader.GetOutput()

# write polydata in file
def writeUnstructuredGrid(ug,filename):
    w = vtk.vtkUnstructuredGridWriter()
    w.SetInputData(ug)
    w.SetFileName(filename)
    w.Write()
