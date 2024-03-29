# -*- coding: utf-8 -*-
"""Custom paraview filter to map line orientation to cells."""
import numpy as np
import vtk
from paraview.util.vtkAlgorithm import smdomain, smproperty, smproxy
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid

# Small constant to prevent division by zero
EPS = 1.0e-10


@smproxy.filter(label="Structure Tensor")
@smproperty.input(name="Cells", port_index=0)
@smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)
class StructureTensorFilter(VTKPythonAlgorithmBase):
    """StructureTensorFilter.

    Attributes
    ----------
    _radius : float
        Search radius

    """

    def __init__(self):
        """Set up the filter."""
        VTKPythonAlgorithmBase.__init__(
            self,
            nInputPorts=1,
            nOutputPorts=1,
            outputType="vtkUnstructuredGrid",
        )
        self._radius = 1.0

    def FillInputPortInformation(self, port, info):
        """Require unstructured grids as input types."""
        info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid")
        return 1

    @smproperty.doublevector(name="Radius", default_values=[1.0])
    @smdomain.doublerange()
    def SetRadius(self, radius):
        """Create an input field for the search radius and save it."""
        if self._radius != radius:
            self._radius = radius
            self.Modified()

    def RequestData(self, request, inInfo, outInfo):
        """Process the request submitted after applying the filter."""
        # Access and wrap input source
        source_vtk = vtkUnstructuredGrid.GetData(inInfo[0])
        source = dsa.WrapDataObject(source_vtk)
        N = source.GetNumberOfCells()

        # Compute volume with filter
        filt = vtk.vtkCellSizeFilter()
        filt.SetInputDataObject(source_vtk)
        filt.SetComputeArea(False)
        filt.SetComputeLength(False)
        filt.SetComputeVertexCount(False)
        filt.Update()
        volume_vtk = filt.GetOutput()

        # Compute centers with filter
        filt = vtk.vtkCellCenters()
        filt.SetInputDataObject(volume_vtk)
        filt.Update()
        centers_vtk = filt.GetOutput()
        centers = dsa.WrapDataObject(centers_vtk)

        # Create a vtkKdTreePointLocator for fast computation of neighborhood
        point_locator = vtk.vtkKdTreePointLocator()
        point_locator.SetDataSet(centers_vtk)
        point_locator.BuildLocator()

        # Set up result arrays
        n = np.zeros(N)
        A = np.nan * np.ones((N, 6))

        # Iterate over cells
        for i in range(N):
            # Find neighbors
            pos = centers.Points[i]
            idList = vtk.vtkIdList()
            point_locator.FindPointsWithinRadius(self._radius, pos, idList)

            # Evaluate neighbors
            N0 = idList.GetNumberOfIds()
            ids = [idList.GetId(j) for j in range(N0)]
            points = centers.Points[ids]
            dist = points - pos
            dist_norm = np.linalg.norm(dist, axis=1) + EPS
            dist = dist / dist_norm
            w = centers.PointData["Volume"][ids]

            # Assign number of neighbors
            n[i] = N0

            # Compute structure tensor
            A_struct = np.einsum("k, ki, kj->ij", w, dist, dist) / np.sum(w)

            # Assign structure tensors
            if N0 > 0:
                # Extract symmetric part
                A[i, 0] = A_struct[0, 0]  # XX
                A[i, 1] = A_struct[1, 1]  # YY
                A[i, 2] = A_struct[2, 2]  # ZZ
                A[i, 3] = A_struct[0, 1]  # XY
                A[i, 4] = A_struct[1, 2]  # YZ
                A[i, 5] = A_struct[0, 2]  # XZ

        # Create output
        output = dsa.WrapDataObject(vtkUnstructuredGrid.GetData(outInfo))
        output.ShallowCopy(volume_vtk)
        output.CellData.append(n, "N")
        output.CellData.append(A, "Structure Tensor")

        return 1
